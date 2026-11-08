"""E2 (pre-registered, descriptive): how much CAN the baseline change the propagated prediction?

For each v2 positive (sQTL: lead-intron usage; isoQTL: lead-isoform proportion):
1. Structural sign invariance. For the multiplicative model,
     d_pi_t = pi_t (f_t - F) / F,  F = sum_u pi_u f_u,
   so the sign is the same for every baseline iff f_t is >= (or <=) every other factor with weight.
   For intron usage within cluster K(i), the sign is baseline-invariant iff the factors of transcripts
   containing i all lie on one side of the factors of transcripts containing a competitor
   (transcripts outside the cluster cancel).
2. Achievable range: 1,000 Dirichlet(1) baselines over the gene's annotated transcripts (seed 7):
   fraction of draws with positive sign, and the 5th-95th percentile of the prediction.
3. Ranking sensitivity: Spearman of signed predictions across positives between two independent
   random baselines (500 pairs), compared with fetal vs adult.
4. Metric ceiling: the largest possible |fetal - adult| in direction accuracy equals the number of
   positives whose predicted sign differs between the two baselines, divided by n.

Output: publication_extension/results/E2_root_cause/{per_variant.tsv, summary.json}
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from propagation import intron_membership, transcript_factors  # noqa: E402
from splice_graph import load_lift37  # noqa: E402

OUT = ROOT / "publication_extension" / "results" / "E2_root_cause"
V2 = ROOT / "results" / "experiment_v2"
N_DRAW, SEED, TOL = 1000, 7, 1e-9


def usage(pi, m, comp, k):
    inc = m @ pi
    return inc[k] / inc[comp[k]].sum()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ev = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    sc = pd.read_csv(ROOT / "results" / "spliceai" / "scores.tsv", sep="\t")
    sc["row"] = np.arange(len(sc))
    key = ["chrom", "pos", "ref", "alt", "strand"]
    ev = ev.merge(sc[key + ["row"]], on=key)
    pos = ev[ev.role == "positive"].reset_index(drop=True)
    pred = pd.read_csv(V2 / "predictions.tsv", sep="\t")
    pp = pred[pred.role == "positive"].reset_index(drop=True)
    assert (pp.variant_id.values == pos.variant_id.values).all()
    tr = np.load(ROOT / "results" / "spliceai" / "tracks.npz")
    W = int(tr["window"])
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", set(pos.ensg))

    rows, draws_all = [], []
    for k, r in enumerate(pos.itertuples()):
        ts = sorted(by_gene[r.ensg], key=lambda t: t.enst)
        f = transcript_factors(ts, r.pos, {x: tr[x][r.row] for x in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}, W)
        D = rng.dirichlet(np.ones(len(ts)), N_DRAW)
        if r.arm == "sQTL":
            introns, index, m, comp = intron_membership(ts, "component")
            s, e = map(int, r.target.split("-"))
            ti = index[(s, e)]
            inside = m[ti] > 0
            others = (m[comp[ti]].sum(axis=0) > 0) & ~inside
            fi, fo = f[inside], f[others]
            invariant = bool(len(fo) == 0 or fi.min() >= fo.max() - TOL or fi.max() <= fo.min() + TOL)
            preds = np.array([usage(d * f / (d @ f), m, comp, ti) - usage(d, m, comp, ti) for d in D])
        else:
            t = [x.enst for x in ts].index(r.target)
            fo = np.delete(f, t)
            invariant = bool(f[t] >= fo.max() - TOL or f[t] <= fo.min() + TOL)
            preds = np.array([d[t] * f[t] / (d @ f) - d[t] for d in D])
        all_one = bool(np.allclose(f, 1.0))
        rows.append(dict(arm=r.arm, ensg=r.ensg, variant_id=r.variant_id, z_obs=r.z_obs, n_tx=len(ts),
                         all_factors_one=all_one, sign_invariant_structural=invariant,
                         frac_draws_positive=float((preds > 0).mean()), pred_p05=float(np.quantile(preds, .05)),
                         pred_p95=float(np.quantile(preds, .95)), pred_absmax=float(np.abs(preds).max()),
                         pred_fetal=pp.loc[k, "pred_fetal"], pred_adult=pp.loc[k, "pred_adult"],
                         pred_generic=pp.loc[k, "pred_generic"], site_maxdelta=pp.loc[k, "site_maxdelta"]))
        draws_all.append(preds)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "per_variant.tsv", sep="\t", index=False)

    summ = {}
    for arm in ["sQTL", "isoQTL"]:
        a = d[d.arm == arm]
        idx = np.flatnonzero(d.arm.values == arm)
        M = np.stack([draws_all[i] for i in idx])        # variants x draws
        rand_rho = [spearmanr(M[:, rng.integers(N_DRAW)], M[:, rng.integers(N_DRAW)]).statistic for _ in range(500)]
        flips_fa = int((np.sign(a.pred_fetal) != np.sign(a.pred_adult)).sum())
        flips_fg = int((np.sign(a.pred_fetal) != np.sign(a.pred_generic)).sum())
        inf = a[a.site_maxdelta >= 0.05]
        summ[arm] = {
            "n_positive": len(a),
            "all_factors_one (no annotated site within window changes)": int(a.all_factors_one.sum()),
            "sign_invariant_by_structure": int(a.sign_invariant_structural.sum()),
            "sign_invariant_by_structure_informative": f"{int(inf.sign_invariant_structural.sum())} of {len(inf)}",
            "sign_flippable_in_principle": int((~a.sign_invariant_structural).sum()),
            "random_baselines_frac_variants_with_both_signs": round(float(((a.frac_draws_positive > 0) & (a.frac_draws_positive < 1)).mean()), 3),
            "observed_sign_flips_fetal_vs_adult": flips_fa,
            "observed_sign_flips_fetal_vs_generic": flips_fg,
            "max_possible_abs_diff_dir_acc_fetal_vs_adult": round(flips_fa / len(a), 4),
            "spearman_signed_pred_between_random_baselines_median": round(float(np.median(rand_rho)), 3),
            "spearman_signed_pred_fetal_vs_adult": round(float(spearmanr(a.pred_fetal, a.pred_adult).statistic), 3),
            "median_ratio_absmax_random_to_abs_fetal": round(float(np.median(a.pred_absmax / (a.pred_fetal.abs() + 1e-12))), 1),
        }
    (OUT / "summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()

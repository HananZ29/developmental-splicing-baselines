"""E5 (pre-registered + deviations D1, D2): test the propagation equation against measured exon-inclusion
changes (MFASS, Chong et al. 2019; HEK293T minigene reporter; NOT brain, NOT developmental).

Two-isoform model of a cassette exon in the reporter: the inclusion isoform uses the exon's acceptor
and donor; the skipping isoform uses neither (both share the reporter's flanking sites, which cancel).
  a = r_acc * r_don,  r = (p_alt + eps)/(p_ref + eps), eps = 0.01
  p' = p a / (p a + 1 - p)            (identical to logit(p') = logit(p) + log a)
Predictors (PREREGISTRATION.md E5): P0 raw site change, P1 measured baseline, P2 constant baseline
(train-exon mean), P3 shuffled baseline, P4 p = 0.5, P5 calibrated logit scale k fit on train exons,
P6 null (measured p, permuted a).
Split: 50% of exons train / 50% test (seed 20261104). Metrics on test exons; bootstrap over exons.
Output: publication_extension/results/E5_mfass/{predictions.tsv, metrics.tsv, contrasts.tsv, report.json}
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import expit, logit
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "publication_extension" / "results" / "E5_mfass"
EPS, SEED, N_BOOT = 0.01, 20261104, 2000
KEY = ["chrom", "pos", "ref", "alt", "strand"]
PRED = ["P0_raw", "P1_prop_measured", "P2_prop_constant", "P3_prop_shuffled", "P4_prop_half", "P5_calibrated", "P6_null_permuted_factor"]


def site_values(v, tr, W, col, kind):
    o = (v[col] - v["pos"]).to_numpy()
    ok = np.abs(o) <= W
    rr, ra = np.full(len(v), np.nan), np.full(len(v), np.nan)
    i = np.flatnonzero(ok)
    rows = v["track_row"].to_numpy()[i]
    rr[i] = tr[f"{kind}_ref"][rows, W + o[i]]
    ra[i] = tr[f"{kind}_alt"][rows, W + o[i]]
    return rr, ra


def check_convention(v, tr, W):
    """REF-only check that probabilities peak at the annotated exon sites (first 100 variants)."""
    out = {}
    for col, kind in [("acceptor", "acc"), ("donor", "don")]:
        offs = []
        for r in v.head(100).itertuples():
            c = W + getattr(r, col) - r.pos
            if 3 <= c < 2 * W - 2:
                win = tr[f"{kind}_ref"][r.track_row, c - 3:c + 4]
                offs.append(int(np.argmax(win)) - 3)
        vals, cnt = np.unique(offs, return_counts=True)
        out[col] = {"modal_offset": int(vals[np.argmax(cnt)]), "frac_offset0": round(float(np.mean(np.array(offs) == 0)), 3), "n": len(offs)}
    return out


def prop(p, a):
    return p * a / (p * a + 1 - p) - p


def metrics(y, x):
    big = np.abs(y) >= 0.1
    s = np.sign(x[big])
    return {"spearman": spearmanr(x, y).statistic, "pearson": pearsonr(x, y).statistic,
            "mae": float(np.mean(np.abs(x - y))) if x.std() > 0 else np.nan,
            "dir_acc_abs_dpsi_ge_0.1": float(np.mean(np.where(s == 0, 0.5, s == np.sign(y[big]))))}


def main():
    rng = np.random.default_rng(SEED)
    v = pd.read_csv(D / "mfass_variants.tsv", sep="\t")
    sc = pd.read_csv(D / "spliceai" / "scores.tsv", sep="\t")
    sc["track_row"] = np.arange(len(sc))
    v = v.merge(sc[KEY + ["track_row", "DS_max"]], on=KEY, how="inner")
    t = np.load(D / "spliceai" / "tracks.npz")
    W = int(t["window"])
    tr = {k: t[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    conv = check_convention(v, tr, W)
    for col in ["acceptor", "donor"]:           # pre-registered rule: shift by the modal REF-peak offset
        v[col] = v[col] + conv[col]["modal_offset"]   # offsets are measured in genomic order
    ar, aa = site_values(v, tr, W, "acceptor", "acc")
    dr, da = site_values(v, tr, W, "donor", "don")
    ok = np.isfinite(ar) & np.isfinite(dr)
    v = v[ok].copy()
    ar, aa, dr, da = ar[ok], aa[ok], dr[ok], da[ok]
    v["p_ref_acc"], v["p_ref_don"] = ar, dr
    v["a"] = (aa + EPS) / (ar + EPS) * (da + EPS) / (dr + EPS)
    p = v.nat_v2_index.clip(1e-4, 1 - 1e-4).to_numpy()
    y = v.v2_dpsi.to_numpy()

    exons = v.exon_id.unique()
    train = set(rng.choice(exons, len(exons) // 2, replace=False))
    v["split"] = np.where(v.exon_id.isin(train), "train", "test")
    tr_m = (v.split == "train").to_numpy()

    v["P0_raw"] = (aa - ar) + (da - dr)
    v["P1_prop_measured"] = prop(p, v.a.to_numpy())
    v["P2_prop_constant"] = prop(np.full(len(v), p[tr_m].mean()), v.a.to_numpy())
    ex_p = v.groupby("exon_id").nat_v2_index.first()
    shuf = dict(zip(ex_p.index, rng.permutation(ex_p.values)))
    v["P3_prop_shuffled"] = prop(v.exon_id.map(shuf).clip(1e-4, 1 - 1e-4).to_numpy(), v.a.to_numpy())
    v["P4_prop_half"] = prop(np.full(len(v), 0.5), v.a.to_numpy())
    la = np.log(v.a.to_numpy())

    def loss(k):
        return np.mean((expit(logit(p[tr_m]) + k * la[tr_m]) - p[tr_m] - y[tr_m]) ** 2)
    k = minimize_scalar(loss, bounds=(0, 5), method="bounded").x
    v["P5_calibrated"] = expit(logit(p) + k * la) - p
    v["P6_null_permuted_factor"] = prop(p, rng.permutation(v.a.to_numpy()))
    v.to_csv(D / "predictions.tsv", sep="\t", index=False)

    test = v[v.split == "test"].reset_index(drop=True)
    ex = test.exon_id.to_numpy()
    uex = np.unique(ex)
    idx = {e: np.flatnonzero(ex == e) for e in uex}
    point = {c: metrics(test.v2_dpsi.to_numpy(), test[c].to_numpy()) for c in PRED}
    draws = {c: [] for c in PRED}
    for _ in range(N_BOOT):
        sel = np.concatenate([idx[e] for e in rng.choice(uex, len(uex))])
        yy = test.v2_dpsi.to_numpy()[sel]
        for c in PRED:
            draws[c].append(metrics(yy, test[c].to_numpy()[sel]))
    rows = []
    for c in PRED:
        dfc = pd.DataFrame(draws[c])
        for mname, val in point[c].items():
            rows.append(dict(predictor=c, metric=mname, estimate=val, ci_low=dfc[mname].quantile(.025), ci_high=dfc[mname].quantile(.975)))
    met = pd.DataFrame(rows)
    con = []
    for a, b in [("P1_prop_measured", "P2_prop_constant"), ("P1_prop_measured", "P3_prop_shuffled"),
                 ("P1_prop_measured", "P0_raw"), ("P5_calibrated", "P1_prop_measured"), ("P1_prop_measured", "P6_null_permuted_factor"),
                 ("P1_prop_measured", "P4_prop_half")]:
        da_, db_ = pd.DataFrame(draws[a]), pd.DataFrame(draws[b])
        for mname in ["pearson", "spearman", "mae", "dir_acc_abs_dpsi_ge_0.1"]:
            d = (da_[mname] - db_[mname]).dropna()
            primary = (a, b) in [("P1_prop_measured", "P2_prop_constant"), ("P1_prop_measured", "P3_prop_shuffled")] and mname in ["pearson", "mae"]
            con.append(dict(contrast=f"{a} - {b}", metric=mname, estimate=point[a][mname] - point[b][mname],
                            ci_low=d.quantile(.025), ci_high=d.quantile(.975),
                            p_boot=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean())),
                            family="primary (Bonferroni 0.0125)" if primary else "secondary"))
    con = pd.DataFrame(con)
    met.to_csv(D / "metrics.tsv", sep="\t", index=False, float_format="%.4g")
    con.to_csv(D / "contrasts.tsv", sep="\t", index=False, float_format="%.4g")
    rep = {"convention_check_ref_only": conv, "n_scored": int(len(v)), "n_test_variants": int(len(test)),
           "n_test_exons": int(len(uex)), "calibrated_k_train": round(float(k), 4),
           "corr_measured_baseline_with_dpsi_test": round(float(spearmanr(test.nat_v2_index, test.v2_dpsi).statistic), 3),
           "frac_test_abs_dpsi_ge_0.1": round(float((test.v2_dpsi.abs() >= 0.1).mean()), 3),
           "point": point}
    (D / "report.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))
    print(con.round(4).to_string())


if __name__ == "__main__":
    main()

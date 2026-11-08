"""E12b (pre-registered): do GSE310131 stage-specific baselines change model predictions?

Propagation exactly as v2 (experiments/08 predict_all; eps 0.01, alpha 1e-3, +-500 nt, component clusters),
for Wen sQTL (147), Wen isoQTL (270) and GTEx sQTL (82) positives.
Conditions: Leung fetal/adult, Patowary pooled, gse_t00 / gse_t04 / gse_t30 (FSM, pooled per stage),
gse_t30_L1..L3 (per lineage), gse_t30_fsm_ism (sensitivity).

E12e: sign flips among |pred| > tau (tau = 1e-12, 1e-6, 1e-3), magnitude ratios, rank agreement.
E12f: descriptive benchmark (direction accuracy with ties |pred| <= 1e-12; Spearman with z),
      1-Mb block bootstrap over positives (B = 2000, seed 11). Not confirmatory (E9: ~6% power).
Output: publication_extension/results/E12_gse310131/
"""
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
PE = ROOT / "publication_extension" / "results"
OUT = PE / "E12_gse310131"
GSE = ["gse_t00", "gse_t04", "gse_t30", "gse_t30_L1", "gse_t30_L2", "gse_t30_L3", "gse_t30_fsm_ism"]
CONDS = ["fetal", "adult", "patowary_fetal"] + GSE
PAIRS = [("gse_t30", "gse_t00"), ("gse_t30", "gse_t04"), ("gse_t04", "gse_t00"), ("gse_t30", "fetal"), ("gse_t30", "adult"),
         ("gse_t30", "patowary_fetal"), ("gse_t00", "fetal"), ("fetal", "adult"), ("gse_t30_L1", "gse_t30_L2"), ("gse_t30_L1", "gse_t30_L3"),
         ("gse_t30_L2", "gse_t30_L3"), ("gse_t30", "gse_t30_fsm_ism")]
CONTRASTS = [("gse_t30", "gse_t00"), ("gse_t30", "gse_t04"), ("gse_t30", "fetal"), ("gse_t30", "adult"), ("fetal", "adult")]
TAUS = [1e-12, 1e-6, 1e-3]
TIE = 1e-12


def load_e34():
    spec = importlib.util.spec_from_file_location("e34", ROOT / "publication_extension" / "scripts" / "E34_stage_swap_and_patowary.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def dir_acc(p, z, tau=TIE):
    s = np.where(np.abs(p) <= tau, 0, np.sign(p))
    return float(np.mean(np.where(s == 0, 0.5, (s == np.sign(z)).astype(float))))


def sensitivity(a, b):
    out = {}
    for tau in TAUS:
        both = (np.abs(a) > tau) & (np.abs(b) > tau)
        out[f"tau={tau:g}"] = {"both_above": int(both.sum()), "sign_flips": int((np.sign(a[both]) != np.sign(b[both])).sum()),
                               "a_above": int((np.abs(a) > tau).sum()), "b_above": int((np.abs(b) > tau).sum())}
    nz = (np.abs(a) > 1e-12) & (np.abs(b) > 1e-12)
    r = np.abs(a[nz]) / np.abs(b[nz])
    out["ratio_quantiles"] = {str(q): round(float(np.quantile(r, q)), 4) for q in [.05, .25, .5, .75, .95]} if nz.any() else {}
    out["frac_ratio_outside_[0.5,2]"] = round(float(((r < .5) | (r > 2)).mean()), 3) if nz.any() else None
    out["frac_ratio_outside_[0.1,10]"] = round(float(((r < .1) | (r > 10)).mean()), 3) if nz.any() else None
    out["spearman_pred"] = round(float(spearmanr(a, b).statistic), 4)
    out["spearman_abs_pred"] = round(float(spearmanr(np.abs(a), np.abs(b)).statistic), 4)
    out["top20_overlap_abs"] = int(len(set(np.argsort(-np.abs(a))[:20]) & set(np.argsort(-np.abs(b))[:20])))
    out["median_abs_diff"] = float(np.median(np.abs(a - b)))
    out["n_abs_diff_ge_1e-3"] = int((np.abs(a - b) >= 1e-3).sum())
    out["n_abs_diff_ge_1e-2"] = int((np.abs(a - b) >= 1e-2).sum())
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    e34 = load_e34()
    e34.CONDS = CONDS
    ev8 = e34.load08()
    ev8.CONDITIONS = CONDS
    base = ev8.load_baselines() | e34.patowary_raw()
    bt = pd.read_csv(OUT / "gse_baseline_counts.tsv.gz", sep="\t")
    raw = base | {c: dict(zip(bt.enst, bt[c])) for c in GSE}

    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from splice_graph import load_lift37
    v2v = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    gtex = pd.read_csv(PE / "E3_gtex" / "eligible_gtex.tsv", sep="\t")
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", set(v2v.ensg) | set(gtex.ensg))
    v2v, tr2, W = e34.with_tracks(v2v, ROOT / "results" / "spliceai" / "scores.tsv", ROOT / "results" / "spliceai" / "tracks.npz")
    pos = v2v[v2v.role == "positive"].copy()
    pw = ev8.predict_all(pos, by_gene, tr2, W, raw)
    pw["arm"] = pw.arm.map({"sQTL": "Wen_sQTL", "isoQTL": "Wen_isoQTL"})
    gt, trg, Wg = e34.with_tracks(gtex, PE / "E3_gtex" / "spliceai" / "scores.tsv", PE / "E3_gtex" / "spliceai" / "tracks.npz")
    pg = ev8.predict_all(gt.assign(arm="sQTL"), by_gene, trg, Wg, raw).assign(arm="GTEx_sQTL")
    allp = pd.concat([pw, pg], ignore_index=True)
    allp.to_csv(OUT / "predictions.tsv", sep="\t", index=False)

    # check: fetal/adult predictions reproduce E6 (same code, same baselines)
    e6 = pd.read_csv(PE / "E6_patowary_donors" / "predictions_per_donor.tsv", sep="\t")
    chk = allp.merge(e6[["arm", "variant_id", "ensg", "pred_fetal", "pred_adult"]], on=["arm", "variant_id", "ensg"], suffixes=("", "_e6"))
    rep = {"reproduces_E6_fetal_adult_max_abs_diff": float(max((chk.pred_fetal - chk.pred_fetal_e6).abs().max(), (chk.pred_adult - chk.pred_adult_e6).abs().max())),
           "n_matched_to_E6": len(chk)}

    rng = np.random.default_rng(11)
    rep["sensitivity"], rep["benchmark"], rep["near_zero"] = {}, {}, {}
    for arm in ["Wen_sQTL", "Wen_isoQTL", "GTEx_sQTL"]:
        p = allp[allp.arm == arm].reset_index(drop=True)
        P = {c: p[f"pred_{c}"].to_numpy() for c in CONDS}
        z = p.z_obs.to_numpy()
        rep["near_zero"][arm] = {c: {"n": len(p), "abs<=1e-12": int((np.abs(P[c]) <= 1e-12).sum()), "abs>1e-3": int((np.abs(P[c]) > 1e-3).sum()),
                                     "abs>1e-2": int((np.abs(P[c]) > 1e-2).sum())} for c in CONDS}
        rep["sensitivity"][arm] = {f"{a} vs {b}": sensitivity(P[a], P[b]) for a, b in PAIRS}
        blk = (p.chrom + ":" + (p.pos // 1_000_000).astype(str)).to_numpy()
        ub = np.unique(blk)
        ix = {u: np.flatnonzero(blk == u) for u in ub}
        draws = [np.concatenate([ix[u] for u in rng.choice(ub, len(ub))]) for _ in range(2000)]
        pt = {c: {"dir_acc_tie1e-12": dir_acc(P[c], z), "spearman": float(spearmanr(P[c], z).statistic)} for c in CONDS}
        con = {}
        for a, b in CONTRASTS:
            for met in ["dir_acc_tie1e-12", "spearman"]:
                f = (lambda q, s: dir_acc(q[s], z[s])) if met.startswith("dir") else (lambda q, s: spearmanr(q[s], z[s]).statistic)
                dd = np.array([f(P[a], s) - f(P[b], s) for s in draws])
                con[f"{a} - {b} | {met}"] = {"estimate": pt[a][met] - pt[b][met], "ci": [float(np.quantile(dd, .025)), float(np.quantile(dd, .975))],
                                             "p_boot": float(min(1, 2 * min((dd <= 0).mean(), (dd >= 0).mean())))}
        rep["benchmark"][arm] = {"n_positives": len(p), "n_blocks": len(ub), "point": pt, "contrasts": con}
    (OUT / "summary_E12b.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()

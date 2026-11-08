"""E15 pre-freeze summary (outcome-blind): population counts, prediction agreement between baselines,
F1 stability and the effect of excluding low-spliced-fraction samples.

Reads only predictions.tsv, predictions_gtex_crossover.tsv, baseline_diagnostics.json and the Wen
cluster->gene map. No QTL statistic is read.
Output: publication_extension/results/E15_devstage/prefreeze_report.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "publication_extension/results/E15_devstage"
CL2GENE = ROOT / "data/external/wen_synapse/leafcutter_clusters_to_genes.txt"


def spearman(a, b):
    m = a.notna() & b.notna()
    return float(a[m].rank().corr(b[m].rank())) if m.sum() > 2 else None


def diffs(a, b):
    m = a.notna() & b.notna()
    a, b = a[m].to_numpy(), b[m].to_numpy()
    d = np.abs(a - b)
    big = np.maximum(np.abs(a), np.abs(b))
    small = np.minimum(np.abs(a), np.abs(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = np.where(small > 0, big / small, np.inf)
    return {"n": int(m.sum()),
            "median_abs_diff": float(np.median(d)), "p90_abs_diff": float(np.percentile(d, 90)),
            "share_abs_diff_>0.01": float(np.mean(d > 0.01)), "share_abs_diff_>0.05": float(np.mean(d > 0.05)),
            "share_>2fold": float(np.mean(fold > 2)),
            "share_opposite_sign": float(np.mean(np.sign(a) * np.sign(b) < 0)),
            "max_achievable_delta_rho(1-spearman)": 1 - float(pd.Series(a).rank().corr(pd.Series(b).rank()))}


def main():
    p = pd.read_csv(RES / "predictions.tsv", sep="\t")
    c2g = pd.read_csv(CL2GENE, sep="\t")
    genes = {c: str(g).split(",") for c, g in zip(c2g.clu, c2g.genes)}
    prim = p.dropna(subset=["pred_F1", "pred_A1", "pred_G", "pred_U"])
    rep = {"population": {
        "units(variant-intron pairs)": int(len(prim)),
        "units_dropped_missing_prediction": int(len(p) - len(prim)),
        "unique_rule_P_variants": int(prim.variant_id_P.nunique()),
        "phenotype_introns": int(prim.phenotype_id.nunique()),
        "clusters": int(prim.cluster.nunique()),
        "genes": int(len({g for c in prim.cluster for g in genes.get(c, [])})),
        "signal_groups": int(prim.signal_group.nunique()),
        "resampling_components_P": int(prim.component_P.nunique()),
        "units_rule_P_equals_rule_D": int((prim.variant_id_P == prim.variant_id_D).sum()),
        "units_DSmax_P_>=0.10": int((prim.DS_max_P >= 0.10).sum()),
    }}
    B = ["F1", "A1", "G", "U"]
    rep["spearman_between_baselines"] = {f"{a}-{b}": spearman(prim[f"pred_{a}"], prim[f"pred_{b}"])
                                         for i, a in enumerate(B) for b in B[i + 1:]}
    rep["differences"] = {f"F1-{b}": diffs(prim.pred_F1, prim[f"pred_{b}"]) for b in ("A1", "G", "U")}
    strong = prim[prim.DS_max_P >= 0.10]
    rep["differences_strong_DSmax>=0.10"] = {f"F1-{b}": diffs(strong.pred_F1, strong[f"pred_{b}"]) for b in ("A1", "G")}
    rep["sensitivity_baselines_vs_primary"] = {
        k: spearman(prim[f"pred_{a}"], prim[f"pred_{b}"]) for k, a, b in [
            ("F1 vs F1_sf(exclude spliced fraction<0.10)", "F1", "F1_sf"),
            ("A1 vs A1_sf", "A1", "A1_sf"),
            ("F1 vs F1_all(+HSB155)", "F1", "F1_all"),
            ("F1 vs F1_std(11 regions)", "F1", "F1_std"),
            ("F1 vs F1_nothin", "F1", "F1_nothin"),
            ("G vs G_nonbrain", "G", "G_nonbrain")]}
    rep["F1_sf_minus_A1_sf_vs_F1_minus_A1"] = {
        "spearman_of_differences": spearman(prim.pred_F1 - prim.pred_A1, prim.pred_F1_sf - prim.pred_A1_sf),
        "spearman_F1_sf-A1_sf": spearman(prim.pred_F1_sf, prim.pred_A1_sf)}
    diag = json.loads((RES / "baseline_diagnostics.json").read_text())
    rep["F1_stability"] = diag["F1_stability"]
    rep["baselines"] = {g: {k: v for k, v in d.items() if k != "samples_per_unit"} for g, d in diag["per_baseline"].items()}
    rep["eligibility_waterfall"] = diag["eligibility_waterfall"]
    rep["exclusion_bias"] = diag["exclusion_bias"]
    gx = RES / "predictions_gtex_crossover.tsv"
    if gx.exists():
        g = pd.read_csv(gx, sep="\t").dropna(subset=["pred_F1", "pred_A1"])
        rep["crossover_population"] = {"units": int(len(g)), "clusters": int(g.cluster.nunique()),
                                       "spearman_F1-A1": spearman(g.pred_F1, g.pred_A1)}
    (RES / "prefreeze_report.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in ("population", "spearman_between_baselines", "differences")}, indent=1))


if __name__ == "__main__":
    main()

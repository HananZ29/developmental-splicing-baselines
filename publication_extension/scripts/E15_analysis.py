"""E15 confirmatory analysis (PREREGISTRATION_E15 sections 6-11). READS TRUTH.

Guard: refuses to run unless publication_extension/developmental_stage/PREREG_E15.sha256 exists and
every listed file still has the recorded SHA-256 (i.e. the pre-registration and code are frozen and
unchanged). Gate 0 must have passed (gate0_power.json verdict PASS).

Inputs: predictions.tsv, predictions_gtex_crossover.tsv (from E15_predict.py), Wen SuSiE file (z),
Wen conditional top variants (bwd_slope, column 18), GTEx Brain_Cortex signifpairs (slope).
Output: publication_extension/results/E15_devstage/E15_results.json
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_stats as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "publication_extension/results/E15_devstage"
LOCK = ROOT / "publication_extension/developmental_stage/PREREG_E15.sha256"
SUSIE = ROOT / "data/external/wen_synapse/mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz"
COND = ROOT / "data/external/wen_synapse/conditional_top_variants.txt"
GTEX = ROOT / "data/external/drive_new/Brain_Cortex.v8.sqtl_signifpairs.txt.gz"
B = 10_000
SEED = 20261106
ALPHA = 0.05
MIN_EFFECT = 0.0228   # MEI, decision D3-B: 0.10*sqrt(2*(1-r)), r = 0.97404 (gate0_mei.json)
FAMILY = ("A1", "G", "U")


def check_lock():
    if not LOCK.exists():
        sys.exit("E15 is not frozen (PREREG_E15.sha256 missing). Refusing to read truth.")
    for line in LOCK.read_text().splitlines():
        if not line.strip():
            continue
        h, p = line.split(None, 1)
        if hashlib.sha256((ROOT / p.strip()).read_bytes()).hexdigest() != h:
            sys.exit(f"Frozen file changed: {p}. Refusing to run.")
    g = json.loads((RES / "gate0_power.json").read_text())
    m = json.loads((RES / "gate0_mei.json").read_text())
    if g["gate0"]["max_typeI_F1-A1_N0_N1"] > 0.05 or m["verdict"] != "PASS" or abs(m["MEI"] - MIN_EFFECT) > 1e-9:
        sys.exit("Gate 0 (D3-B) did not pass; E15 stops before reading truth (PREREG section 10).")


def susie_z(keys):
    """keys: set of (phenotype_id, cs_id, variant_id) -> z"""
    out = {}
    with gzip.open(SUSIE, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        iz = h.index("z")
        for line in f:
            p = line.rstrip("\n").split("\t")
            k = (p[0], p[6], p[1])
            if k in keys:
                out[k] = float(p[iz])
    return out


def evaluate(df, truth_col, comp_col, pred_suffix="", others=FAMILY, ref="F1", rng=None):
    cols = [f"pred_{ref}{pred_suffix}"] + [f"pred_{o}{pred_suffix}" for o in others]
    d = df.dropna(subset=cols + [truth_col])
    if len(d) < 10:
        return {"n_units": int(len(d)), "note": "too few units"}
    R = {ref: S.ranks(d[cols[0]].to_numpy())}
    for o, c in zip(others, cols[1:]):
        R[o] = S.ranks(d[c].to_numpy())
    ry = S.ranks(d[truth_col].to_numpy())
    comps = d[comp_col].astype("category").cat.codes.to_numpy()
    W = S.component_weights(comps, B, rng)
    point, boot = S.contrasts(ry, R, ref=ref, others=others, W=W)
    better, worse, p, adj = S.holm(boot, point, ALPHA)
    rho = {b: float(np.corrcoef(ry, R[b])[0, 1]) for b in R}
    rho_ci = {b: [float(x) for x in np.quantile(S.wcorr(W, ry, R[b]), [0.025, 0.975])] for b in R}
    res = {"n_units": int(len(d)), "n_components": int(len(np.unique(comps))), "rho": rho, "rho_95CI": rho_ci,
           "contrasts": {}}
    order = sorted(p, key=p.get)
    holm_level = {k: ALPHA / (len(order) - i) for i, k in enumerate(order)}
    for o in others:
        ci95 = [float(x) for x in np.quantile(boot[o], [0.025, 0.975])]
        lv = holm_level[o]
        ci_holm = [float(x) for x in np.quantile(boot[o], [lv / 2, 1 - lv / 2])]
        if o in better:
            dec = "Improvement" + (" (practically meaningful)" if point[o] >= MIN_EFFECT else "")
        elif o in worse:
            dec = "Worse"
        elif ci95[0] > -MIN_EFFECT and ci95[1] < MIN_EFFECT:
            dec = "No meaningful difference"
        else:
            dec = "Inconclusive"
        res["contrasts"][f"{ref}-{o}"] = {"delta_rho": float(point[o]), "p_boot": float(p[o]), "p_holm": float(adj[o]),
                                         "ci95": ci95, "ci_holm_level": ci_holm, "holm_level": lv,
                                         "decision": dec}
    return res


def main():
    check_lock()
    rng = np.random.default_rng(SEED)
    pr = pd.read_csv(RES / "predictions.tsv", sep="\t")
    keys = set(zip(pr.phenotype_id, pr.cs_id, pr.variant_id_P)) | set(zip(pr.phenotype_id, pr.cs_id, pr.variant_id_D))
    z = susie_z(keys)
    pr["truth_P"] = [z.get(k, np.nan) for k in zip(pr.phenotype_id, pr.cs_id, pr.variant_id_P)]
    pr["truth_D"] = [z.get(k, np.nan) for k in zip(pr.phenotype_id, pr.cs_id, pr.variant_id_D)]
    out = {"units_missing_truth_P": int(pr.truth_P.isna().sum())}

    # ---- primary
    out["primary"] = evaluate(pr, "truth_P", "component_P", rng=rng)
    F1A1 = out["primary"]["contrasts"]["F1-A1"]["decision"]
    F1G = out["primary"]["contrasts"]["F1-G"]["decision"]

    # ---- sensitivity analyses (section 11)
    sens = {}
    sens["rule_D"] = evaluate(pr, "truth_D", "component_D", pred_suffix="_D", rng=rng)
    sens["strong_DSmax>=0.10"] = evaluate(pr[pr.DS_max_P >= 0.10], "truth_P", "component_P", rng=rng)
    sens["weak_DSmax<0.10"] = evaluate(pr[pr.DS_max_P < 0.10], "truth_P", "component_P", rng=rng)
    one = pr.sort_values(["signal_group", "ph_start"]).groupby("signal_group").head(1)
    sens["one_unit_per_signal_group"] = evaluate(one, "truth_P", "component_P", rng=rng)
    for name, f, a, g in [("F1_all", "F1_all", "A1", "G"), ("G_nonbrain", "F1", "A1", "G_nonbrain"),
                          ("standard_regions", "F1_std", "A1_std", "G"), ("spliced_fraction>=0.10", "F1_sf", "A1_sf", "G"),
                          ("no_thinning", "F1_nothin", "A1_nothin", "G_nothin")]:
        d = pr.rename(columns={f"pred_{f}": "pred_xF", f"pred_{a}": "pred_xA", f"pred_{g}": "pred_xG"})
        d = d.assign(pred_xU=pr.pred_U)
        sens[name] = evaluate(d, "truth_P", "component_P", ref="xF", others=("xA", "xG", "xU"), rng=rng)
    for eps in (0.001, 0.1):
        d = pr.rename(columns={f"pred_{b}_eps{eps}": f"pred_e{b}" for b in ("F1", "A1", "G", "U")})
        sens[f"eps={eps}"] = evaluate(d, "truth_P", "component_P", ref="eF1", others=("eA1", "eG", "eU"), rng=rng)
    # secondary truth: conditional bwd_slope (col 18) matched by phenotype and variant position
    cond = pd.read_csv(COND, sep=r"\s+", header=None)
    cond = cond.rename(columns={0: "phe_id", 8: "var_chr", 9: "var_from", 17: "bwd_slope"})
    cmap = {(r.phe_id, "chr" + str(r.var_chr).replace("chr", ""), int(r.var_from)): r.bwd_slope for r in cond.itertuples()}
    pr["truth_cond"] = [cmap.get((ph, c, p), np.nan) for ph, c, p in zip(pr.phenotype_id, pr.chrom, pr.pos_P)]
    sens["secondary_truth_conditional_bwd_slope"] = evaluate(pr, "truth_cond", "component_P", rng=rng)
    # analysis 8: by F1-vs-A1 baseline difference tertile
    q = pr.cluster_max_abs_usage_diff_F1_A1.quantile([1 / 3, 2 / 3]).to_numpy()
    t = np.digitize(pr.cluster_max_abs_usage_diff_F1_A1, q)
    sens["by_baseline_difference_tertile"] = {f"T{k + 1}": evaluate(pr[t == k], "truth_P", "component_P", others=("A1",), rng=rng)
                                              for k in range(3)}
    out["sensitivity"] = sens

    # ---- robustness rule (section 9)
    def sign(res, k="F1-A1"):
        c = res.get("contrasts", {}).get(k)
        return np.sign(c["delta_rho"]) if c else 0
    rev = [k for k in ("rule_D", "one_unit_per_signal_group")
           if sign(sens[k]) != 0 and sign(sens[k]) != sign(out["primary"])]

    # ---- crossover (secondary, descriptive)
    gx = RES / "predictions_gtex_crossover.tsv"
    if gx.exists():
        g = pd.read_csv(gx, sep="\t")
        want = set(zip(g.phenotype_id, g.variant_id_P))
        sl = {}
        with gzip.open(GTEX, "rt") as f:
            h = f.readline().rstrip("\n").split("\t")
            i_s = h.index("slope")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if (p[1], p[0]) in want:
                    sl[(p[1], p[0])] = float(p[i_s])
        g["truth_gtex"] = [sl.get((ph, v), np.nan) * (-1 if fl else 1)
                           for ph, v, fl in zip(g.phenotype_id, g.variant_id_P, g.allele_flipped_P)]
        adult = evaluate(g, "truth_gtex", "component_P", others=("A1",), rng=rng)
        pre = out["primary"]["contrasts"]["F1-A1"]["delta_rho"]
        ad = adult.get("contrasts", {}).get("F1-A1", {}).get("delta_rho")
        out["crossover"] = {"adult_truth": adult, "prenatal_minus_adult_F1A1_delta":
                            (pre - ad) if ad is not None else None,
                            "note": "Different variant sets per arm; descriptive only (PREREG 11.7)."}

    # ---- E15-level conclusion (section 9)
    if F1A1.startswith("Improvement") and F1G.startswith("Improvement"):
        concl = "Fetal baseline improved junction-level agreement with prenatal sQTL association effects over adult and generic baselines."
    elif F1A1.startswith("Improvement"):
        concl = "Fetal baseline improved junction-level agreement over the adult baseline only."
    else:
        concl = f"No fetal improvement declared (F1-A1: {F1A1}; F1-G: {F1G})."
    if rev and (F1A1.startswith("Improvement")):
        concl = "Inconclusive: primary result not robust (sign reverses under " + ", ".join(rev) + ")."
    out["conclusion"] = concl
    out["claim_limits"] = ("Junction-level, association-based. Not a validation of full-length isoform predictions; "
                           "not causal.")
    (RES / "E15_results.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"primary": out["primary"]["contrasts"], "conclusion": concl}, indent=1))


if __name__ == "__main__":
    main()

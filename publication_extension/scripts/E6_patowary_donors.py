"""E6 (pre-registered, exploratory): per-donor fetal baselines from Patowary et al. 2024 Data S2B.

1. Integrity: sheet B dimensions; per-ENST library sums vs J. Dirawi's pooled table.
2. Baseline agreement (smoothed, on evaluation genes): technical (libraries within donor-zone),
   biological (donors), zones, and stage (Leung fetal vs adult) for comparison.
3. Propagation (v2 settings) with patowary_d209 / d334 / d336 on Wen sQTL, Wen isoQTL and GTEx sQTL:
   sign disagreements, metrics per donor, donor - adult contrasts (block bootstrap, B = 2000).

Output: publication_extension/results/E6_patowary_donors/
"""
import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from propagation import smooth_baseline  # noqa: E402
from splice_graph import load_lift37  # noqa: E402

OUT = ROOT / "publication_extension" / "results" / "E6_patowary_donors"
PE = ROOT / "publication_extension" / "results"
DONORS = ["209", "334", "336"]


def load_e34():
    spec = importlib.util.spec_from_file_location("e34", ROOT / "publication_extension" / "scripts" / "E34_stage_swap_and_patowary.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    b = pd.read_excel(ROOT / "data" / "external" / "patowary" / "DataS2.xlsx", sheet_name="B", engine="openpyxl")
    libs = [c for c in b.columns if c[:3] in DONORS]
    integ = {"rows": len(b), "libraries": len(libs), "library_names": libs,
             "reads_per_library": {c: int(b[c].sum()) for c in libs}}
    b["enst"] = b.transcript_id.astype(str).str.split(".").str[0]
    enst = b[b.enst.str.startswith("ENST")].copy()
    integ["enst_rows"] = len(enst)
    gz = [c for c in libs if c.endswith("GZ")]
    cp = [c for c in libs if c.endswith("CP")]
    enst["fetal"], enst["gz"], enst["cp"] = enst[libs].sum(1), enst[gz].sum(1), enst[cp].sum(1)
    j = pd.read_csv(ROOT / "data" / "external" / "drive_new" / "patowary_shared_transcript_counts.tsv.gz", sep="\t")
    j["enst"] = j.transcript.str.split(".").str[0]
    mrg = j.merge(enst[["enst", "fetal", "gz", "cp"]], on="enst", how="left")
    integ["pooled_table_rows"] = len(j)
    integ["pooled_rows_found_in_S2B"] = int(mrg.fetal.notna().sum())
    integ["max_abs_diff_fetal_gz_cp"] = [float((mrg.fetal_reads - mrg.fetal).abs().max()),
                                         float((mrg.gz_reads - mrg.gz).abs().max()), float((mrg.cp_reads - mrg.cp).abs().max())]

    raw = {}
    for d in DONORS:
        raw[f"patowary_d{d}"] = dict(zip(enst.enst, enst[[c for c in libs if c.startswith(d)]].sum(1)))
        for z in ["GZ", "CP"]:
            raw[f"d{d}_{z}"] = dict(zip(enst.enst, enst[[c for c in libs if c.startswith(d) and c.endswith(z)]].sum(1)))
    for c in libs:
        raw[f"lib_{c}"] = dict(zip(enst.enst, enst[c]))

    e34 = load_e34()
    ev8 = e34.load08()
    base = ev8.load_baselines() | e34.patowary_raw()
    raw_all = base | raw

    # --- agreement on evaluation genes
    v2v = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    gtex = pd.read_csv(PE / "E3_gtex" / "eligible_gtex.tsv", sep="\t")
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", set(v2v.ensg) | set(gtex.ensg))
    names = {g: sorted(t.enst for t in ts) for g, ts in by_gene.items()}

    def vec(c):
        return np.concatenate([smooth_baseline(raw_all[c], names[g]) for g in names])

    def r(a, c):
        return float(np.corrcoef(vec(a), vec(c))[0, 1])
    agree = {"technical (4 libraries within donor-zone), median r": float(np.median(
        [r(f"lib_{a}", f"lib_{c}") for d in DONORS for z in ["GZ", "CP"]
         for a, c in itertools.combinations([x for x in libs if x.startswith(d) and x.endswith(z)], 2)])),
        "biological (donor vs donor), r": {f"{a}~{c}": round(r(f"patowary_d{a}", f"patowary_d{c}"), 3) for a, c in itertools.combinations(DONORS, 2)},
        "zone (GZ vs CP within donor), r": {d: round(r(f"d{d}_GZ", f"d{d}_CP"), 3) for d in DONORS},
        "donor vs Leung fetal, r": {d: round(r(f"patowary_d{d}", "fetal"), 3) for d in DONORS},
        "donor vs Leung adult, r": {d: round(r(f"patowary_d{d}", "adult"), 3) for d in DONORS},
        "Leung fetal vs Leung adult, r": round(r("fetal", "adult"), 3),
        "nonzero_weight_fraction": {d: round(float(np.mean([raw_all[f'patowary_d{d}'].get(t, 0) > 0 for g in names for t in names[g]])), 4) for d in DONORS}}

    # --- propagation with per-donor baselines
    conds = ["fetal", "adult", "patowary_fetal"] + [f"patowary_d{d}" for d in DONORS]
    e34.CONDS = conds
    ev8.CONDITIONS = conds
    e34.CONTRASTS = [(f"patowary_d{d}", "adult") for d in DONORS] + [("fetal", "adult")] + \
                    [(f"patowary_d{a}", f"patowary_d{c}") for a, c in itertools.combinations(DONORS, 2)]
    v2v, tr2, W = e34.with_tracks(v2v, ROOT / "results" / "spliceai" / "scores.tsv", ROOT / "results" / "spliceai" / "tracks.npz")
    pos = v2v[v2v.role == "positive"].copy()
    p_wen = ev8.predict_all(pos, by_gene, tr2, W, raw_all)
    p_wen["arm"] = p_wen.arm.map({"sQTL": "Wen_sQTL", "isoQTL": "Wen_isoQTL"})
    gt, trg, Wg = e34.with_tracks(gtex, PE / "E3_gtex" / "spliceai" / "scores.tsv", PE / "E3_gtex" / "spliceai" / "tracks.npz")
    p_gt = ev8.predict_all(gt.assign(arm="sQTL"), by_gene, trg, Wg, raw_all).assign(arm="GTEx_sQTL")
    allp = pd.concat([p_wen, p_gt], ignore_index=True)
    allp.to_csv(OUT / "predictions_per_donor.tsv", sep="\t", index=False)

    rng = np.random.default_rng(11)
    rows, flips = [], {}
    for arm in ["Wen_sQTL", "Wen_isoQTL", "GTEx_sQTL"]:
        p = allp[allp.arm == arm].reset_index(drop=True)
        point = e34.metrics(p)
        B, _ = e34.boot(p, rng)
        rows += e34.summarize(arm, point, B)
        s = {c: np.sign(p[f"pred_{c}"]) for c in conds}
        flips[arm] = {"n": len(p), **{f"{a} vs {c}": int((s[a] != s[c]).sum()) for a, c in itertools.combinations(conds, 2)
                                      if a.startswith("patowary_d") or c.startswith("patowary_d")}}
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "metrics_with_ci.tsv", sep="\t", index=False, float_format="%.4g")
    rep = {"integrity": integ, "agreement": agree, "sign_disagreements": flips}
    (OUT / "report.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))
    print(tab[tab.metric.isin(["dir_acc_all", "spearman_all"])].round(4).to_string())


if __name__ == "__main__":
    main()

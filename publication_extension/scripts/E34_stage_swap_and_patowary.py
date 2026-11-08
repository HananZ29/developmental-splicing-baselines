"""E3 + E4 (pre-registered): stage-swap test with adult GTEx sQTL, and the Patowary fetal baseline.

Propagation is exactly v2 (experiments/08 predict_all, component clusters, eps=0.01, alpha=1e-3).
Conditions: v2's five + patowary_fetal, patowary_gz, patowary_cp (pooled FL reads per ENST, derived
from Patowary et al. 2024 Data S2 by J. Dirawi; GENCODE-matched transcripts only).

Arms:
  Wen_sQTL (fetal, v2 positives + controls), Wen_isoQTL (fetal), GTEx_sQTL (adult cortex, positives only)
Primary (E3): interaction I = (fetal - adult)_Wen_sQTL - (fetal - adult)_GTEx_sQTL for dir_acc_all
and spearman_all; independent 1-Mb block bootstraps per arm, B = 2000; Bonferroni alpha = 0.025.
Exploratory (E4): patowary_* contrasts.

Output: publication_extension/results/E34/
"""
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from propagation import smooth_baseline  # noqa: E402
from splice_graph import load_lift37  # noqa: E402

OUT = ROOT / "publication_extension" / "results" / "E34"
PE = ROOT / "publication_extension" / "results"
N_BOOT, SEED = 2000, 11
CONDS = ["fetal", "adult", "generic", "gtex_cortex", "uniform", "patowary_fetal", "patowary_gz", "patowary_cp"]
CONTRASTS = [("fetal", "adult"), ("fetal", "generic"), ("patowary_fetal", "adult"), ("patowary_fetal", "fetal"),
             ("patowary_gz", "patowary_cp"), ("fetal", "uniform")]
KEY = ["chrom", "pos", "ref", "alt", "strand"]
DSC = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL", "DS_max"]


def load08():
    os.environ.setdefault("EXPERIMENT_OUT", str(OUT / "_unused"))
    spec = importlib.util.spec_from_file_location("ev8", ROOT / "experiments" / "08_propagate_evaluate.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.CONDITIONS = CONDS
    return m


def patowary_raw():
    t = pd.read_csv(ROOT / "data" / "external" / "drive_new" / "patowary_shared_transcript_counts.tsv.gz", sep="\t")
    t["enst"] = t.transcript.str.split(".").str[0]
    return {f"patowary_{k}": dict(zip(t.enst, t[f"{k}_reads"].fillna(0))) for k in ["fetal", "gz", "cp"]}


def with_tracks(v, sc_path, tr_path):
    sc = pd.read_csv(sc_path, sep="\t")
    sc["track_row"] = np.arange(len(sc))
    v = v.merge(sc[KEY + ["track_row"] + DSC], on=KEY, how="inner")
    tr = np.load(tr_path)
    return v, {k: tr[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}, int(tr["window"])


def dir_acc(p, z):
    s = np.sign(p)
    return float(np.mean(np.where(s == 0, 0.5, (s == np.sign(z)).astype(float))))


def metrics(pos):
    inf = pos[pos.site_maxdelta >= 0.05]
    out = {}
    for c in CONDS:
        out[f"{c}|dir_acc_all"] = dir_acc(pos[f"pred_{c}"], pos.z_obs)
        out[f"{c}|spearman_all"] = spearmanr(pos[f"pred_{c}"], pos.z_obs).statistic
        if len(inf) >= 5:
            out[f"{c}|dir_acc_informative"] = dir_acc(inf[f"pred_{c}"], inf.z_obs)
            out[f"{c}|spearman_informative"] = spearmanr(inf[f"pred_{c}"], inf.z_obs).statistic
    return out


def boot(pos, rng):
    block = (pos.chrom + ":" + (pos.pos // 1_000_000).astype(str)).to_numpy()
    ub = np.unique(block)
    idx = {b: np.flatnonzero(block == b) for b in ub}
    rows = [metrics(pos.iloc[np.concatenate([idx[b] for b in rng.choice(ub, len(ub))])]) for _ in range(N_BOOT)]
    return pd.DataFrame(rows), len(ub)


def summarize(arm, point, B):
    rows = []
    for k, v in point.items():
        c, mname = k.split("|")
        b = B[k].dropna()
        rows.append(dict(arm=arm, metric=mname, condition=c, estimate=v, ci_low=b.quantile(.025), ci_high=b.quantile(.975)))
    for a, c in CONTRASTS:
        for mname in ["dir_acc_all", "spearman_all", "dir_acc_informative", "spearman_informative"]:
            if f"{a}|{mname}" not in point:
                continue
            d = (B[f"{a}|{mname}"] - B[f"{c}|{mname}"]).dropna()
            rows.append(dict(arm=arm, metric=mname, condition=f"{a} - {c}", estimate=point[f"{a}|{mname}"] - point[f"{c}|{mname}"],
                             ci_low=d.quantile(.025), ci_high=d.quantile(.975),
                             p_boot=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    ev8 = load08()
    raw = ev8.load_baselines() | patowary_raw()

    # --- E4 coverage/agreement of Patowary vs Leung on the v2 evaluation genes
    v2v = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    gtex = pd.read_csv(PE / "E3_gtex" / "eligible_gtex.tsv", sep="\t")
    genes = set(v2v.ensg) | set(gtex.ensg)
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", genes)
    names = {g: sorted(t.enst for t in ts) for g, ts in by_gene.items()}
    allt = [t for g in names for t in names[g]]
    vec = {c: np.concatenate([smooth_baseline(raw[c], names[g]) for g in names]) for c in CONDS if c != "uniform"}
    cov = {"evaluation_genes": len(names), "evaluation_transcripts": len(allt),
           "nonzero_weight_fraction": {c: round(float(np.mean([raw[c].get(t, 0) > 0 for t in allt])), 4) for c in vec},
           "pearson_smoothed": {f"{a}~{b}": round(float(np.corrcoef(vec[a], vec[b])[0, 1]), 3)
                                for a, b in [("patowary_fetal", "fetal"), ("patowary_fetal", "adult"), ("fetal", "adult"),
                                             ("patowary_fetal", "generic"), ("patowary_gz", "patowary_cp"),
                                             ("patowary_fetal", "gtex_cortex")]},
           "patowary_total_fetal_reads_matched": float(sum(raw["patowary_fetal"].values()))}

    # --- predictions
    v2v, tr2, W = with_tracks(v2v, ROOT / "results" / "spliceai" / "scores.tsv", ROOT / "results" / "spliceai" / "tracks.npz")
    flag = pd.read_csv(ROOT / "results" / "variants" / "controls_flagged.tsv", sep="\t").set_index(["arm", "variant_id", "ensg"])["gene_sig_in_other_arm"]
    flag = flag[~flag.index.duplicated()]
    v2v["strict_control"] = (v2v.role == "control") & ~pd.Series(flag.reindex(list(zip(v2v.arm, v2v.variant_id, v2v.ensg))).fillna(True).to_numpy(dtype=bool))
    p_wen = ev8.predict_all(v2v, by_gene, tr2, W, raw)
    gt, trg, Wg = with_tracks(gtex, PE / "E3_gtex" / "spliceai" / "scores.tsv", PE / "E3_gtex" / "spliceai" / "tracks.npz")
    gt = gt.assign(arm="sQTL")           # use the intron-usage branch of predict_all
    p_gtex = ev8.predict_all(gt, by_gene, trg, Wg, raw).assign(arm="GTEx_sQTL")
    p_wen["arm"] = p_wen.arm.map({"sQTL": "Wen_sQTL", "isoQTL": "Wen_isoQTL"})
    allp = pd.concat([p_wen, p_gtex], ignore_index=True)
    allp.to_csv(OUT / "predictions_all_conditions.tsv", sep="\t", index=False)

    # --- metrics with bootstrap per arm (positives only)
    tables, boots, n_blocks = [], {}, {}
    for arm in ["Wen_sQTL", "Wen_isoQTL", "GTEx_sQTL"]:
        pos = allp[(allp.arm == arm) & (allp.role == "positive")].reset_index(drop=True)
        point = metrics(pos)
        B, nb = boot(pos, rng)
        boots[arm], n_blocks[arm] = B, nb
        tables += summarize(arm, point, B)
    tab = pd.DataFrame(tables)

    # --- primary E3 interaction
    inter = []
    for base_f in ["fetal", "patowary_fetal"]:
        for mname in ["dir_acc_all", "spearman_all"]:
            dw = boots["Wen_sQTL"][f"{base_f}|{mname}"] - boots["Wen_sQTL"][f"adult|{mname}"]
            dg = boots["GTEx_sQTL"][f"{base_f}|{mname}"] - boots["GTEx_sQTL"][f"adult|{mname}"]
            d = (dw - dg).dropna()
            pw = tab[(tab.arm == "Wen_sQTL") & (tab.metric == mname) & (tab.condition == f"{base_f} - adult")].estimate.iloc[0]
            pg = tab[(tab.arm == "GTEx_sQTL") & (tab.metric == mname) & (tab.condition == f"{base_f} - adult")].estimate.iloc[0]
            inter.append(dict(fetal_baseline=base_f, metric=mname, fetal_minus_adult_Wen=pw, fetal_minus_adult_GTEx=pg,
                              interaction=pw - pg, ci_low=d.quantile(.025), ci_high=d.quantile(.975),
                              p_boot=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean())),
                              role="primary (Bonferroni 0.025)" if base_f == "fetal" else "exploratory"))
    inter = pd.DataFrame(inter)

    # --- shared-variant check: same gene + same intron in both cohorts
    w = allp[(allp.arm == "Wen_sQTL") & (allp.role == "positive")]
    g = allp[(allp.arm == "GTEx_sQTL")]
    wv = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t").query("arm=='sQTL' and role=='positive'")
    shared = wv.merge(gtex, on=["ensg", "target"], suffixes=("_wen", "_gtex"))
    shared_info = {"same_gene_same_intron": len(shared), "same_lead_variant": int((shared.pos_wen == shared.pos_gtex).sum()),
                   "z_sign_concordant_fetal_vs_adult_cohort": int((np.sign(shared.z_obs_wen) == np.sign(shared.z_obs_gtex)).sum())}

    tab.to_csv(OUT / "metrics_with_ci.tsv", sep="\t", index=False, float_format="%.4g")
    inter.to_csv(OUT / "stage_swap_interaction.tsv", sep="\t", index=False, float_format="%.4g")
    rep = {"coverage_agreement": cov, "n_positives": {a: int(((allp.arm == a) & (allp.role == "positive")).sum()) for a in boots},
           "n_blocks": n_blocks, "shared_wen_gtex": shared_info,
           "sign_flips_vs_adult": {a: {c: int((np.sign(allp[(allp.arm == a) & (allp.role == 'positive')][f'pred_{c}']) !=
                                               np.sign(allp[(allp.arm == a) & (allp.role == 'positive')]['pred_adult'])).sum())
                                       for c in ["fetal", "patowary_fetal", "generic", "uniform"]} for a in boots}}
    (OUT / "report.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))
    print(inter.round(4).to_string())
    print(tab[tab.condition.str.contains(" - ")].round(4).to_string())


if __name__ == "__main__":
    main()

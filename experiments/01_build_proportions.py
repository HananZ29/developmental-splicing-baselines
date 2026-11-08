"""Build fetal/adult cortex isoform-proportion tables and check whether they are fit for purpose.

Usage: python experiments/01_build_proportions.py
Outputs: results/proportions/*.tsv.gz and results/proportions/diagnostics.json
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from isoform_proportions import build, load_gtf_transcripts, load_sqanti, long_counts  # noqa: E402

RAW = ROOT / "data" / "raw"
GTF = ROOT / "data" / "external" / "gencode" / "gencode.v31.annotation.gtf.gz"
OUT = ROOT / "results" / "proportions"
MIN_GENE_READS = 20  # pooled FL reads per gene per stage for comparisons


def _pair_r(ca, cb, min_reads):
    """Pearson r of within-gene ENST proportions between two count vectors, genes >=min_reads in both."""
    ta = ca.groupby(level="gene_id").transform("sum")
    tb = cb.groupby(level="gene_id").transform("sum")
    known = ca.index.get_level_values("feature") != "NOVEL"
    n_known = pd.Series(known, index=ca.index).groupby(level="gene_id").transform("sum")
    ok = (ta >= min_reads) & (tb >= min_reads) & known & (n_known > 1)
    return round(float(np.corrcoef((ca / ta)[ok], (cb / tb)[ok])[0, 1]), 3), int(ok.sum())


def replicate_agreement(counts, samples, min_reads=MIN_GENE_READS):
    """Pairwise between-sample r, multi-ENST genes with >=min_reads in both samples."""
    return {f"{a}~{b}": _pair_r(counts[a], counts[b], min_reads)[0] for a, b in itertools.combinations(samples, 2)}


def split_half(counts, samples, min_reads=MIN_GENE_READS):
    """Mean r over all 2-vs-3 sample splits within a stage: the noise ceiling for cross-stage r."""
    rs = []
    for half in itertools.combinations(samples, 2):
        rest = [s for s in samples if s not in half]
        rs.append(_pair_r(counts[list(half)].sum(1), counts[rest].sum(1), min_reads)[0])
    return {"mean_r": round(float(np.mean(rs)), 3), "min_r": float(min(rs)), "max_r": float(max(rs))}


def overdispersion(tab, samples, raw_counts, min_reads=20):
    """Ratio of observed between-sample variance of proportions to binomial sampling variance.

    ~1 means replicates differ only by read sampling; >1 means extra biological/technical variation.
    """
    ratios = []
    for (g, f), row in tab.iterrows():
        if f == "NOVEL" or row["gene_reads"] < min_reads * len(samples):
            continue
        n = raw_counts.loc[g][samples].sum()
        k = raw_counts.loc[(g, f)][samples]
        ok = n >= min_reads
        if ok.sum() < 3:
            continue
        p = k[ok].sum() / n[ok].sum()
        if p in (0, 1):
            continue
        obs = ((k[ok] / n[ok] - p) ** 2 * n[ok]).sum() / (ok.sum() - 1)
        ratios.append(obs / (p * (1 - p)))
    r = np.array(ratios)
    return {"n_features": int(r.size), "median_phi": round(float(np.median(r)), 2),
            "frac_phi_gt_2": round(float((r > 2).mean()), 3)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gtf = load_gtf_transcripts(GTF)
    diag, tabs = {}, {}
    for stage in ["Fetal", "Adult"]:
        for mode in ["exclude", "assign"]:
            tab, samples = build(RAW / f"{stage}CTX_classification.tsv.gz", gtf, ism_mode=mode)
            tab.to_csv(OUT / f"{stage.lower()}_props_ism-{mode}.tsv.gz", sep="\t")
            tabs[(stage, mode)] = (tab, samples)
        tab, samples = tabs[(stage, "exclude")]
        c, _ = load_sqanti(RAW / f"{stage}CTX_classification.tsv.gz")
        counts = long_counts(c, samples, gtf, "exclude")
        genes = tab.groupby(level="gene_id")
        d = {
            "samples": samples,
            "reads_per_sample_retained": counts.sum().astype(int).to_dict(),
            "n_genes": int(genes.ngroups),
            "n_genes_ge_min_reads": int((genes["gene_reads"].first() >= MIN_GENE_READS).sum()),
            "n_known_enst": int((tab.index.get_level_values("feature") != "NOVEL").sum()),
            "novel_read_share": round(float(tab.xs("NOVEL", level="feature")["reads"].sum() / tab["reads"].sum()), 4)
            if "NOVEL" in tab.index.get_level_values("feature") else 0.0,
            "replicate_r": replicate_agreement(counts, samples),
            "split_half_r": split_half(counts, samples),
            "overdispersion": overdispersion(tab, samples, counts),
        }
        d["replicate_r_median"] = float(np.nanmedian([v for v in d["replicate_r"].values() if v is not None]))
        diag[stage] = d

    # ISM sensitivity: pooled proportion of known ENST with vs without ISM reads
    for stage in ["Fetal", "Adult"]:
        a = tabs[(stage, "exclude")][0]
        b = tabs[(stage, "assign")][0]
        j = a.join(b, lsuffix="_ex", rsuffix="_as", how="inner")
        j = j[(j.index.get_level_values("feature") != "NOVEL") & (j["gene_reads_ex"] >= MIN_GENE_READS)]
        diag[stage]["ism_sensitivity_r"] = round(float(np.corrcoef(j["prop_pooled_ex"], j["prop_pooled_as"])[0, 1]), 3)
        diag[stage]["ism_sensitivity_median_abs_diff"] = round(float((j["prop_pooled_ex"] - j["prop_pooled_as"]).abs().median()), 4)

    # Fetal vs adult on ENST, genes covered in both stages; ENSTs absent in one stage get proportion 0
    f, a = tabs[("Fetal", "exclude")][0], tabs[("Adult", "exclude")][0]
    fg = f.groupby(level="gene_id")["gene_reads"].first()
    ag = a.groupby(level="gene_id")["gene_reads"].first()
    shared = fg.index[(fg >= MIN_GENE_READS)].intersection(ag.index[ag >= MIN_GENE_READS])
    fa = pd.concat({"fetal": f["prop_pooled"], "adult": a["prop_pooled"]}, axis=1)
    fa = fa[fa.index.get_level_values("gene_id").isin(shared)].fillna(0.0)
    fa["fetal_gene_reads"] = fa.index.get_level_values("gene_id").map(fg)
    fa["adult_gene_reads"] = fa.index.get_level_values("gene_id").map(ag)
    fa["delta_adult_minus_fetal"] = fa["adult"] - fa["fetal"]
    fa.to_csv(OUT / "fetal_vs_adult_pooled.tsv.gz", sep="\t")
    known = fa[fa.index.get_level_values("feature") != "NOVEL"]
    n_iso = known.groupby(level="gene_id").size()
    multi = known[known.index.get_level_values("gene_id").isin(n_iso[n_iso > 1].index)]
    diag["fetal_vs_adult"] = {
        "min_gene_reads_each_stage": MIN_GENE_READS,
        "shared_genes": int(len(shared)),
        "shared_genes_multi_enst": int((n_iso > 1).sum()),
        "enst_rows": int(len(known)),
        "enst_only_fetal_detected": int(((known["fetal"] > 0) & (known["adult"] == 0)).sum()),
        "enst_only_adult_detected": int(((known["adult"] > 0) & (known["fetal"] == 0)).sum()),
        "pearson_r_multi_isoform": round(float(np.corrcoef(multi["fetal"], multi["adult"])[0, 1]), 3),
        "frac_abs_delta_gt_0.2": round(float((multi["delta_adult_minus_fetal"].abs() > 0.2).mean()), 3),
    }
    (OUT / "diagnostics.json").write_text(json.dumps(diag, indent=2))
    print(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()

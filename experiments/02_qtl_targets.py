"""Harmonize isoQTL / sQTL permutation results to unversioned Ensembl IDs and intersect with the
fetal/adult isoform-proportion tables.

Usage: python experiments/02_qtl_targets.py
Outputs: results/qtl/isoqtl_targets.tsv.gz, results/qtl/sqtl_targets.tsv.gz, results/qtl/summary.json

The `slope` column is a QTLtools regression slope on a per-feature standardized phenotype (the implied
residual SD is ~0.045 for every feature), so it is reported as a signed z (slope/slope_se) and is not
interpreted as a change in isoform proportion. Allele orientation (which allele the sign refers to)
requires the hg19 .bim file that is not yet available locally.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROPS = ROOT / "results" / "proportions"
OUT = ROOT / "results" / "qtl"
FDR = 0.05


def load_isoqtl():
    q = pd.read_csv(RAW / "devBrain_isoQTL_ALL_grouped_permutation_70HCP_FDR0.05 (1).txt.gz", sep="\t")
    feat = q["gene_id"].str.split(":", expand=True)
    q["transcript_id"] = feat[0].str.split(".").str[0]
    q["gene"] = feat[1].str.split(".").str[0]
    return q


def load_sqtl():
    q = pd.read_csv(RAW / "sQTL.mixed.40hcp.group.perm.genes.txt", sep="\t")
    f = q["gene_id"].str.split(":", expand=True)
    q["chrom_hg19"], q["intron_start"], q["intron_end"], q["cluster"] = f[0], f[1].astype(int), f[2].astype(int), f[3]
    q["gene"] = f[4].str.split(".").str[0]
    q["gene_name"] = f[5]
    return q


def add_stats(q):
    q["z"] = q["slope"] / q["slope_se"]
    q["n_individuals"] = (q["ma_count"] / (2 * q["maf"])).round().astype(int)
    q["significant"] = q["qval"] < FDR
    return q


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fa = pd.read_csv(PROPS / "fetal_vs_adult_pooled.tsv.gz", sep="\t")
    fa = fa[fa["feature"] != "NOVEL"].assign(transcript_id=lambda d: d["feature"].str.split(".").str[0])
    fetal = pd.read_csv(PROPS / "fetal_props_ism-exclude.tsv.gz", sep="\t")
    fetal_tx = set(fetal.loc[fetal["feature"] != "NOVEL", "feature"].str.split(".").str[0])
    prop_genes = set(fa["gene_id"])

    iso = add_stats(load_isoqtl())
    iso = iso.merge(fa[["transcript_id", "fetal", "adult", "delta_adult_minus_fetal",
                        "fetal_gene_reads", "adult_gene_reads"]], on="transcript_id", how="left")
    iso.to_csv(OUT / "isoqtl_targets.tsv.gz", sep="\t", index=False)

    sq = add_stats(load_sqtl())
    sq["gene_in_proportion_table"] = sq["gene"].isin(prop_genes)
    sq.to_csv(OUT / "sqtl_targets.tsv.gz", sep="\t", index=False)

    sig = iso[iso["significant"]]
    summary = {
        "isoQTL": {
            "rows": len(iso), "genes": int(iso["gene"].nunique()), "significant_fdr05": int(len(sig)),
            "n_individuals": sorted(iso["n_individuals"].unique().tolist()),
            "sig_lead_isoform_detected_in_fetal_sqanti": int(sig["transcript_id"].isin(fetal_tx).sum()),
            "sig_lead_isoform_in_fetal_vs_adult_table": int(sig["fetal"].notna().sum()),
            "sig_in_table_with_abs_delta_gt_0.1": int((sig["delta_adult_minus_fetal"].abs() > 0.1).sum()),
            "abs_z_median_sig": round(float(sig["z"].abs().median()), 2),
        },
        "sQTL": {
            "rows": len(sq), "genes": int(sq["gene"].nunique()), "significant_fdr05": int(sq["significant"].sum()),
            "n_individuals": sorted(sq["n_individuals"].unique().tolist()),
            "sig_genes_in_fetal_vs_adult_table": int(sq.loc[sq["significant"], "gene"][sq["gene_in_proportion_table"]].nunique()),
        },
        "isoQTL_sQTL_shared_sig_genes": int(len(set(sig["gene"]) & set(sq.loc[sq["significant"], "gene"]))),
        "isoQTL_sQTL_same_lead_variant_sig": int(len(
            set(zip(sig["gene"], sig["variant_id"])) & set(zip(sq.loc[sq["significant"], "gene"], sq.loc[sq["significant"], "variant_id"])))),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

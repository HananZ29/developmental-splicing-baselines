"""Select eligible evaluation variants and distance-matched controls, independently of any prediction.

Gene universe G: autosomal genes with >=20 pooled FL reads in BOTH fetal and adult SQANTI and >=2
GENCODE v31lift37 transcripts (so all three baselines can be defined on the same isoform set).

Positives (lead association variants, NOT fine-mapped causal variants):
- sQTL arm:  FDR<0.05 sGene in G; lead intron annotated in the gene's v31lift37 transcripts and with
             >=1 annotated competitor intron sharing its donor or acceptor; lead variant an SNV
             within WIN nt of an annotated splice site of the gene.
- isoQTL arm: FDR<0.05 isoGene in G; lead isoform present in v31lift37; same variant criteria.

Controls (presumed-null, not confirmed negatives): SNVs from the mapping .bim lying within WIN nt of
an annotated splice site of a gene in G whose QTL q-value is >= 0.2 in the same arm (tested, no
evidence). For each positive, K controls are drawn with the same distance-to-nearest-splice-site bin.

Alleles: hg19 REF = A2 for 100% of .bim SNVs (results/build_checks.json), so ALT = A1 = the QTL effect
allele; observed z is therefore per ALT allele.

Usage: python experiments/06_select_variants.py
Outputs: results/variants/eligible_variants.tsv, attrition.json
"""
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from splice_graph import load_lift37, site_table  # noqa: E402

RAW, EXT = ROOT / "data" / "raw", ROOT / "data" / "external"
OUT = ROOT / "results" / "variants"
WIN, K, MIN_READS, NULL_Q = 500, 4, 20, 0.2
BINS = [-1, 2, 10, 50, 200, WIN]
SEED = 7


def gene_universe():
    fa = pd.read_csv(ROOT / "results" / "proportions" / "fetal_vs_adult_pooled.tsv.gz", sep="\t")
    ok = fa[(fa["fetal_gene_reads"] >= MIN_READS) & (fa["adult_gene_reads"] >= MIN_READS)]
    return set(ok["gene_id"])


def nearest_site(sites_by_gene, ensg, pos):
    s = sites_by_gene.get(ensg)
    if s is None:
        return np.inf
    i = np.searchsorted(s, pos)
    d = [abs(pos - s[j]) for j in (i - 1, i) if 0 <= j < len(s)]
    return min(d) if d else np.inf


def intron_competitors(transcripts):
    """{intron: n_annotated_introns sharing its start or end (incl. itself)}"""
    introns = {i for t in transcripts for i in t.introns}
    by_start, by_end = {}, {}
    for s, e in introns:
        by_start.setdefault(s, set()).add((s, e))
        by_end.setdefault(e, set()).add((s, e))
    return {i: len(by_start[i[0]] | by_end[i[1]]) for i in introns}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    G0 = gene_universe()
    by_gene = load_lift37(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", G0)
    G = {g for g, ts in by_gene.items() if len(ts) >= 2 and ts[0].chrom.replace("chr", "").isdigit()}
    sites = site_table({g: by_gene[g] for g in G})
    sites_by_gene = {g: np.sort(d["pos"].to_numpy()) for g, d in sites.groupby("ensg")}
    gene_info = {g: (by_gene[g][0].chrom, by_gene[g][0].strand) for g in G}
    att = {"genes_fetal_adult_ge20": len(G0), "genes_G_autosomal_ge2_tx": len(G)}

    bim = pd.read_csv(zipfile.ZipFile(EXT / "devbrain" / "test.bim.zip").open("test.bim"), sep="\t", header=None,
                      names=["chr", "id", "cm", "pos", "A1", "A2"], dtype={"chr": str})
    bim = bim[(bim["A1"].str.len() == 1) & (bim["A2"].str.len() == 1)]
    bim_by_id = bim.drop_duplicates("id").set_index("id")

    rows = []

    def add_positive(arm, ensg, vid, q, z, target, extra):
        if vid not in bim_by_id.index:
            return "not_snv"
        b = bim_by_id.loc[vid]
        chrom, strand = gene_info[ensg]
        if "chr" + str(b["chr"]) != chrom:
            return "chrom_mismatch"
        d = nearest_site(sites_by_gene, ensg, int(b["pos"]))
        if d > WIN:
            return "far_from_splice_site"
        rows.append(dict(arm=arm, role="positive", ensg=ensg, variant_id=vid, chrom=chrom, pos=int(b["pos"]),
                         ref=b["A2"], alt=b["A1"], strand=strand, dist_site=int(d), qval=q, z_obs=z,
                         target=target, **extra))
        return "eligible"

    # sQTL arm
    sq = pd.read_csv(RAW / "sQTL.mixed.40hcp.group.perm.genes.txt", sep="\t")
    p = sq["gene_id"].str.split(":", expand=True)
    sq["ensg"] = p[4].str.split(".").str[0]
    sq["intron"] = list(zip(p[1].astype(int) + 1, p[2].astype(int) - 1))
    sq["z"] = sq["slope"] / sq["slope_se"]
    comp = {g: intron_competitors(by_gene[g]) for g in G}
    a = {"sig_rows": int((sq["qval"] < 0.05).sum())}
    sig = sq[(sq["qval"] < 0.05) & sq["ensg"].isin(G)]
    a["sig_in_G"] = len(sig)
    reasons = {}
    for _, r in sig.iterrows():
        nc = comp[r["ensg"]].get(r["intron"], 0)
        why = "intron_not_annotated" if nc == 0 else "no_competitor_intron" if nc < 2 else \
            add_positive("sQTL", r["ensg"], r["variant_id"], r["qval"], r["z"], f"{r['intron'][0]}-{r['intron'][1]}",
                         {"n_competitors": nc})
        reasons[why] = reasons.get(why, 0) + 1
    a["filter_outcomes"] = reasons
    att["sQTL"] = a

    # isoQTL arm
    iso = pd.read_csv(RAW / "devBrain_isoQTL_ALL_grouped_permutation_70HCP_FDR0.05 (1).txt.gz", sep="\t")
    f = iso["gene_id"].str.split(":", expand=True)
    iso["enst"], iso["ensg"] = f[0].str.split(".").str[0], f[1].str.split(".").str[0]
    iso["z"] = iso["slope"] / iso["slope_se"]
    tx_in = {t.enst for g in G for t in by_gene[g]}
    a = {"sig_rows": int((iso["qval"] < 0.05).sum())}
    sig = iso[(iso["qval"] < 0.05) & iso["ensg"].isin(G)]
    a["sig_in_G"] = len(sig)
    reasons = {}
    for _, r in sig.iterrows():
        why = "lead_isoform_not_in_v31lift37" if r["enst"] not in tx_in else \
            add_positive("isoQTL", r["ensg"], r["variant_id"], r["qval"], r["z"], r["enst"], {})
        reasons[why] = reasons.get(why, 0) + 1
    a["filter_outcomes"] = reasons
    att["isoQTL"] = a

    pos = pd.DataFrame(rows)
    pos["dist_bin"] = pd.cut(pos["dist_site"], BINS, labels=False)

    # controls: SNVs near splice sites of tested genes with no QTL evidence in that arm
    null_genes = {"sQTL": set(sq.loc[sq["qval"] >= NULL_Q, "ensg"]) & G - set(sq.loc[sq["qval"] < 0.05, "ensg"]),
                  "isoQTL": set(iso.loc[iso["qval"] >= NULL_Q, "ensg"]) & G - set(iso.loc[iso["qval"] < 0.05, "ensg"])}
    lead_ids = set(sq["variant_id"]) | set(iso["variant_id"])
    bim_chr = {c: d for c, d in bim.groupby("chr")}
    ctrl_rows = []
    for arm, genes in null_genes.items():
        pool = []
        for g in genes:
            chrom, strand = gene_info[g]
            s = sites_by_gene[g]
            d = bim_chr.get(chrom.replace("chr", ""))
            if d is None:
                continue
            d = d[(d["pos"] >= s[0] - WIN) & (d["pos"] <= s[-1] + WIN) & ~d["id"].isin(lead_ids)]
            for _, b in d.iterrows():
                dist = nearest_site(sites_by_gene, g, int(b["pos"]))
                if dist <= WIN:
                    pool.append(dict(arm=arm, role="control", ensg=g, variant_id=b["id"], chrom=chrom,
                                     pos=int(b["pos"]), ref=b["A2"], alt=b["A1"], strand=strand, dist_site=int(dist)))
        pool = pd.DataFrame(pool)
        pool["dist_bin"] = pd.cut(pool["dist_site"], BINS, labels=False)
        att[arm]["null_genes"] = len(genes)
        att[arm]["control_pool"] = len(pool)
        need = pos[pos["arm"] == arm]["dist_bin"].value_counts()
        for b, n in need.items():
            cand = pool[pool["dist_bin"] == b]
            take = min(len(cand), K * n)
            ctrl_rows.append(cand.iloc[rng.choice(len(cand), take, replace=False)])
            att[arm].setdefault("controls_by_bin", {})[int(b)] = {"positives": int(n), "controls": int(take)}

    out = pd.concat([pos] + ctrl_rows, ignore_index=True)
    out.to_csv(OUT / "eligible_variants.tsv", sep="\t", index=False)
    for arm in ["sQTL", "isoQTL"]:
        d = out[out["arm"] == arm]
        att[arm]["final"] = {"positives": int((d["role"] == "positive").sum()),
                             "controls": int((d["role"] == "control").sum()),
                             "genes_positive": int(d.loc[d["role"] == "positive", "ensg"].nunique()),
                             "genes_control": int(d.loc[d["role"] == "control", "ensg"].nunique())}
    att["unique_variant_strand_to_score"] = int(out[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().shape[0])
    (OUT / "attrition.json").write_text(json.dumps(att, indent=2, default=int))
    print(json.dumps(att, indent=2, default=int))


if __name__ == "__main__":
    main()

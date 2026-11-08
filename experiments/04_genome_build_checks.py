"""Establish the genome build of every coordinate-bearing input by sequence evidence, not ID overlap.

Tests against the UCSC hg19 FASTA (MD5-verified):
- variants (.bim, QTL lead variants): the hg19 base at `pos` must equal A1 or A2 for nearly all SNVs.
- GTF annotations: introns of multi-exon transcripts must start with GT and end with AG on the
  transcript strand at ~98-99% if the coordinates are hg19, and near chance if they are not.
- sQTL leafcutter introns: same GT-AG test; also identifies the coordinate convention.
- SQANTI classification: CDS_genomic_start must fall inside the associated transcript's span.

Usage: python experiments/04_genome_build_checks.py
Output: results/build_checks.json
"""
import gzip
import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external"
RAW = ROOT / "data" / "raw"
RNG = np.random.default_rng(0)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_gtf_exons(path, max_tx=None):
    tx = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
    return tx


def gtag_rate(fa, introns):
    """introns: iterable of (chrom, first_intron_base, last_intron_base, strand), 1-based inclusive."""
    ok = n = 0
    for c, s, e, strand in introns:
        if c not in fa or e - s < 4:
            continue
        seq = fa[c][s - 1:e].seq.upper()
        if strand == "-":
            seq = revcomp(seq)
        n += 1
        ok += seq[:2] == "GT" and seq[-2:] == "AG"
    return round(ok / n, 4) if n else None, n


def gtf_check(fa, path, k=4000):
    tx = read_gtf_exons(path)
    keys = [t for t, v in tx.items() if len(v[2]) > 1]
    pick = RNG.choice(len(keys), size=min(k, len(keys)), replace=False)
    introns = []
    for i in pick:
        c, strand, ex = tx[keys[i]]
        ex = sorted(ex)
        introns += [(c, a[1] + 1, b[0] - 1, strand) for a, b in zip(ex, ex[1:])]
    rate, n = gtag_rate(fa, introns)
    return {"gt_ag_rate_hg19": rate, "introns_tested": n}


def variant_check(fa, df, label):
    snv = df[(df["A1"].str.len() == 1) & (df["A2"].str.len() == 1)]
    snv = snv.sample(min(len(snv), 200000), random_state=0)
    ref = [fa[f"chr{c}"][p - 1:p].seq.upper() if f"chr{c}" in fa else "N" for c, p in zip(snv["chr"], snv["pos"])]
    ref = pd.Series(ref, index=snv.index)
    in_pair = ((ref == snv["A1"]) | (ref == snv["A2"])).mean()
    return {"snvs_tested": int(len(snv)), "hg19_ref_in_A1_A2": round(float(in_pair), 5),
            "hg19_ref_equals_A2": round(float((ref == snv["A2"]).mean()), 5),
            "hg19_ref_equals_A1": round(float((ref == snv["A1"]).mean()), 5), "set": label}


def main():
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    out = {"hg19_chr1_length": len(fa["chr1"])}
    bim = pd.read_csv(zipfile.ZipFile(EXT / "devbrain" / "test.bim.zip").open("test.bim"), sep="\t", header=None,
                      names=["chr", "id", "cm", "pos", "A1", "A2"], dtype={"chr": str})
    out["bim_all"] = variant_check(fa, bim, "test.bim random SNVs")
    iso = pd.read_csv(RAW / "devBrain_isoQTL_ALL_grouped_permutation_70HCP_FDR0.05 (1).txt.gz", sep="\t")
    sq = pd.read_csv(RAW / "sQTL.mixed.40hcp.group.perm.genes.txt", sep="\t")
    lead = bim[bim["id"].isin(set(iso["variant_id"]) | set(sq["variant_id"]))]
    out["bim_qtl_leads"] = variant_check(fa, lead, "isoQTL+sQTL lead variants")

    for v in ["v24lift37", "v31lift37", "v33lift37", "v26", "v31"]:
        out[f"gencode.{v}"] = gtf_check(fa, EXT / "gencode" / f"gencode.{v}.annotation.gtf.gz")

    p = sq["gene_id"].str.split(":", expand=True)
    strand_by_gene = {}
    with gzip.open(EXT / "gencode" / "gencode.v24lift37.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if "\tgene\t" in line:
                f = line.split("\t")
                strand_by_gene[re.search(r'gene_id "([^".]+)', f[8]).group(1)] = f[6]
    strand = p[4].str.split(".").str[0].map(strand_by_gene)
    keep = strand.notna()
    s, e = p[1].astype(int), p[2].astype(int)
    for name, (a, b) in {"start..end": (s, e), "start+1..end-1": (s + 1, e - 1), "start..end-1": (s, e - 1),
                         "start+1..end": (s + 1, e)}.items():
        rate, n = gtag_rate(fa, zip("chr" + p[0][keep], a[keep], b[keep], strand[keep]))
        out[f"sqtl_introns_hg19[{name}]"] = {"gt_ag_rate": rate, "n": n}

    # SQANTI: ORF genomic start should sit inside the associated transcript's span in its native build
    spans = {}
    for v in ["v31", "v31lift37"]:
        tx = read_gtf_exons(EXT / "gencode" / f"gencode.{v}.annotation.gtf.gz")
        spans[v] = {t.split("_")[0]: (c, min(x[0] for x in ex), max(x[1] for x in ex)) for t, (c, _, ex) in tx.items()}
    for stage in ["Fetal", "Adult"]:
        c = pd.read_csv(RAW / f"{stage}CTX_classification.tsv.gz", sep="\t", low_memory=False)
        c = c[(c["structural_category"] == "full-splice_match") & c["CDS_genomic_start"].notna()]
        res = {}
        for v, sp in spans.items():
            inside = [(t in sp and sp[t][0] == ch and sp[t][1] <= g <= sp[t][2])
                      for t, ch, g in zip(c["associated_transcript"], c["chrom"], c["CDS_genomic_start"])]
            res[f"cds_start_inside_tx_{v}"] = round(float(np.mean(inside)), 4)
        out[f"sqanti_{stage}"] = res | {"n_fsm_coding": int(len(c))}

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "build_checks.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

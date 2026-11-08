"""Per-gene isoform proportions from SQANTI-classified PacBio full-length (FL) read counts.

Decisions (see docs/analysis_notes.md):
- Isoforms are keyed by GENCODE v31 ENST (the release the SQANTI files were classified against),
  never by PB id, whose numbering is independent per dataset.
- Numerator: FL reads of all full-splice-match (FSM) PB isoforms assigned to the ENST, summed.
  FSM PBs sharing an ENST differ only in TSS/TTS ends, so they are the same splice isoform.
- Incomplete-splice-match (ISM) reads are excluded from the primary definition (fragments are
  compatible with more than one transcript); `ism_mode="assign"` adds them to their associated ENST.
- Denominator: all retained reads of the gene, including novel (NIC/NNC) isoforms, which are pooled
  into one `NOVEL` row per gene so the known-ENST proportions are not inflated.
- Genes are keyed by ENSG via the v31 GTF; fusion/antisense/intergenic/genic isoforms are dropped.
"""
import gzip
import re

import numpy as np
import pandas as pd

KNOWN = "full-splice_match"
FRAGMENT = "incomplete-splice_match"
NOVEL = ("novel_in_catalog", "novel_not_in_catalog")


def load_gtf_transcripts(path):
    """ENST (versioned) -> ENSG (unversioned), gene_name, transcript_type."""
    rows = []
    pat = re.compile(r'(gene_id|transcript_id|gene_name|transcript_type) "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.split("\t", 3)[2] != "transcript":
                continue
            a = dict(pat.findall(line))
            rows.append((a["transcript_id"], a["gene_id"].split(".")[0], a["gene_name"], a["transcript_type"]))
    return pd.DataFrame(rows, columns=["transcript_id", "gene_id", "gene_name", "transcript_type"])


def load_sqanti(class_path):
    c = pd.read_csv(class_path, sep="\t", low_memory=False)
    fl = [x for x in c.columns if x.startswith("FL.")]
    c = c.rename(columns={x: x[3:] for x in fl})
    return c, [x[3:] for x in fl]


def long_counts(c, samples, gtf, ism_mode="exclude"):
    """Collapse PB isoforms to (gene_id, feature) rows of per-sample FL counts."""
    tx2gene = gtf.set_index("transcript_id")["gene_id"]
    name2gene = gtf.drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
    # gene names shared by several ENSG cannot be resolved for novel isoforms
    ambiguous = gtf.groupby("gene_name")["gene_id"].nunique()
    name2gene = name2gene[ambiguous.reindex(name2gene.index) == 1]

    cats = [KNOWN, *NOVEL] + ([FRAGMENT] if ism_mode == "assign" else [])
    c = c[c["structural_category"].isin(cats)].copy()
    known = c["structural_category"].isin([KNOWN, FRAGMENT])
    c["gene_id"] = np.where(known, c["associated_transcript"].map(tx2gene), c["associated_gene"].map(name2gene))
    c["feature"] = np.where(known, c["associated_transcript"], "NOVEL")
    c = c.dropna(subset=["gene_id"])
    return c.groupby(["gene_id", "feature"])[samples].sum()


def proportions(counts):
    """Per-sample and pooled within-gene proportions; gene totals kept for depth filtering."""
    gene_tot = counts.groupby(level="gene_id").transform("sum")
    per_sample = counts / gene_tot.replace(0, np.nan)
    pooled_counts = counts.sum(axis=1)
    pooled_tot = pooled_counts.groupby(level="gene_id").transform("sum")
    out = per_sample.add_prefix("prop_")
    out["reads"] = pooled_counts
    out["gene_reads"] = pooled_tot
    out["prop_pooled"] = pooled_counts / pooled_tot
    out["n_samples_gene_detected"] = (gene_tot > 0).sum(axis=1)
    return out


def build(class_path, gtf, ism_mode="exclude"):
    c, samples = load_sqanti(class_path)
    return proportions(long_counts(c, samples, gtf, ism_mode)), samples

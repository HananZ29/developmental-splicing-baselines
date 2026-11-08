"""E15 pre-registration audit (no outcomes): annotation harmonization and isoform identifiability.

Inputs: GENCODE v10 (BrainSpan MRF), v19 (PsychENCODE DER-10), v31lift37 (our pipeline, src/splice_graph),
v33lift37 (Wen isoQTL / Salmon) GTFs, all hg19/GRCh37 coordinates.

Definitions:
- intron chain: the ordered tuple of introns (first and last intronic base, 1-based) of a multi-exon
  transcript, plus strand and chromosome.
- An "exact match" between annotations requires an identical intron chain. Single-exon transcripts are
  matched only if their exon coordinates are identical (reported separately).
- Matching is one-to-one only when the chain is unique in both annotations. Chains shared by several
  transcripts in one annotation (differing only in TSS/TES) are reported as "collapsible": they are
  indistinguishable by splicing and are merged into one splice isoform, never silently paired.
- Junction-level identifiability (within one annotation, per gene): two multi-exon transcripts are
  indistinguishable from junction reads if their intron chains are identical. A transcript is
  "junction-nested" if its intron set is a strict subset of another transcript's set (its reads are
  all compatible with the longer transcript too). EM proportions are weakly identified for those.

Output: publication_extension/results/E15_devstage/annotation_audit.json
"""
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / "data" / "external" / "gencode"
OUT = ROOT / "publication_extension" / "results" / "E15_devstage"
TAGS = {"v10": "gencode.v10.annotation.gtf.gz", "v19": "gencode.v19.annotation.gtf.gz",
        "v31lift37": "gencode.v31lift37.annotation.gtf.gz", "v33lift37": "gencode.v33lift37.annotation.gtf.gz"}
RX = {k: re.compile(k + r' "([^"]+)"') for k in ("transcript_id", "gene_id", "transcript_type", "gene_type")}


def load(path):
    ex = defaultdict(list)
    meta = {}
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("#"):
                continue
            p = l.rstrip("\n").split("\t")
            if p[2] != "exon":
                continue
            a = p[8]
            t = RX["transcript_id"].search(a).group(1).split("_")[0]
            ex[t].append((int(p[3]), int(p[4])))
            if t not in meta:
                tt = RX["transcript_type"].search(a)
                meta[t] = (p[0], p[6], RX["gene_id"].search(a).group(1).split(".")[0].split("_")[0],
                           tt.group(1) if tt else "NA")
    tx = {}
    for t, e in ex.items():
        e.sort()
        chrom, strand, gene, ttype = meta[t]
        introns = tuple((e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1))
        tx[t] = dict(chrom=chrom, strand=strand, gene=gene, type=ttype, introns=introns,
                     single=(e[0] if len(e) == 1 else None))
    return tx


def chain_key(r):
    return (r["chrom"], r["strand"], r["introns"]) if r["introns"] else None


def identifiability(tx):
    by_gene = defaultdict(list)
    for t, r in tx.items():
        if r["introns"]:
            by_gene[r["gene"]].append(r)
    n_multi = dup = nested = 0
    genes_multi_iso = genes_any_issue = 0
    for g, rs in by_gene.items():
        if len(rs) < 2:
            n_multi += len(rs)
            continue
        genes_multi_iso += 1
        chains = [frozenset(r["introns"]) for r in rs]
        seen = defaultdict(int)
        for c in chains:
            seen[c] += 1
        issue = False
        for c in chains:
            n_multi += 1
            if seen[c] > 1:
                dup += 1
                issue = True
            elif any(c < d for d in chains):
                nested += 1
                issue = True
        genes_any_issue += issue
    return {"multi_exon_transcripts": n_multi, "identical_chain_within_gene": dup,
            "junction_nested_strict_subset": nested,
            "genes_with_>=2_multiexon_tx": genes_multi_iso, "genes_with_any_unidentifiable_pair": genes_any_issue}


def compare(src, dst):
    dchains = defaultdict(list)
    for t, r in dst.items():
        k = chain_key(r)
        if k:
            dchains[k].append(t)
    schains = defaultdict(list)
    for t, r in src.items():
        k = chain_key(r)
        if k:
            schains[k].append(t)
    exact_unique = exact_collapsible = none = same_id = id_diff_chain = 0
    for t, r in src.items():
        k = chain_key(r)
        if not k:
            continue
        hits = dchains.get(k, [])
        if not hits:
            none += 1
        elif len(hits) == 1 and len(schains[k]) == 1:
            exact_unique += 1
        else:
            exact_collapsible += 1
        if hits and t.split(".")[0] in {h.split(".")[0] for h in hits}:
            same_id += 1
    # stable-ID lookups whose chain differs (the dangerous case for ID-only matching)
    dst_by_stable = {t.split(".")[0]: r for t, r in dst.items()}
    for t, r in src.items():
        d = dst_by_stable.get(t.split(".")[0])
        if d is not None and r["introns"] and chain_key(d) != chain_key(r):
            id_diff_chain += 1
    n = sum(1 for r in src.values() if r["introns"])
    return {"multi_exon_src": n, "exact_chain_unique_1to1": exact_unique,
            "exact_chain_but_collapsible_TSS_TES_variants": exact_collapsible,
            "no_exact_chain_in_dst": none, "same_stable_id_and_same_chain": same_id,
            "same_stable_id_but_DIFFERENT_chain": id_diff_chain}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    A = {k: load(G / v) for k, v in TAGS.items()}
    rep = {"identifiability": {k: identifiability(v) for k, v in A.items()}}
    rep["to_v31lift37"] = {k: compare(A[k], A["v31lift37"]) for k in ("v10", "v19", "v33lift37")}
    (OUT / "annotation_audit.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()

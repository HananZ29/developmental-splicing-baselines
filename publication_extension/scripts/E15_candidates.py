"""E15 candidate SNVs for SpliceAI scoring (outcome-blind superset of the eligible population).

Uses only: Wen intron IDs, Wen fine-mapping columns 0-7 (phenotype_id, variant_id, chr, pos, ref, alt,
cs_id, cs_index) - never pip, z, posterior_* or cs_log10bf - GENCODE v31lift37 and the GTEx junction
presence list. Coverage rules (E2/E3) are applied later, after the BrainSpan extraction completes, so
this set is a superset of what will be analysed.

A candidate is a biallelic SNV in a credible set whose phenotype intron is an annotated intron of a
cluster passing E1 (>= 2 annotated introns, full-rank pattern matrix) and E4 (all annotated introns in
GTEx), lying within +-50 nt of an annotated splice site of that cluster, with REF matching hg19.
Strand = strand of the annotated transcripts using the cluster's introns (clusters with conflicting
strands are excluded and counted).

Output: publication_extension/results/E15_devstage/candidates.tsv (one row per CS member candidate)
        publication_extension/results/E15_devstage/to_score.tsv (unique chrom,pos,ref,alt,strand)
"""
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_baseline_diagnostics as D  # noqa: E402

SUSIE = D.ROOT / "data/external/wen_synapse/mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz"
FA = D.ROOT / "data/external/genome/hg19.fa"
WIN = 50


def strands_by_intron():
    s = defaultdict(set)
    with gzip.open(D.GTF, "rt") as f:
        ex = defaultdict(list)
        meta = {}
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "exon":
                continue
            t = p[8].split('transcript_id "')[1].split('"')[0]
            ex[t].append((int(p[3]), int(p[4])))
            meta[t] = (p[0], p[6])
    for t, e in ex.items():
        e.sort()
        c, st = meta[t]
        for i in range(len(e) - 1):
            s[(c, e[i][1] + 1, e[i + 1][0] - 1)].add(st)
    return s


def main():
    W = D.wen_universe()
    tx = D.transcript_introns()
    ann_set = set().union(*tx.values())
    keys = list(zip(W.chrom.tolist(), W.start.tolist(), W.end.tolist()))
    ann = np.array([k in ann_set for k in keys])
    K = D.Clusters(W, ann)
    pats = D.patterns_for_clusters(W, K, tx)
    e1 = (K.tot(K.ann) >= 2) & np.array([D.full_rank(p, n) for p, n in zip(pats, K.size)])
    g = pd.read_csv(D.BDIR / "gtex_tissue_counts.tsv.gz", sep="\t", usecols=["chrom", "start", "end"])
    present = set(zip(g.chrom, g.start, g.end))
    gp = np.array([k in present for k in keys])
    e4 = K.tot((K.ann * ~gp).astype(np.int64)) == 0
    strand_of = strands_by_intron()
    st = defaultdict(int)
    ok_cluster, cl_strand, cl_sites, intron_cluster = set(), {}, {}, {}
    for ci in np.flatnonzero(e1 & e4):
        a, n = K.starts[ci], K.size[ci]
        idx = [i for i in range(a, a + n) if ann[i]]
        ss = set().union(*(strand_of[keys[i]] for i in idx))
        if len(ss) != 1:
            st["clusters_conflicting_strand"] += 1
            continue
        s = ss.pop()
        cid = K.ids[ci]
        ok_cluster.add(cid)
        cl_strand[cid] = s
        cl_sites[cid] = sorted({keys[i][1] - 1 for i in idx} | {keys[i][2] + 1 for i in idx})
        for i in idx:
            intron_cluster[keys[i]] = cid
    st["clusters_E1"] = int(e1.sum())
    st["clusters_E1_E4"] = int((e1 & e4).sum())
    st["clusters_E1_E4_single_strand"] = len(ok_cluster)
    fa = Fasta(str(FA))
    rows = []
    with gzip.open(SUSIE, "rt") as f:
        f.readline()
        for line in f:
            p = line.split("\t", 8)
            ph, vid, c, pos, ref, alt, cs, csi = p[0], p[1], p[2], int(p[3]), p[4], p[5], p[6], p[7]
            st["susie_rows"] += 1
            pc, pa, pb = ph.split(":")[:3]
            k = ("chr" + pc.replace("chr", ""), int(pa) + 1, int(pb) - 1)
            cid = intron_cluster.get(k)
            if cid is None:
                continue
            st["rows_phenotype_in_E1_E4_cluster_annotated_intron"] += 1
            if len(ref) != 1 or len(alt) != 1:
                st["rows_not_SNV"] += 1
                continue
            sites = cl_sites[cid]
            dist = min(abs(pos - x) for x in sites)
            if dist > WIN:
                continue
            chrom = "chr" + c.replace("chr", "")
            if str(fa[chrom][pos - 1]).upper() != ref.upper():
                st["rows_REF_mismatch_hg19"] += 1
                continue
            rows.append((ph, cs, csi, vid, chrom, pos, ref, alt, cid, cl_strand[cid], dist, k[0], k[1], k[2]))
    df = pd.DataFrame(rows, columns=["phenotype_id", "cs_id", "cs_index", "variant_id", "chrom", "pos", "ref", "alt",
                                     "cluster", "strand", "dist_to_site", "ph_chrom", "ph_start", "ph_end"])
    out = D.RES
    df.to_csv(out / "candidates.tsv", sep="\t", index=False)
    u = df[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().sort_values(["chrom", "pos"])
    u.to_csv(out / "to_score.tsv", sep="\t", index=False)
    st["candidate_rows"] = len(df)
    st["candidate_credible_sets"] = int(df[["phenotype_id", "cs_id"]].drop_duplicates().shape[0])
    st["unique_snvs_to_score"] = len(u)
    (out / "candidates_stats.json").write_text(json.dumps(dict(st), indent=1))
    print(json.dumps(dict(st), indent=1))


if __name__ == "__main__":
    main()

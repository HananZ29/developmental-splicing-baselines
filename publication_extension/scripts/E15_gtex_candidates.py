"""E15 secondary crossover arm: candidate SNVs from GTEx v8 Brain_Cortex significant sQTL pairs.

Reads ONLY columns 1-2 (variant_id, phenotype_id). No slope, p-value or other statistic is read.
- GTEx phenotype intron chr:a:b (hg38) = intron first base a+1, last base b-1 (verified: 81.3% of the
  1,937 unique introns match GENCODE v26 at this offset; no other offset matches > 0.1%).
- Intron and variant are lifted to hg19 with the same rule as the G baseline (unique, same
  chromosome/strand, length unchanged). The intron must equal an annotated intron of an E15-eligible
  cluster (cluster_eligibility.tsv.gz).
- Variant: biallelic SNV within +-50 nt of an annotated splice site of that cluster. Alleles are
  oriented to hg19: if GTEx REF equals hg19 base -> as is; if GTEx ALT equals hg19 base -> swapped
  (truth slope sign is flipped in the analysis); otherwise dropped.
Unit (crossover) = one GTEx phenotype intron; representative variant chosen by rule P (max SpliceAI
DS_max) or D (closest to a splice site) among that intron's candidates.

Outputs (publication_extension/results/E15_devstage/): gtex_candidates.tsv, gtex_to_score.tsv
(SNVs not already scored), gtex_candidates_stats.json
"""
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from pyliftover import LiftOver

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_baseline_diagnostics as D  # noqa: E402
import E15_candidates as C  # noqa: E402

SQTL = D.ROOT / "data/external/drive_new/Brain_Cortex.v8.sqtl_signifpairs.txt.gz"
CHAIN = D.ROOT / "data/external/genome/hg38ToHg19.over.chain.gz"


def lift1(lo, c, pos):
    r = lo.convert_coordinate(c, pos - 1)
    if not r or len(r) != 1 or r[0][2] != "+":
        return None
    return r[0][0], r[0][1] + 1


def main():
    st = Counter()
    W = D.wen_universe()
    tx = D.transcript_introns()
    ann_set = set().union(*tx.values())
    keys = list(zip(W.chrom.tolist(), W.start.tolist(), W.end.tolist()))
    ann = np.array([k in ann_set for k in keys])
    K = D.Clusters(W, ann)
    el = pd.read_csv(D.BDIR / "cluster_eligibility.tsv.gz", sep="\t")
    elig = set(el.cluster[el.eligible])
    strand_of = C.strands_by_intron()
    intron_cluster, cl_sites, cl_strand = {}, {}, {}
    for ci, cid in enumerate(K.ids):
        if cid not in elig:
            continue
        a, n = K.starts[ci], K.size[ci]
        idx = [i for i in range(a, a + n) if ann[i]]
        ss = set().union(*(strand_of[keys[i]] for i in idx))
        if len(ss) != 1:
            continue
        cl_strand[cid] = ss.pop()
        cl_sites[cid] = sorted({keys[i][1] - 1 for i in idx} | {keys[i][2] + 1 for i in idx})
        for i in idx:
            intron_cluster[keys[i]] = cid
    lo = LiftOver(str(CHAIN))
    fa = Fasta(str(C.FA))
    rows = []
    with gzip.open(SQTL, "rt") as f:
        f.readline()
        for line in f:
            vid, ph = line.split("\t", 2)[:2]
            st["pairs"] += 1
            c, a, b = ph.split(":")[:3]
            s38, e38 = int(a) + 1, int(b) - 1
            ls, le = lift1(lo, c, s38), lift1(lo, c, e38)
            if not ls or not le or ls[0] != le[0] or le[1] - ls[1] != e38 - s38:
                st["intron_liftover_fail"] += 1
                continue
            k = (ls[0], ls[1], le[1])
            cid = intron_cluster.get(k)
            if cid is None:
                st["intron_not_in_eligible_cluster"] += 1
                continue
            vc, vp, ref, alt = vid.split("_")[:4]
            if len(ref) != 1 or len(alt) != 1:
                st["not_SNV"] += 1
                continue
            lv = lift1(lo, vc, int(vp))
            if not lv or lv[0] != k[0]:
                st["variant_liftover_fail"] += 1
                continue
            pos = lv[1]
            dist = min(abs(pos - x) for x in cl_sites[cid])
            if dist > C.WIN:
                continue
            base = str(fa[lv[0]][pos - 1]).upper()
            if base == ref:
                r19, a19, flipped = ref, alt, 0
            elif base == alt:
                r19, a19, flipped = alt, ref, 1
            else:
                st["allele_mismatch_hg19"] += 1
                continue
            rows.append((ph, vid, lv[0], pos, r19, a19, flipped, cid, cl_strand[cid], dist, k[0], k[1], k[2]))
    df = pd.DataFrame(rows, columns=["phenotype_id", "variant_id", "chrom", "pos", "ref", "alt", "allele_flipped",
                                     "cluster", "strand", "dist_to_site", "ph_chrom", "ph_start", "ph_end"])
    df.to_csv(D.RES / "gtex_candidates.tsv", sep="\t", index=False)
    done = pd.read_csv(D.RES / "spliceai/scores.tsv", sep="\t")[["chrom", "pos", "ref", "alt", "strand"]]
    u = df[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates()
    u = u.merge(done, how="left", indicator=True)
    u = u[u._merge == "left_only"].drop(columns="_merge").sort_values(["chrom", "pos"])
    u.to_csv(D.RES / "gtex_to_score.tsv", sep="\t", index=False)
    st["candidate_rows"] = len(df)
    st["candidate_introns(units)"] = int(df.phenotype_id.nunique())
    st["snvs_needing_spliceai"] = len(u)
    (D.RES / "gtex_candidates_stats.json").write_text(json.dumps(dict(st), indent=1))
    print(json.dumps(dict(st), indent=1))


if __name__ == "__main__":
    main()

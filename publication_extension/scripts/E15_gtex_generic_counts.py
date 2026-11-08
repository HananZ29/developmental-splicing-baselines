"""E15 generic baseline (G): GTEx v8 junction read counts on the Wen intron universe, summed per tissue.

No QTL outcome is read. Inputs:
- GTEx v8 STAR junction GCT (hg38; junction ID chr_<first intron base>_<last intron base>, 1-based)
- GTEx v8 sample attributes (SMTSD = tissue)
- Wen all.introns.tested.tsv (hg19; Wen chr:a:b = intron first base a+1, last base b-1; verified)
- UCSC hg38ToHg19 chain

Liftover rule (R5): both intron ends lift to a single hg19 position on the same chromosome and strand,
and intron length is unchanged; otherwise the junction is dropped from G (counted).

Also measured: the fraction of GTEx junctions that are annotated in GENCODE v26 (GTEx v8's annotation),
so we know whether unannotated Wen introns can ever receive reads in G.

Output (data/processed/E15_baselines/):
- gtex_tissue_counts.tsv.gz: chrom, start, end (hg19, 1-based intron), one column per tissue (summed reads)
- gtex_tissue_samples.tsv: tissue, n_samples, is_brain
- gtex_liftover_stats.json
"""
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
GCT = ROOT / "data/external/gtex_junctions/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
ATTR = ROOT / "data/external/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
WEN = ROOT / "data/external/wen_synapse/all.introns.tested.tsv"
CHAIN = ROOT / "data/external/genome/hg38ToHg19.over.chain.gz"
V26 = ROOT / "data/external/gencode/gencode.v26.annotation.gtf.gz"
OUT = ROOT / "data/processed/E15_baselines"


def wen_universe():
    s = set()
    with open(WEN) as f:
        for line in f:
            c, a, b = line.split("\t")[0].split(":")[:3]
            s.add(("chr" + c.replace("chr", ""), int(a) + 1, int(b) - 1))
    return s


def annotated_introns(path):
    ex = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "exon":
                continue
            t = p[8].split('transcript_id "')[1].split('"')[0]
            ex.setdefault(t, [p[0], []])[1].append((int(p[3]), int(p[4])))
    s = set()
    for c, e in ex.values():
        e.sort()
        for i in range(len(e) - 1):
            s.add((c, e[i][1] + 1, e[i + 1][0] - 1))
    return s


def lift(lo, c, s, e):
    a = lo.convert_coordinate(c, s - 1)
    b = lo.convert_coordinate(c, e - 1)
    if not a or not b or len(a) != 1 or len(b) != 1:
        return None, "unmapped_or_multi"
    (ca, pa, sa, _), (cb, pb, sb, _) = a[0], b[0]
    if ca != cb or sa != "+" or sb != "+":
        return None, "chrom_or_strand_change"
    s19, e19 = pa + 1, pb + 1
    if e19 - s19 != e - s:
        return None, "length_change"
    return (ca, s19, e19), "ok"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    W = wen_universe()
    ann26 = annotated_introns(V26)
    lo = LiftOver(str(CHAIN))
    tissue_of = {}
    with open(ATTR) as f:
        h = f.readline().rstrip("\n").split("\t")
        i_id, i_t = h.index("SAMPID"), h.index("SMTSD")
        for line in f:
            p = line.rstrip("\n").split("\t")
            tissue_of[p[i_id]] = p[i_t]
    st = Counter()
    rows = {}
    with gzip.open(GCT, "rt") as f:
        f.readline(); f.readline()
        samples = f.readline().rstrip("\n").split("\t")[2:]
        tissues = sorted({tissue_of[s] for s in samples})
        tidx = np.array([tissues.index(tissue_of[s]) for s in samples])
        seen_ids = set()
        collided = set()
        for line in f:
            jid, rest = line.split("\t", 1)
            if jid in seen_ids:                 # same junction listed again for another gene
                st["duplicate_junction_rows_skipped"] += 1
                continue
            seen_ids.add(jid)
            c, s, e = jid.rsplit("_", 2)
            s, e = int(s), int(e)
            st["gtex_junctions"] += 1
            st["gtex_annotated_v26"] += (c, s, e) in ann26
            key, why = lift(lo, c, s, e)
            st["liftover_" + why] += 1
            if key is None or key not in W:
                continue
            counts = np.array(rest.split("\t", 1)[1].split("\t"), dtype=np.int64)
            agg = np.bincount(tidx, weights=counts, minlength=len(tissues)).astype(np.int64)
            if key in rows:
                st["collisions_after_liftover"] += 1
                collided.add(key)
            else:
                rows[key] = agg
    # distinct hg38 junctions lifting onto one hg19 intron: ambiguous, dropped from G (E4 then
    # excludes the cluster)
    for k in collided:
        rows.pop(k, None)
    st["hg19_introns_dropped_ambiguous_liftover"] = len(collided)
    st["wen_introns"] = len(W)
    st["wen_introns_present_in_gtex"] = len(rows)
    st["wen_introns_with_reads_in_gtex"] = sum(1 for v in rows.values() if v.sum() > 0)
    with gzip.open(OUT / "gtex_tissue_counts.tsv.gz", "wt") as fo:
        fo.write("chrom\tstart\tend\t" + "\t".join(tissues) + "\n")
        for k in sorted(rows):
            fo.write(f"{k[0]}\t{k[1]}\t{k[2]}\t" + "\t".join(map(str, rows[k])) + "\n")
    n_per = Counter(tissue_of[s] for s in samples)
    with open(OUT / "gtex_tissue_samples.tsv", "w") as fo:
        fo.write("tissue\tn_samples\tis_brain\n")
        for t in tissues:
            fo.write(f"{t}\t{n_per[t]}\t{int(t.startswith('Brain'))}\n")
    st["n_samples"] = len(samples)
    st["n_tissues"] = len(tissues)
    (OUT / "gtex_liftover_stats.json").write_text(json.dumps(dict(st), indent=1))
    print(json.dumps(dict(st), indent=1))


if __name__ == "__main__":
    main()

"""Validate the SpliceAI implementation and the splice-site coordinate convention on hg19.

Positive controls: for random annotated introns (GENCODE v31lift37, protein-coding, autosomes), mutate
the invariant intronic G of the donor GT (+1) or of the acceptor AG (-1) to A on the transcript strand.
Expectations if coordinates, strand handling and model loading are right:
- the REF donor/acceptor probability peaks exactly at our annotated donor/acceptor position (offset 0)
- the canonical-dinucleotide mutation gives large DS_DL / DS_AL (published median >0.9 for such sites)

Usage: python experiments/05_validate_spliceai.py
Output: results/spliceai_validation.json
"""
import tensorflow  # noqa: F401  (load the TF runtime before other native libraries)
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from splice_graph import load_lift37  # noqa: E402
from spliceai_runner import SpliceAIScorer  # noqa: E402

EXT = ROOT / "data" / "external"
COMP = str.maketrans("ACGT", "TGCA")


def protein_coding_genes(path, n, rng):
    genes = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if "\tgene\t" in line and 'gene_type "protein_coding"' in line and re.match(r"chr\d+\t", line):
                genes.add(re.search(r'gene_id "([^".]+)', line).group(1))
    return set(rng.choice(sorted(genes), n, replace=False))


def main():
    rng = np.random.default_rng(1)
    gtf = EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz"
    by_gene = load_lift37(gtf, protein_coding_genes(gtf, 60, rng))
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    W = 500
    scorer = SpliceAIScorer(fa, str(EXT / "spliceai"), window=W)
    variants = []
    for g, ts in by_gene.items():
        t = max(ts, key=lambda t: len(t.exons))
        if len(t.introns) < 3:
            continue
        s, e = t.introns[len(t.introns) // 2]
        if t.strand == "+":
            donor_g, acc_g, donor_site, acc_site = s, e, s - 1, e + 1
        else:
            donor_g, acc_g, donor_site, acc_site = e, s, e + 1, s - 1
        for kind, gpos, site in [("donor", donor_g, donor_site), ("acceptor", acc_g, acc_site)]:
            ref = fa[t.chrom][gpos - 1:gpos].seq.upper()
            want = "G" if t.strand == "+" else "C"
            if ref != want:
                continue  # non-canonical site
            alt = "A" if t.strand == "+" else "T"
            variants.append(dict(chrom=t.chrom, pos=gpos, ref=ref, alt=alt, strand=t.strand, kind=kind, site=site))
    res = []
    for i in range(0, len(variants), 8):
        res += scorer.score(variants[i:i + 8])
    offsets = {"donor": [], "acceptor": []}
    ds = {"donor": [], "acceptor": []}
    for v, r in zip(variants, res):
        track = r["don_ref"] if v["kind"] == "donor" else r["acc_ref"]
        centre = W + (v["site"] - v["pos"])
        win = track[centre - 5:centre + 6]
        offsets[v["kind"]].append(int(np.argmax(win)) - 5)
        ds[v["kind"]].append(r["DS_DL"] if v["kind"] == "donor" else r["DS_AL"])
        v["ref_prob_at_site"] = float(track[centre])
    out = {}
    for k in ["donor", "acceptor"]:
        o, d = np.array(offsets[k]), np.array(ds[k])
        p = np.array([v["ref_prob_at_site"] for v in variants if v["kind"] == k])
        out[k] = {"n": int(len(o)), "frac_peak_at_offset0": round(float((o == 0).mean()), 3),
                  "offset_counts": {int(a): int(b) for a, b in zip(*np.unique(o, return_counts=True))},
                  "median_ref_prob_at_site": round(float(np.median(p)), 3),
                  "median_loss_score": round(float(np.median(d)), 3), "frac_loss_gt_0.5": round(float((d > 0.5).mean()), 3)}
    out["strands"] = {s: sum(v["strand"] == s for v in variants) for s in "+-"}
    (ROOT / "results" / "spliceai_validation.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

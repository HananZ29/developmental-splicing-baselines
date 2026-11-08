"""hg19 transcript structures (GENCODE v31lift37) as per-gene splice graphs.

A transcript is its ordered list of introns; each intron is (donor, acceptor) in hg19 1-based
coordinates of the *exonic* bases flanking the intron:
  + strand: donor = last base of upstream exon (intron_start - 1), acceptor = first base of
            downstream exon (intron_end + 1)
  - strand: donor = intron_end + 1, acceptor = intron_start - 1 (transcript orientation)
`intron` keys are (intron_start, intron_end), the first and last intronic bases, which equals the
leafcutter convention start+1..end-1 verified by experiments/04_genome_build_checks.py.
"""
import gzip
import re
from dataclasses import dataclass, field

import pandas as pd

ATTR = re.compile(r'(gene_id|transcript_id|gene_type|transcript_type|gene_name) "([^"]+)"')


@dataclass
class Transcript:
    enst: str
    ensg: str
    chrom: str
    strand: str
    exons: list = field(default_factory=list)

    @property
    def introns(self):
        ex = sorted(self.exons)
        return [(a[1] + 1, b[0] - 1) for a, b in zip(ex, ex[1:])]

    def sites(self):
        """{(position, 'donor'|'acceptor')} for all splice sites of the transcript."""
        out = set()
        for s, e in self.introns:
            if self.strand == "+":
                out |= {(s - 1, "donor"), (e + 1, "acceptor")}
            else:
                out |= {(e + 1, "donor"), (s - 1, "acceptor")}
        return out


def load_lift37(path, genes=None):
    """Return {ensg: [Transcript, ...]} with unversioned ids (lift37 '_k' suffixes stripped)."""
    tx = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            ensg = a["gene_id"].split(".")[0]
            if genes is not None and ensg not in genes:
                continue
            enst = a["transcript_id"].split(".")[0]
            t = tx.get(enst)
            if t is None:
                t = tx[enst] = Transcript(enst, ensg, f[0], f[6])
            t.exons.append((int(f[3]), int(f[4])))
    by_gene = {}
    for t in tx.values():
        by_gene.setdefault(t.ensg, []).append(t)
    return by_gene


def site_table(by_gene):
    rows = [(g, t.chrom, t.strand, pos, kind) for g, ts in by_gene.items() for t in ts for pos, kind in t.sites()]
    return pd.DataFrame(rows, columns=["ensg", "chrom", "strand", "pos", "kind"]).drop_duplicates()

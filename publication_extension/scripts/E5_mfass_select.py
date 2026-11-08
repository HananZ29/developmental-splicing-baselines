"""E5 (pre-registered): choose MFASS variants for the measured-effect test of the propagation equation.

Rules (PREREGISTRATION.md, E5):
- category == "mutant", non-missing v2_dpsi and nat_v2_index
- hg19 base at snp_position equals ref_allele (+ strand) or its complement (decided globally by the
  higher match rate); variants that match neither are dropped
- random sample of 3,000 (numpy seed 20261104), drawn before any scoring

Exon sites (1-based, hg19; tile = intron1 + exon + intron2, start..end inclusive):
  + strand: acceptor = start + intron1_len, donor = acceptor + exon_len - 1
  - strand: acceptor = end - intron1_len,   donor = acceptor - exon_len + 1
(this convention is checked against REF SpliceAI peaks after scoring, using REF tracks only)

Output: publication_extension/results/E5_mfass/mfass_variants.tsv, selection.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension" / "results" / "E5_mfass"
N, SEED = 3000, 20261104
COMP = str.maketrans("ACGT", "TGCA")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(ROOT / "data" / "external" / "mfass" / "snv_data_clean.txt", sep="\t", low_memory=False)
    info = {"rows": len(m)}
    m = m[(m.category == "mutant") & m.v2_dpsi.notna() & m.nat_v2_index.notna() & m.snp_position.notna()].copy()
    m = m[m.ref_allele.str.len().eq(1) & m.alt_allele.str.len().eq(1)]
    info["mutants_with_measurements"] = len(m)
    fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))
    m["pos"] = m.snp_position.astype(int)
    m["hg19_base"] = [fa[c][p - 1:p].seq.upper() for c, p in zip(m.chr, m.pos)]
    plus = (m.hg19_base == m.ref_allele).mean()
    minus = (m.hg19_base == m.ref_allele.str.translate(COMP)).mean()
    info["ref_match_plus_strand"], info["ref_match_complement"] = round(float(plus), 4), round(float(minus), 4)
    if plus >= minus:
        m["ref"], m["alt"] = m.ref_allele, m.alt_allele
    else:
        m["ref"], m["alt"] = m.ref_allele.str.translate(COMP), m.alt_allele.str.translate(COMP)
    m = m[m.hg19_base == m.ref]
    info["ref_verified"] = len(m)
    m["acceptor"] = np.where(m.strand == "+", m.start + m.intron1_len, m.end - m.intron1_len).astype(int)
    m["donor"] = np.where(m.strand == "+", m.acceptor + m.exon_len - 1, m.acceptor - m.exon_len + 1).astype(int)
    m["exon_id"] = m.ensembl_id
    rng = np.random.default_rng(SEED)
    s = m.iloc[np.sort(rng.choice(len(m), min(N, len(m)), replace=False))].copy()
    s["chrom"] = s.chr
    cols = ["id", "exon_id", "chrom", "pos", "ref", "alt", "strand", "acceptor", "donor", "exon_len", "label",
            "rel_position", "nat_v2_index", "v2_index", "v2_dpsi", "strong_lof"]
    s[cols].to_csv(OUT / "mfass_variants.tsv", sep="\t", index=False)
    info |= {"sampled": len(s), "exons_sampled": int(s.exon_id.nunique()),
             "nat_v2_index_quantiles": s.nat_v2_index.quantile([0, .1, .5, .9, 1]).round(3).to_dict(),
             "v2_dpsi_quantiles": s.v2_dpsi.quantile([0, .05, .5, .95, 1]).round(3).to_dict(),
             "abs_dpsi_ge_0.1": int((s.v2_dpsi.abs() >= 0.1).sum()), "strong_lof": int(s.strong_lof.fillna(False).sum()),
             "label_counts": s.label.value_counts().to_dict()}
    (OUT / "selection.json").write_text(json.dumps(info, indent=2, default=str))
    print(json.dumps(info, indent=2, default=str))


if __name__ == "__main__":
    main()

"""MMSplice (Cheng et al. 2019, Genome Biol; mmsplice 2.4.0 from PyPI, MIT) exon-level scoring.

For each (variant, exon) pair:
- Build the overhanged exon sequence: 100 nt intron + exon + 100 nt intron, MMSplice's default
  overhang, in transcript orientation from hg19. Do this for REF and ALT.
- Obtain the five modular scores with MMSplice().predict_on_seq.
- Apply MMSplice's own linear model: dlogit_psi = predict_deltaLogitPsi(X_ref, X_alt).
- Downstream, the predicted delta-PSI for a baseline psi is
  mmsplice.utils.delta_logit_PSI_to_delta_PSI(dlogit_psi, psi).

Variants farther than 100 nt from the exon are outside MMSplice's input. They get dlogit_psi = 0 and
in_range = False, as the VCF dataloader would not score them.

Input table (TSV): chrom, pos, ref, alt, strand, exon_start, exon_end.
- hg19, 1-based inclusive exon coordinates.
- REF/ALT on the + strand.
- SNVs only.

Usage: <mmsplice venv python> -E score_mmsplice.py <in.tsv> <out.tsv>

Note: mmsplice is installed without cyvcf2, which does not build on Windows. Only the
sequence-level API is used, so the VCF loader is not needed.
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from pyfaidx import Fasta  # noqa: E402
from mmsplice import MMSplice  # noqa: E402
from mmsplice.utils import predict_deltaLogitPsi  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OH = 100
COMP = str.maketrans("ACGTN", "TGCAN")
MOD = ["acceptor_intron", "acceptor", "exon", "donor", "donor_intron"]


def seqs(fa, r):
    s0 = r.exon_start - 1 - OH                  # 0-based start of the overhanged window
    e0 = r.exon_end + OH                        # 0-based exclusive end
    ref = fa[r.chrom][s0:e0].seq.upper()
    i = r.pos - 1 - s0
    if not (0 <= i < len(ref)):
        return None
    if ref[i] != r.ref:
        raise ValueError(f"REF mismatch {r.chrom}:{r.pos}")
    alt = ref[:i] + r.alt + ref[i + 1:]
    if r.strand == "-":
        ref, alt = ref.translate(COMP)[::-1], alt.translate(COMP)[::-1]
    return ref, alt


def main(inp, out):
    v = pd.read_csv(inp, sep="\t")
    fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))
    m = MMSplice()
    rows = []
    for r in v.itertuples(index=False):
        s = seqs(fa, r)
        rec = dict(chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt, strand=r.strand,
                   exon_start=r.exon_start, exon_end=r.exon_end, in_range=s is not None)
        if s is None:
            rec["mmsplice_dlogit_psi"] = 0.0
        else:
            xr = np.asarray(m.predict_on_seq(s[0], overhang=(OH, OH))).reshape(1, -1)
            xa = np.asarray(m.predict_on_seq(s[1], overhang=(OH, OH))).reshape(1, -1)
            rec["mmsplice_dlogit_psi"] = float(predict_deltaLogitPsi(xr, xa)[0])
            rec.update({f"ref_{k}": float(x) for k, x in zip(MOD, xr[0])})
            rec.update({f"alt_{k}": float(x) for k, x in zip(MOD, xa[0])})
        rows.append(rec)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"scored {len(rows)}; in range {sum(x['in_range'] for x in rows)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])

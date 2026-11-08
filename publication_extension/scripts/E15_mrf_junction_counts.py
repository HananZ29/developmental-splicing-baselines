"""E15 baseline construction step 1: per-sample splice-junction read counts from BrainSpan MRF entries.

Input: data/external/brainspan_mrf/neocortex/<entry>.zipentry plus manifest.tsv (from
E15_fetch_brainspan_mrf.py). Each line of an MRF file is one read or read pair:
- mates are separated by '|';
- each mate is a list of alignment blocks separated by ',';
- each block is chr:strand:start:end:qstart:qend (1-based, GRCh37).

Rules (fixed before any use):
- A junction is the gap between consecutive blocks of one mate, on the same chromosome, with genomic
  gap >= MIN_INTRON nt. Shorter gaps are treated as deletions and ignored.
- Junction coordinates are the intron's first and last bases, 1-based inclusive (STAR SJ.out
  convention): intron_start = block_end + 1, intron_end = next_block_start - 1.
- Reads are counted once per junction per mate. Mates of one pair crossing the same junction count
  once (a set per read line).
- Libraries are treated as unstranded. Strand is assigned later from the GT-AG/CT-AC motif and the
  annotation, not from the MRF read strand.
- The CRC-32 of the decompressed MRF must equal the zip central-directory CRC, or the sample is
  rejected.

Output: data/processed/E15_junctions/<entry>.junctions.tsv.gz (chrom, intron_start, intron_end, reads)
and <entry>.stats.json. Resumable.

Usage: python -I E15_mrf_junction_counts.py [entry_name ...]   (default: all entries in the manifest)
"""
import gzip
import json
import sys
import zlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "external" / "brainspan_mrf" / "neocortex"
OUT = ROOT / "data" / "processed" / "E15_junctions"
MIN_INTRON = 20
CH = 16 << 20


def lines_from_entry(path, crc_expected):
    with path.open("rb") as fh:
        lh = fh.read(30)
        n, e = int.from_bytes(lh[26:28], "little"), int.from_bytes(lh[28:30], "little")
        fh.read(n + e)
        dec, crc, buf = zlib.decompressobj(-15), 0, b""
        while True:
            chunk = fh.read(CH)
            x = dec.decompress(chunk) if chunk else dec.flush()
            crc = zlib.crc32(x, crc)
            buf += x
            parts = buf.split(b"\n")
            buf = parts.pop()
            yield from parts
            if not chunk:
                break
        if buf:
            yield buf
    if f"{crc:08x}" != crc_expected:
        raise ValueError(f"CRC mismatch for {path.name}")


def count(entry, crc):
    J, n_reads, n_spliced = Counter(), 0, 0
    it = lines_from_entry(SRC / (entry + ".zipentry"), crc)
    header = next(it)
    if not header.startswith(b"AlignmentBlocks"):
        raise ValueError(f"unexpected MRF header {header[:60]!r}")
    for line in it:
        if not line:
            continue
        n_reads += 1
        if b"," not in line:
            continue
        seen = set()
        for mate in line.split(b"|"):
            bl = mate.split(b",")
            if len(bl) < 2:
                continue
            prev = bl[0].split(b":")
            for b in bl[1:]:
                cur = b.split(b":")
                if cur[0] == prev[0]:
                    s, e = int(prev[3]) + 1, int(cur[2]) - 1
                    if e - s + 1 >= MIN_INTRON:
                        seen.add((prev[0], s, e))
                prev = cur
        if seen:
            n_spliced += 1
            for k in seen:
                J[k] += 1
    return J, n_reads, n_spliced


def main(names):
    OUT.mkdir(parents=True, exist_ok=True)
    man = [l.split("\t") for l in (SRC / "manifest.tsv").read_text().splitlines()[1:]]
    crcs = {m[1]: m[3] for m in man}
    for entry in names or list(crcs):
        out = OUT / (entry + ".junctions.tsv.gz")
        if out.exists():
            continue
        J, n, ns = count(entry, crcs[entry])
        tmp = out.with_suffix(".part")
        with gzip.open(tmp, "wt") as fh:
            fh.write("chrom\tintron_start\tintron_end\treads\n")
            for (c, s, e), r in sorted(J.items()):
                fh.write(f"{c.decode()}\t{s}\t{e}\t{r}\n")
        tmp.replace(out)
        stats = {"entry": entry, "reads_or_pairs": n, "spliced_reads": ns, "junctions": len(J),
                 "junctions_ge5_reads": sum(1 for v in J.values() if v >= 5), "crc_ok": True,
                 "min_intron": MIN_INTRON}
        (OUT / (entry + ".stats.json")).write_text(json.dumps(stats, indent=1))
        print(json.dumps(stats), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])

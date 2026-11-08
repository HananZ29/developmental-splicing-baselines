"""Verify the four Wen et al. Synapse files for E15 (integrity + column names only; no values are read).

Checks per file: exact byte size vs the Synapse Dataset sumFileSizes (DATA_PROVENANCE_E15.md),
gzip integrity for .gz, line count, and the header line. Prints our MD5 so it can be compared
with the MD5 shown on the Synapse file page. Effect-size / PIP columns are never parsed.

Usage: python -I E15_verify_wen_files.py
"""
import gzip
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "external" / "wen_synapse"
EXPECTED = {
    "mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz": ("syn51183885", 19113819),
    "conditional_top_variants.txt": ("syn51183855", 5888581),
    "leafcutter_clusters_to_genes.txt": ("syn51183854", 3476917),
    "all.introns.tested.tsv": ("syn51183859", 13471910),
}
HEADERED = {"mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz", "leafcutter_clusters_to_genes.txt"}


def main():
    for name, (syn, size) in EXPECTED.items():
        p = D / name
        if not p.exists():
            print(f"{syn} {name}: MISSING")
            continue
        raw = p.read_bytes()
        md5 = hashlib.md5(raw).hexdigest()
        ok = "OK" if len(raw) == size else f"SIZE MISMATCH ({len(raw)} != {size})"
        opener = gzip.open if name.endswith(".gz") else open
        try:
            with opener(p, "rt") as f:
                header = f.readline().rstrip("\n")
                n = 1 + sum(1 for _ in f)
            integ = "readable"
        except Exception as e:  # truncated gzip etc.
            header, n, integ = "", 0, f"UNREADABLE: {e}"
        fields = header.split("\t") if "\t" in header else header.split()
        print(f"{syn} {name}: size {ok}; {integ}; {n} lines; md5 {md5}")
        # Files without a header row would otherwise print a data row (values must not be shown pre-freeze).
        if name in HEADERED:
            print(f"  columns: {fields}")
        else:
            print(f"  no header row; {len(fields)} fields per line (values not shown)")


if __name__ == "__main__":
    main()

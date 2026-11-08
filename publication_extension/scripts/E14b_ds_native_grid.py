"""E14 step 3 (PREREGISTRATION_E14.md §E14d): native DeltaSplice (exon mode) baseline grids.

For each selected variant, run DeltaSplice exon mode at
- baselines PSI in GRID (acceptor SSU = donor SSU = PSI, the authors' convention);
- in `check` mode, also at the exact off-grid baselines CHECK.
Record dacc = acc_alt - acc_ref, ddon = don_alt - don_ref and dpsi = (dacc + ddon) / 2.

Interpolation (PCHIP on logit PSI) and the accuracy test are done in the analysis script, not here.

Exon sites: acceptor = exon_start (+ strand) or exon_end (- strand); donor the other end
(1-based hg19), shifted by the DeltaSplice modal offset.
- check mode: offset 0, measured on training variants during setup.
- fresh mode: the modal offset from the DeltaSplice reference-free tracks of T-fresh (E11 rule), read
  from deltasplice_fresh/.

Modes:
- check: 20 random variants from the E13c TRAINING families (seed 20261106).
- fresh: 40 random T-fresh families x 5 random variants (seed 20261106).

Output: results/E14_comparators/ds_native_<mode>.tsv. Resumable: variants already written are skipped.

Usage: <deltasplice venv python> -E E14b_ds_native_grid.py <check|fresh>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_deltasplice import DeltaSplice  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension" / "results" / "E14_comparators"
SEED = 20261106
GRID = [0.05, 0.2, 0.5, 0.8, 0.95]
CHECK = [0.35, 0.65]
KEY = ["chrom", "pos", "ref", "alt", "strand"]


def sites(v):
    plus = v.strand == "+"
    v = v.copy()
    v["acceptor"] = np.where(plus, v.exon_start_hg19, v.exon_end_hg19).astype(int)
    v["donor"] = np.where(plus, v.exon_end_hg19, v.exon_start_hg19).astype(int)
    return v


def modal_offset_from_tracks(v, tdir):
    """E11 rule: the most common argmax offset (within +-3) of the REF track around the annotated site,
    over the first 150 variants."""
    t = np.load(tdir / "tracks.npz")
    sc = pd.read_csv(tdir / "scores.tsv", sep="\t")
    W = int(t["window"])
    row = {tuple(r): i for i, r in enumerate(sc[KEY].itertuples(index=False))}
    out = {}
    for col, kind in [("acceptor", "acc"), ("donor", "don")]:
        offs = []
        for r in v.head(150).itertuples():
            i = row.get((r.chrom, r.pos, r.ref, r.alt, r.strand))
            c = W + getattr(r, col) - r.pos
            if i is not None and 3 <= c < 2 * W - 2:
                offs.append(int(np.argmax(t[f"{kind}_ref"][i, c - 3:c + 4])) - 3)
        vals, cnt = np.unique(offs, return_counts=True)
        out[col] = int(vals[np.argmax(cnt)])
    return out


def select(mode):
    if mode == "check":
        e13 = pd.read_csv(ROOT / "publication_extension" / "results" / "E13_compass" / "subset_variants.tsv", sep="\t")
        v = e13[e13.split == "train"].drop_duplicates(KEY).sample(20, random_state=SEED)
        return sites(v), {"acceptor": 0, "donor": 0}, GRID + CHECK
    f = pd.read_csv(OUT / "fresh_variants.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    fams = rng.choice(sorted(f.event_id.unique()), 40, replace=False)
    v = f[f.event_id.isin(fams)].drop_duplicates(KEY)
    v = v.groupby("event_id", group_keys=False).apply(lambda g: g.sample(min(len(g), 5), random_state=SEED))
    v = sites(v)
    off = modal_offset_from_tracks(v, OUT / "deltasplice_fresh")
    return v, off, GRID


def main(mode):
    v, off, psis = select(mode)
    path = OUT / f"ds_native_{mode}.tsv"
    done = set()
    if path.exists():
        d0 = pd.read_csv(path, sep="\t")
        done = set(map(tuple, d0[KEY].astype(str).to_numpy()))
    else:
        path.write_text("\t".join(KEY + ["event_id", "acc_site", "don_site", "psi", "dacc", "ddon", "dpsi"]) + "\n")
    ds = DeltaSplice(ROOT / "tools" / "src" / "DeltaSplice")
    print(f"{mode}: {len(v)} variants, offsets {off}, baselines {psis}", flush=True)
    for k, r in enumerate(v.itertuples()):
        if tuple(map(str, (r.chrom, r.pos, r.ref, r.alt, r.strand))) in done:
            continue
        a, d = int(r.acceptor + off["acceptor"]), int(r.donor + off["donor"])
        rows = []
        for p in psis:
            e = ds.exon_delta(dict(chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt, strand=r.strand), a, d, p, p)
            da, dd = e["acc_alt"] - e["acc_ref"], e["don_alt"] - e["don_ref"]
            rows.append([r.chrom, r.pos, r.ref, r.alt, r.strand, r.event_id, a, d, p, da, dd, (da + dd) / 2])
        with path.open("a") as fh:
            for x in rows:
                fh.write("\t".join(map(str, x)) + "\n")
        print(f"{k + 1}/{len(v)}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])

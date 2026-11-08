"""Resumable SpliceAI scoring for any variant table (same scorer and settings as experiments/07).

Usage: python publication_extension/scripts/score_spliceai.py <variants.tsv> <out_dir>
The table needs chrom, pos, ref, alt, strand (REF/ALT on the + strand, hg19). Outputs scores.tsv and
tracks.npz in out_dir; row order of tracks = scores.tsv.
"""
import tensorflow  # noqa: F401  (load TF runtime first on Windows)
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from spliceai_runner import SpliceAIScorer  # noqa: E402

W, BATCH = 500, 16
KEY = ["chrom", "pos", "ref", "alt", "strand"]
DS = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]
TR = ["acc_ref", "acc_alt", "don_ref", "don_alt"]


def main(table, out):
    out = Path(out)
    parts = out / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(table, sep="\t")
    todo = v[KEY].drop_duplicates().sort_values(["chrom", "pos"]).reset_index(drop=True)
    fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))
    scorer = SpliceAIScorer(fa, str(ROOT / "data" / "external" / "spliceai"), window=W)
    t0 = time.time()
    for start in range(0, len(todo), BATCH):
        part = parts / f"part_{start:06d}.npz"
        if part.exists():
            continue
        chunk = todo.iloc[start:start + BATCH]
        res = scorer.score(chunk.to_dict("records"))
        keep = [i for i, r in enumerate(res) if r is not None]
        np.savez_compressed(part, idx=chunk.index.to_numpy()[keep],
                            **{k: np.stack([res[i][k] for i in keep]) for k in TR},
                            **{k: np.array([res[i][k] for i in keep]) for k in DS})
        print(f"{start + len(chunk)}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    ps = [np.load(p) for p in sorted(parts.glob("part_*.npz"))]
    idx = np.concatenate([p["idx"] for p in ps])
    sc = todo.loc[idx].reset_index(drop=True)
    for k in DS:
        sc[k] = np.concatenate([p[k] for p in ps])
    sc["DS_max"] = sc[["DS_AG", "DS_AL", "DS_DG", "DS_DL"]].max(axis=1)
    sc.to_csv(out / "scores.tsv", sep="\t", index=False)
    np.savez_compressed(out / "tracks.npz", window=W, **{k: np.concatenate([p[k] for p in ps]) for k in TR})
    print(f"scored {len(sc)} of {len(todo)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])

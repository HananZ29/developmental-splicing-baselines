"""Score every eligible variant with SpliceAI on hg19 (resumable).

Each unique (chrom, pos, ref, alt, gene strand) is scored once. REF is checked against hg19 and a
mismatch aborts the run. Outputs:
  results/spliceai/scores.tsv   standard DS/DP scores per variant-strand
  results/spliceai/tracks.npz   REF/ALT acceptor and donor tracks (pos-W .. pos+W), row order = scores.tsv

Usage: python experiments/07_run_spliceai.py
"""
import tensorflow  # noqa: F401  (load the TF runtime before other native libraries)
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from spliceai_runner import SpliceAIScorer  # noqa: E402

EXT = ROOT / "data" / "external"
OUT = ROOT / "results" / "spliceai"
PARTS = OUT / "parts"
W, BATCH_VARIANTS = 500, 16
KEY = ["chrom", "pos", "ref", "alt", "strand"]


def main():
    PARTS.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    todo = v[KEY].drop_duplicates().sort_values(["chrom", "pos"]).reset_index(drop=True)
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    scorer = SpliceAIScorer(fa, str(EXT / "spliceai"), window=W)
    t0 = time.time()
    for start in range(0, len(todo), BATCH_VARIANTS):
        part = PARTS / f"part_{start:06d}.npz"
        if part.exists():
            continue
        chunk = todo.iloc[start:start + BATCH_VARIANTS]
        res = scorer.score(chunk.to_dict("records"))
        keep = [i for i, r in enumerate(res) if r is not None]
        np.savez_compressed(
            part, idx=chunk.index.to_numpy()[keep],
            **{k: np.stack([res[i][k] for i in keep]) for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]},
            **{k: np.array([res[i][k] for i in keep]) for k in
               ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]})
        done = start + len(chunk)
        print(f"{done}/{len(todo)} variants  {time.time() - t0:.0f}s", flush=True)

    parts = [np.load(p) for p in sorted(PARTS.glob("part_*.npz"))]
    idx = np.concatenate([p["idx"] for p in parts])
    scores = todo.loc[idx].reset_index(drop=True)
    for k in ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]:
        scores[k] = np.concatenate([p[k] for p in parts])
    scores["DS_max"] = scores[["DS_AG", "DS_AL", "DS_DG", "DS_DL"]].max(axis=1)
    scores.to_csv(OUT / "scores.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "tracks.npz", window=W,
                        **{k: np.concatenate([p[k] for p in parts]) for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]})
    print(f"scored {len(scores)} of {len(todo)} variant-strands; DS_max>=0.2: {(scores['DS_max'] >= 0.2).sum()}")


if __name__ == "__main__":
    main()

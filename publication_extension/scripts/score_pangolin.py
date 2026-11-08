"""Resumable Pangolin scoring with the same interface and outputs as score_spliceai.py.

Pangolin (Zeng & Li 2022, Genome Biol; github.com/tkzeng/Pangolin, weights final.{1-3}.{0,2,4,6}.3.v2)
predicts, for each of four tissues (heart, liver, brain, testis), the probability that a position is
a splice site (acceptor and donor are not distinguished). Following pangolin.py, each tissue is the mean
of its three models, reading output channel [1, 4, 7, 10][tissue].

Input: 10,000 nt of flanking context + (2W+1) nt, in transcript orientation for the given strand,
the same sequence extraction as src/spliceai_runner.py (hg19, REF checked). Tracks are returned in
genomic (+ strand) order: position i corresponds to pos - W + i.

Outputs in out_dir:
- scores.tsv;
- tracks.npz with acc_ref/acc_alt/don_ref/don_alt, all equal to the tissue-mean P(splice) so that the
  E11/E13 site-ratio code can be reused (site identity comes from the exon coordinate), plus per-tissue
  tracks t_ref/t_alt of shape (n, 4, 2W+1).
- Score columns: PG_gain / PG_loss, Pangolin's own +-50 nt max increase / max decrease over tissues;
  unmasked.

Usage: <pangolin venv python> -E score_pangolin.py <variants.tsv> <out_dir> <pangolin_repo_dir>
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[2]
W, BATCH, CONTEXT, DIST = 500, 8, 10000, 50
KEY = ["chrom", "pos", "ref", "alt", "strand"]
TISSUES = ["heart", "liver", "brain", "testis"]
COMP = str.maketrans("ACGTN", "TGCAN")
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot(seq):
    x = np.zeros((4, len(seq)), dtype=np.float32)
    for i, c in enumerate(seq):
        j = IDX.get(c, -1)
        if j >= 0:
            x[j, i] = 1.0
    return x


def load_models(repo):
    sys.path.insert(0, str(repo))
    from pangolin.model import Pangolin, L, W as PW, AR  # noqa: E402
    models = []
    for i in [0, 2, 4, 6]:
        for j in range(1, 4):
            m = Pangolin(L, PW, AR)
            m.load_state_dict(torch.load(Path(repo) / "pangolin" / "models" / f"final.{j}.{i}.3.v2",
                                         map_location="cpu"))
            m.eval()
            models.append(m)
    return models


def seqs_for(fa, v):
    half = CONTEXT // 2 + W
    start = v["pos"] - 1 - half
    if start < 0 or v["pos"] + half > len(fa[v["chrom"]]):
        return None
    s = fa[v["chrom"]][start:v["pos"] + half].seq.upper()
    if s[half] != v["ref"]:
        raise ValueError(f"REF mismatch at {v['chrom']}:{v['pos']}")
    a = s[:half] + v["alt"] + s[half + 1:]
    if v["strand"] == "-":
        s, a = s.translate(COMP)[::-1], a.translate(COMP)[::-1]
    return s, a


@torch.no_grad()
def predict(models, seqs):
    x = torch.from_numpy(np.stack([one_hot(s) for s in seqs]))
    out = np.zeros((len(seqs), 4, 2 * W + 1), dtype=np.float32)
    for t in range(4):
        acc = []
        for m in models[3 * t:3 * t + 3]:
            y = m(x)[:, [1, 4, 7, 10][t], :].numpy()
            acc.append(y)
        out[:, t, :] = np.mean(acc, axis=0)
    return out


def main(table, out, repo):
    out = Path(out)
    parts = out / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(max(1, torch.get_num_threads()))
    v = pd.read_csv(table, sep="\t")
    todo = v[KEY].drop_duplicates().sort_values(["chrom", "pos"]).reset_index(drop=True)
    fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))
    models = load_models(repo)
    t0 = time.time()
    lo, hi = W - DIST, W + DIST + 1
    for start in range(0, len(todo), BATCH):
        part = parts / f"part_{start:06d}.npz"
        if part.exists():
            continue
        chunk = todo.iloc[start:start + BATCH]
        seqs, keep = [], []
        for i, r in enumerate(chunk.to_dict("records")):
            s = seqs_for(fa, r)
            if s is not None:
                seqs += list(s)
                keep.append(i)
        y = predict(models, seqs)
        ref, alt = y[0::2], y[1::2]
        strands = chunk["strand"].to_numpy()[keep]
        for k, sd in enumerate(strands):
            if sd == "-":
                ref[k], alt[k] = ref[k][:, ::-1], alt[k][:, ::-1]
        d = (alt - ref)[:, :, lo:hi]
        np.savez_compressed(part, idx=chunk.index.to_numpy()[keep], t_ref=ref, t_alt=alt,
                            PG_gain=d.max(axis=(1, 2)), PG_loss=-d.min(axis=(1, 2)))
        print(f"{start + len(chunk)}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    ps = [np.load(p) for p in sorted(parts.glob("part_*.npz"))]
    idx = np.concatenate([p["idx"] for p in ps])
    sc = todo.loc[idx].reset_index(drop=True)
    for k in ["PG_gain", "PG_loss"]:
        sc[k] = np.concatenate([p[k] for p in ps])
    sc["PG_max"] = sc[["PG_gain", "PG_loss"]].max(axis=1)
    sc.to_csv(out / "scores.tsv", sep="\t", index=False)
    t_ref = np.concatenate([p["t_ref"] for p in ps])
    t_alt = np.concatenate([p["t_alt"] for p in ps])
    mr, ma = t_ref.mean(axis=1), t_alt.mean(axis=1)
    np.savez_compressed(out / "tracks.npz", window=W, tissues=np.array(TISSUES), t_ref=t_ref, t_alt=t_alt,
                        acc_ref=mr, acc_alt=ma, don_ref=mr, don_alt=ma)
    print(f"scored {len(sc)} of {len(todo)}")


if __name__ == "__main__":
    main(*sys.argv[1:4])

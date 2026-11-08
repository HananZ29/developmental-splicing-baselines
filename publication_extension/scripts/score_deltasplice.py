"""DeltaSplice scoring with the same interface and outputs as score_spliceai.py.

DeltaSplice (Xu et al. 2024, Genome Res 34:1052; github.com/chaolinzhanglab/DeltaSplice @71c2fd16,
MIT) predicts acceptor and donor splice-site usage (SSU) over a 5,000-nt window from 35,000 nt of
sequence. This wrapper reproduces the sequence handling of the authors' pred_deltassu.py:
- seq_start = pos0 - (EL+CL)//2, with the REF base at the centre;
- strand reverse-complement;
- channels 1 = acceptor, 2 = donor;
- the 5-model ensemble mean of DeltaSplice_models/model.ckpt-{0..4}.

Two departures, neither of which changes the model:
1. Runs on CPU. Module.cuda and Tensor.cuda are made no-ops before the package is imported, and the
   weights are loaded with map_location="cpu".
2. Reads hg19 with pyfaidx instead of pyfasta.

Modes:
- tracks (default): reference-free prediction (use_ref=False), as in the authors' no-exon mode.
  tracks.npz has acc_ref/acc_alt/don_ref/don_alt over +-W in genomic (+ strand) order, where
  position i corresponds to pos - W + i (1-based). scores.tsv has DS_* style max |delta| within
  +-50 nt.
- exon_delta(): reference-informed DeltaSSU for a given exon with measured reference usage at its
  acceptor and donor. Use it from a separate script once the coordinate convention is fixed
  (see the setup notes).

Usage: <deltasplice venv python> -E score_deltasplice.py <variants.tsv> <out_dir> <deltasplice_repo>
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta

torch.nn.Module.cuda = lambda self, *a, **k: self
torch.Tensor.cuda = lambda self, *a, **k: self

ROOT = Path(__file__).resolve().parents[2]
W, DIST, BATCH = 500, 50, 4
KEY = ["chrom", "pos", "ref", "alt", "strand"]


class DeltaSplice:
    def __init__(self, repo):
        repo = Path(repo).resolve()
        sys.path.insert(0, str(repo))
        cwd = os.getcwd()
        os.chdir(repo / "deltasplice")  # constant.py uses relative paths
        try:
            import copy
            from deltasplice import constant as C
            self.C = C
            self.models = []
            for p in C.default_model_paths:
                m = copy.deepcopy(C.model)
                m.load_state_dict(torch.load(p, map_location="cpu"))
                m.eval()
                self.models.append(m)
        finally:
            os.chdir(cwd)
        self.EL, self.CL = C.EL, C.CL
        self.fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))

    def _encode(self, seq):
        C = self.C
        return C.IN_MAP[[C.SeqTable[c] for c in seq]][:, :4]

    def _seqs(self, chrom, pos1, ref, alt, strand):
        pos0 = pos1 - 1
        seq_start = pos0 - (self.EL + self.CL) // 2
        seq_end = seq_start + self.EL + self.CL
        n = len(self.fa[chrom])
        s = self.fa[chrom][max(seq_start, 0):min(seq_end, n)].seq.upper()
        s = "N" * max(0, -seq_start) + s + "N" * max(0, seq_end - n)
        if s[pos0 - seq_start] != ref:
            raise ValueError(f"REF mismatch at {chrom}:{pos1}")
        a = s[:pos0 - seq_start] + alt + s[pos0 - seq_start + 1:]
        s = "".join(c if c in "ACGTN" else "N" for c in s)
        a = "".join(c if c in "ACGTN" else "N" for c in a)
        if strand == "-":
            rep = self.C.repdict
            s, a = "".join(rep[c] for c in s)[::-1], "".join(rep[c] for c in a)[::-1]
        return s, a, seq_start

    def _run(self, s, a, refmat=None):
        d = {"X": torch.tensor(self._encode(s))[None], "mutX": torch.tensor(self._encode(a))[None]}
        if refmat is not None:
            d["single_pred_psi"] = torch.tensor(refmat)[None]
        pred = [m.predict(d, use_ref=refmat is not None) for m in self.models]
        r = sum(p["single_pred_psi"] for p in pred) / len(pred)
        mut = sum(p["mutY"] for p in pred) / len(pred)
        return r[0], mut[0]  # (CL, 3) each, transcript orientation

    def tracks(self, v):
        s, a, _ = self._seqs(v["chrom"], v["pos"], v["ref"], v["alt"], v["strand"])
        r, m = self._run(s, a)
        if v["strand"] == "-":
            r, m = r[::-1], m[::-1]
        c = self.CL // 2  # the variant's index in the output window (genomic order)
        r, m = r[c - W:c + W + 1], m[c - W:c + W + 1]
        return {"acc_ref": r[:, 1], "acc_alt": m[:, 1], "don_ref": r[:, 2], "don_alt": m[:, 2]}

    def exon_delta(self, v, acc_pos1, don_pos1, acc_ssu, don_ssu):
        """Reference-informed prediction, as in the authors' exon mode.

        - acc_pos1 / don_pos1 are 1-based genomic positions in the model's own output convention (the
          positions at which the reference-free acc/don tracks peak for this exon).
        - acc_ssu / don_ssu are the measured reference usages (floored at 1e-3 as in the authors' code).

        Returns predicted reference and mutant usage at both sites, averaged over the 5 models.
        """
        s, a, seq_start = self._seqs(v["chrom"], v["pos"], v["ref"], v["alt"], v["strand"])
        L = self.EL + self.CL
        refmat = np.zeros((L, 3))
        ia, idn = acc_pos1 - 1 - seq_start, don_pos1 - 1 - seq_start
        refmat[ia, 1] = max(acc_ssu, 1e-3)
        refmat[idn, 2] = max(don_ssu, 1e-3)
        if v["strand"] == "-":
            refmat = refmat[::-1]
        refmat[:, 0] = 1 - refmat[:, 1:].sum(-1)
        refmat = refmat[self.EL // 2:self.EL // 2 + self.CL].copy()
        r, m = self._run(s, a, refmat)
        if v["strand"] == "-":
            r, m = r[::-1], m[::-1]
        oa, od = ia - self.EL // 2, idn - self.EL // 2  # genomic-order output indices
        return {"acc_ref": float(r[oa, 1]), "acc_alt": float(m[oa, 1]),
                "don_ref": float(r[od, 2]), "don_alt": float(m[od, 2])}


def main(table, out, repo):
    out = Path(out)
    parts = out / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(table, sep="\t")
    todo = v[KEY].drop_duplicates().sort_values(["chrom", "pos"]).reset_index(drop=True)
    ds = DeltaSplice(repo)
    t0 = time.time()
    lo, hi = W - DIST, W + DIST + 1
    TR = ["acc_ref", "acc_alt", "don_ref", "don_alt"]
    for start in range(0, len(todo), BATCH):
        part = parts / f"part_{start:06d}.npz"
        if part.exists():
            continue
        chunk = todo.iloc[start:start + BATCH]
        res = [ds.tracks(r) for r in chunk.to_dict("records")]
        arr = {k: np.stack([x[k] for x in res]).astype(np.float32) for k in TR}
        da = (arr["acc_alt"] - arr["acc_ref"])[:, lo:hi]
        dd = (arr["don_alt"] - arr["don_ref"])[:, lo:hi]
        np.savez_compressed(part, idx=chunk.index.to_numpy(), **arr,
                            DS_AG=da.max(1), DS_AL=-da.min(1), DS_DG=dd.max(1), DS_DL=-dd.min(1))
        print(f"{start + len(chunk)}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    ps = [np.load(p) for p in sorted(parts.glob("part_*.npz"))]
    sc = todo.loc[np.concatenate([p["idx"] for p in ps])].reset_index(drop=True)
    for k in ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]:
        sc[k] = np.concatenate([p[k] for p in ps])
    sc["DS_max"] = sc[["DS_AG", "DS_AL", "DS_DG", "DS_DL"]].max(axis=1)
    sc.to_csv(out / "scores.tsv", sep="\t", index=False)
    np.savez_compressed(out / "tracks.npz", window=W, **{k: np.concatenate([p[k] for p in ps]) for k in TR})
    print(f"scored {len(sc)} of {len(todo)}")


if __name__ == "__main__":
    main(*sys.argv[1:4])

"""SpliceAI (Jaganathan et al. 2019, v1.3.1 weights) scoring of SNVs on hg19, without the
`spliceai` package (its build fails on this Windows/Python 3.12 setup).

For each variant we predict acceptor/donor probabilities for REF and ALT sequence over a window of
+-W bases around the variant (input = 10,000 nt flanking context + 2W+1), in transcript orientation
for the given gene strand, averaging the five ensemble models, as in the reference implementation.
Arrays are returned in genomic (+ strand) order so position i corresponds to pos - W + i.

The standard SpliceAI delta scores (DS_AG, DS_AL, DS_DG, DS_DL within +-50 nt and their positions)
are computed without the annotation masking the reference CLI can optionally apply.
"""
import os

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

import numpy as np

CONTEXT = 10000
COMP = str.maketrans("ACGTN", "TGCAN")
ONEHOT = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot(seq):
    x = np.zeros((len(seq), 4), dtype=np.float32)
    idx = [ONEHOT.get(c, -1) for c in seq]
    for i, j in enumerate(idx):
        if j >= 0:
            x[i, j] = 1.0
    return x


def load_models(model_dir):
    import keras
    return [keras.models.load_model(os.path.join(model_dir, f"spliceai{i}.h5"), compile=False) for i in range(1, 6)]


class SpliceAIScorer:
    def __init__(self, fasta, model_dir, window=500, batch=8):
        self.fa, self.W, self.batch = fasta, window, batch
        self.models = load_models(model_dir)

    def _seqs(self, chrom, pos, ref, alt, strand):
        half = CONTEXT // 2 + self.W
        start = pos - 1 - half
        if start < 0 or pos + half > len(self.fa[chrom]):
            return None
        s = self.fa[chrom][start:pos + half].seq.upper()
        if s[half] != ref:
            raise ValueError(f"REF mismatch at {chrom}:{pos}: genome={s[half]} given={ref}")
        a = s[:half] + alt + s[half + 1:]
        if strand == "-":
            s, a = s.translate(COMP)[::-1], a.translate(COMP)[::-1]
        return s, a

    def _predict(self, seqs):
        x = np.stack([one_hot(s) for s in seqs])
        # direct calls avoid Keras predict() retracing on varying batch sizes
        y = np.mean([np.concatenate([m(x[i:i + self.batch], training=False).numpy()
                                     for i in range(0, len(x), self.batch)]) for m in self.models], axis=0)
        return y  # (n, 2W+1, 3): neither, acceptor, donor

    def score(self, variants):
        """variants: list of dicts with chrom, pos, ref, alt, strand. Returns list of result dicts
        (None for variants too close to contig ends)."""
        seqs, keep = [], []
        for i, v in enumerate(variants):
            s = self._seqs(v["chrom"], v["pos"], v["ref"], v["alt"], v["strand"])
            if s is not None:
                seqs += list(s)
                keep.append(i)
        out = [None] * len(variants)
        if not seqs:
            return out
        y = self._predict(seqs)
        for k, i in enumerate(keep):
            r, a = y[2 * k], y[2 * k + 1]
            if variants[i]["strand"] == "-":
                r, a = r[::-1], a[::-1]
            out[i] = summarize(r, a, self.W)
        return out


def summarize(ref, alt, W, dist=50):
    """Keep full acceptor/donor tracks plus standard +-50 nt delta scores (genomic orientation)."""
    acc_r, don_r, acc_a, don_a = ref[:, 1], ref[:, 2], alt[:, 1], alt[:, 2]
    lo, hi = W - dist, W + dist + 1
    d_acc, d_don = (acc_a - acc_r)[lo:hi], (don_a - don_r)[lo:hi]
    rel = np.arange(-dist, dist + 1)
    return {
        "acc_ref": acc_r.astype(np.float32), "acc_alt": acc_a.astype(np.float32),
        "don_ref": don_r.astype(np.float32), "don_alt": don_a.astype(np.float32),
        "DS_AG": float(d_acc.max()), "DP_AG": int(rel[d_acc.argmax()]),
        "DS_AL": float(-d_acc.min()), "DP_AL": int(rel[d_acc.argmin()]),
        "DS_DG": float(d_don.max()), "DP_DG": int(rel[d_don.argmax()]),
        "DS_DL": float(-d_don.min()), "DP_DL": int(rel[d_don.argmin()]),
    }

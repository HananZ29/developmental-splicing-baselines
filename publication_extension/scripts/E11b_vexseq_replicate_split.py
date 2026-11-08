"""E11b (post-hoc robustness check, deviation D6): replicate-split version of E11.

Concern: in E11 the measured dPSI = PSI_variant - PSI_reference uses the SAME replicates as the baseline
PSI_reference that the matched predictor M uses. Noise in PSI_reference then enters both the predictor and
(with opposite sign) the target, which could favour M over the swapped (S) or constant (C) baseline
artefactually (regression to the mean).

Fix: for each rotation r in {1,2,3}, the predictor baseline is the reference-exon PSI from replicate r only;
the target dPSI is computed from the other two replicates only (variant and reference both). Contrasts are
averaged over the three rotations; exon-cluster bootstrap (B = 2000) uses the same exon draws for all
rotations. SpliceAI factor a and the permutation for null N are taken unchanged from E11.
Output: publication_extension/results/E11_vexseq/replicate_split.json
"""
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / "data" / "external" / "vexseq" / "geo"
OUT = ROOT / "publication_extension" / "results" / "E11_vexseq"
CELLS = ["K562", "HepG2"]
SEED, N_BOOT = 20261105, 2000


def reps(cell):
    fs = [f for f in sorted(glob.glob(str(V / f"GSM*_{cell}_rep*_PSI.tsv.gz"))) if "KD" not in f and "barcode" not in f]
    return [pd.read_csv(f, sep="\t").set_index("variant").PSI / 100 for f in fs]


def prop(p, a):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return p * a / (p * a + 1 - p) - p


def pear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    return pearsonr(x, y).statistic if len(x) > 2 and np.std(x) > 0 and np.std(y) > 0 else np.nan


def main():
    v = pd.read_csv(OUT / "predictions.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    perm_a = np.random.default_rng(20261104).permutation(v.a.to_numpy())  # same as E11 (first draw of that rng)
    a = v.a.to_numpy()
    ref_id = v.chrom + "_" + v.ctrl_pos.astype(int).astype(str) + "_Control_Original"
    R = {c: reps(c) for c in CELLS}
    rot = []
    for r in range(3):
        d = {}
        for c in CELLS:
            o = [i for i in range(3) if i != r]
            base = R[c][r].reindex(ref_id).to_numpy()
            var_m = np.nanmean([R[c][i].reindex(v.variant).to_numpy() for i in o], axis=0)
            ref_m = np.nanmean([R[c][i].reindex(ref_id).to_numpy() for i in o], axis=0)
            d[c] = {"base": base, "y": var_m - ref_m}
        for c in CELLS:
            oc = [x for x in CELLS if x != c][0]
            d[c]["M"] = prop(d[c]["base"], a)
            d[c]["S"] = prop(d[oc]["base"], a)
            d[c]["C"] = prop(np.full(len(v), np.nanmean(pd.Series(d[c]["base"]).groupby(v.exon).first())), a)
            d[c]["N"] = prop(d[c]["base"], perm_a)
        rot.append(d)

    ex = v.exon.to_numpy()
    ue = np.unique(ex)
    wi = {e: np.flatnonzero(ex == e) for e in ue}
    big = set(v.loc[(v.K562_base - v.HepG2_base).abs() >= 0.1, "exon"])

    def stats(ix):
        out = {k: [] for k in ["M", "S", "C", "N", "between_M", "between_N"]}
        for d in rot:
            y = np.concatenate([d[c]["y"][ix] for c in CELLS])
            for k in ["M", "S", "C", "N"]:
                out[k].append(pear(np.concatenate([d[c][k][ix] for c in CELLS]), y))
            md = d["K562"]["y"][ix] - d["HepG2"]["y"][ix]
            out["between_M"].append(pear(d["K562"]["M"][ix] - d["HepG2"]["M"][ix], md))
            out["between_N"].append(pear(d["K562"]["N"][ix] - d["HepG2"]["N"][ix], md))
        return {k: float(np.nanmean(x)) for k, x in out.items()}

    def run(exons):
        pt = stats(np.concatenate([wi[e] for e in exons]))
        B = pd.DataFrame([stats(np.concatenate([wi[e] for e in rng.choice(exons, len(exons))])) for _ in range(N_BOOT)])

        def con(x, y):
            dd = (B[x] - B[y]).dropna()
            return {"estimate": pt[x] - pt[y], "ci": [float(dd.quantile(.025)), float(dd.quantile(.975))],
                    "p_boot": float(min(1, 2 * min((dd <= 0).mean(), (dd >= 0).mean())))}
        return {"point": pt, "M - S": con("M", "S"), "M - C": con("M", "C"), "between M - N": con("between_M", "between_N")}

    rep = {"design": "baseline from replicate r, target dPSI from the other two replicates; averaged over r = 1..3",
           "n_variants": len(v), "n_exons": len(ue),
           "all_exons": run(ue), "baseline_diff_ge_0.1_exons": run(np.array(sorted(big)))}
    (OUT / "replicate_split.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()

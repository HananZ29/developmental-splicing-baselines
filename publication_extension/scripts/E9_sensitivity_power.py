"""E9: sensitivity of predictions to fetal vs adult baselines, near-zero sign validity, and the power of
the direction/rank benchmark to detect a baseline that is, by construction, correct.

A. Per-variant fetal (Leung) vs adult on the 147 Wen sQTL positives (v2 predictions):
   - numerically zero predictions (|pred| < 1e-12) and negligible ones (< 1e-3)
   - sign flips counted only where both predictions exceed a tolerance tau
   - magnitude ratio |fetal| / |adult|, ranks
B. Direction accuracy with ties for |pred| < tau (tau = 0 as in v2, 1e-12, 1e-6, 1e-3), fetal - adult,
   with the same 1-Mb block bootstrap (B = 2000, seed 11).
C. Oracle power simulation (S = 500, seed 20261104): the truth is generated from the FETAL
   prediction, y = d_fetal + sigma * N(0,1), with sigma chosen so the fetal direction accuracy matches
   the observed 0.70. For each simulation: fetal - adult in direction accuracy (exact McNemar on
   discordant pairs, one-sided alpha 0.05), Spearman (block bootstrap, 200 draws, one-sided 95%) and,
   for contrast, Pearson on the raw magnitudes (bootstrap) -- what a magnitude-measuring benchmark
   could see.
Output: publication_extension/results/E9_sensitivity/
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import binomtest, norm, spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension" / "results" / "E9_sensitivity"
TAUS = [0.0, 1e-12, 1e-6, 1e-3]


def dir_acc(p, z, tau=0.0):
    s = np.where(np.abs(p) <= tau, 0, np.sign(p))
    return float(np.mean(np.where(s == 0, 0.5, (s == np.sign(z)).astype(float))))


def blocks(df):
    b = (df.chrom + ":" + (df.pos // 1_000_000).astype(str)).to_numpy()
    ub = np.unique(b)
    return ub, {u: np.flatnonzero(b == u) for u in ub}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = pd.read_csv(ROOT / "results" / "experiment_v2" / "predictions.tsv", sep="\t")
    w = p[(p.arm == "sQTL") & (p.role == "positive")].reset_index(drop=True)
    f, a, z = w.pred_fetal.to_numpy(), w.pred_adult.to_numpy(), w.z_obs.to_numpy()
    res = {"n": len(w)}

    # A
    A = {}
    for tau in [1e-12, 1e-6, 1e-3, 1e-2]:
        both = (np.abs(f) > tau) & (np.abs(a) > tau)
        A[f"tau={tau:g}"] = {"fetal_|pred|<=tau": int((np.abs(f) <= tau).sum()), "adult_|pred|<=tau": int((np.abs(a) <= tau).sum()),
                             "both_above_tau": int(both.sum()),
                             "sign_flips_among_both_above_tau": int((np.sign(f[both]) != np.sign(a[both])).sum())}
    nz = (np.abs(f) > 1e-12) & (np.abs(a) > 1e-12)
    ratio = np.abs(f[nz]) / np.abs(a[nz])
    A["magnitude_ratio_fetal_over_adult_quantiles(nonzero)"] = {str(q): round(float(np.quantile(ratio, q)), 4) for q in [.05, .25, .5, .75, .95]}
    A["frac_ratio_outside_[0.5,2]"] = round(float(((ratio < 0.5) | (ratio > 2)).mean()), 3)
    A["frac_ratio_outside_[0.1,10]"] = round(float(((ratio < 0.1) | (ratio > 10)).mean()), 3)
    A["spearman_fetal_vs_adult_pred"] = round(float(spearmanr(f, a).statistic), 4)
    A["spearman_abs_fetal_vs_abs_adult"] = round(float(spearmanr(np.abs(f), np.abs(a)).statistic), 4)
    A["top20_overlap_by_abs_pred"] = int(len(set(np.argsort(-np.abs(f))[:20]) & set(np.argsort(-np.abs(a))[:20])))
    res["A_prediction_sensitivity"] = A

    # B
    rng = np.random.default_rng(11)
    ub, idx = blocks(w)
    draws = [np.concatenate([idx[u] for u in rng.choice(ub, len(ub))]) for _ in range(2000)]
    B = {}
    for tau in TAUS:
        est = dir_acc(f, z, tau) - dir_acc(a, z, tau)
        d = np.array([dir_acc(f[s], z[s], tau) - dir_acc(a[s], z[s], tau) for s in draws])
        B[f"tau={tau:g}"] = {"fetal": dir_acc(f, z, tau), "adult": dir_acc(a, z, tau), "fetal_minus_adult": est,
                             "ci": [float(np.quantile(d, .025)), float(np.quantile(d, .975))],
                             "p_boot": float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean())))}
    res["B_direction_accuracy_with_tolerance"] = B

    # C oracle power
    rs = np.random.default_rng(20261104)
    target = 0.70

    def expected_acc(sig):
        return float(np.mean(norm.cdf(np.abs(f) / sig)))
    lo_f = 1e-12
    sigma = brentq(lambda s: expected_acc(s) - target, 1e-9, 10)
    S = 500
    det_dir = det_rho = det_pear = 0
    diffs = {"dir": [], "rho": [], "pearson": []}
    for _ in range(S):
        y = f + sigma * rs.standard_normal(len(f))
        cf, ca = np.sign(f) == np.sign(y), np.sign(a) == np.sign(y)
        b_, c_ = int((cf & ~ca).sum()), int((~cf & ca).sum())
        pd_ = binomtest(b_, b_ + c_, 0.5, alternative="greater").pvalue if b_ + c_ else 1.0
        det_dir += pd_ < 0.05
        diffs["dir"].append((cf.mean() - ca.mean()))
        bd_r, bd_p = [], []
        for s in draws[:200]:
            bd_r.append(spearmanr(f[s], y[s]).statistic - spearmanr(a[s], y[s]).statistic)
            bd_p.append(np.corrcoef(f[s], y[s])[0, 1] - np.corrcoef(a[s], y[s])[0, 1])
        det_rho += np.quantile(bd_r, 0.05) > 0
        det_pear += np.quantile(bd_p, 0.05) > 0
        diffs["rho"].append(spearmanr(f, y).statistic - spearmanr(a, y).statistic)
        diffs["pearson"].append(np.corrcoef(f, y)[0, 1] - np.corrcoef(a, y)[0, 1])
    res["C_oracle_power"] = {
        "truth": "y = pred_fetal + sigma*N(0,1) (fetal baseline is correct by construction)",
        "sigma": sigma, "expected_fetal_dir_acc": round(expected_acc(sigma), 4), "S": S,
        "power_direction_accuracy_mcnemar": det_dir / S,
        "power_spearman_bootstrap": det_rho / S,
        "power_pearson_on_magnitudes_bootstrap": det_pear / S,
        "mean_fetal_minus_adult": {k: round(float(np.mean(v)), 4) for k, v in diffs.items()},
        "max_possible_dir_diff": round(float((np.sign(f) != np.sign(a)).mean()), 4),
    }
    (OUT / "summary.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()

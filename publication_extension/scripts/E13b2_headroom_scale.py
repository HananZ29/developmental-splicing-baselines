"""E13b2 (secondary, descriptive; added before any model comparison): put the E13b headroom on an
absolute scale. Point estimates only (no bootstrap).

Reports, for COMPASS single SNVs passing QC (same definitions as E13b):
  var(y)            variance of measured dPSI (target replicate) pooled over cells
  var(D*)           reliable variance of between-cell dPSI differences = rel * var(D_obs)
  unexplained       H * var(D*) = reliable between-cell variance NOT explained by baseline-only logit additivity
  as fractions of var(y), and the RMS size of each component in PSI units.
Also the same with D_L averaged over both input replicates (reduces noise in the baseline-only predictor,
so H is less inflated by predictor noise).
Output: publication_extension/results/E13_compass/headroom_scale.json
"""
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension" / "results" / "E13_compass"
spec = importlib.util.spec_from_file_location("e13b", ROOT / "publication_extension" / "scripts" / "E13b_compass_context_decomposition.py")
e13b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e13b)
CELLS, clip, pear = e13b.CELLS, e13b.clip, e13b.pear


def main():
    d = pd.read_csv(ROOT / "data" / "external" / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt = d[d.snp == "none"].drop_duplicates("event_id").set_index("event_id")
    s = d[(d.seq_type == "single") & d.event_id.isin(wt.index)].reset_index(drop=True)
    W = wt.loc[s.event_id]
    K = len(CELLS)

    def arr(df, col):
        return df[col].to_numpy(float)
    psi = {r: np.stack([arr(s, f"{c}_rep{r}_psi_raw") for c in CELLS], 1) for r in (1, 2)}
    wps = {r: np.stack([arr(W, f"{c}_rep{r}_psi_raw") for c in CELLS], 1) for r in (1, 2)}
    q = np.ones((len(s), K), bool)
    for j, c in enumerate(CELLS):
        for r in (1, 2):
            q[:, j] &= (arr(s, f"{c}_rep{r}_included") + arr(s, f"{c}_rep{r}_excluded")) >= 20
            q[:, j] &= (arr(W, f"{c}_rep{r}_included") + arr(W, f"{c}_rep{r}_excluded")) >= 20
    pooled = np.stack([arr(s, f"{c}_pooled_psi_raw") for c in CELLS], 1)
    wpool = np.stack([arr(W, f"{c}_pooled_psi_raw") for c in CELLS], 1)
    dlog_pool = np.where(q, logit(clip(pooled)) - logit(clip(wpool)), np.nan)
    dps = {r: np.where(q, psi[r] - wps[r], np.nan) for r in (1, 2)}
    out = {}
    Dobs, Din, DL, DLavg, Y = [], [], [], [], []
    for a, b in itertools.combinations(range(K), 2):
        rest = [k for k in range(K) if k not in (a, b)]
        use = q[:, a] & q[:, b] & (q[:, rest].sum(1) >= 1)
        dbar = np.nanmean(dlog_pool[:, rest], 1)
        Lr = {}
        for r in (1, 2):
            pa, pb = clip(wps[r][:, a]), clip(wps[r][:, b])
            Lr[r] = (expit(logit(pa) + dbar) - pa) - (expit(logit(pb) + dbar) - pb)
        for rin, rtg in [(1, 2), (2, 1)]:
            Dobs.append((dps[rtg][:, a] - dps[rtg][:, b])[use])
            Din.append((dps[rin][:, a] - dps[rin][:, b])[use])
            DL.append(Lr[rin][use])
            # averaged-baseline predictor: average of the two input replicates' baselines is NOT independent of
            # the target replicate, so use it only as an approximate lower bound on H (flagged)
            DLavg.append(((Lr[1] + Lr[2]) / 2)[use])
    for j in range(K):
        Y.append(np.concatenate([dps[2][q[:, j], j], dps[1][q[:, j], j]]))
    Dobs, Din, DL, DLavg, Y = map(np.concatenate, (Dobs, Din, DL, DLavg, Y))
    ok = np.isfinite(Dobs) & np.isfinite(Din) & np.isfinite(DL)
    Dobs, Din, DL, DLavg = Dobs[ok], Din[ok], DL[ok], DLavg[ok]
    rel = pear(Din, Dobs)
    r = pear(DL, Dobs)
    H = 1 - r ** 2 / rel
    vy, vd = float(np.nanvar(Y)), float(np.var(Dobs))
    vstar = rel * vd
    out = {"var_dPSI_y": vy, "rms_dPSI_y": vy ** .5,
           "var_Dobs_between_cell_diff": vd, "reliability": rel, "var_Dstar_reliable": vstar, "rms_Dstar": vstar ** .5,
           "H": H, "unexplained_reliable_var": H * vstar, "rms_unexplained": (H * vstar) ** .5,
           "unexplained_as_frac_of_var_dPSI": H * vstar / vy,
           "explained_by_baseline_only_as_frac_of_var_dPSI": (1 - H) * vstar / vy,
           "mean_abs_Dobs": float(np.mean(np.abs(Dobs))), "frac_pairs_abs_Dobs_ge_0.1": float(np.mean(np.abs(Dobs) >= 0.1)),
           "H_with_replicate_averaged_baseline_(shares_noise_with_target;_approx_lower_bound)": 1 - pear(DLavg, Dobs) ** 2 / rel}
    (OUT / "headroom_scale.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

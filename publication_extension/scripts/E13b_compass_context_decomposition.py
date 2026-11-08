"""E13b (pre-registered): how much of the between-cell-line difference in measured variant effects is
explained solely by the cell line's own reference PSI under logit additivity, and how much is left for a
learned context interaction?

COMPASS (GSE307247) single SNVs with a WT row; cells HEK, HeLa, K562, MCF7, HMC3; QC >= 20 reads in rep1
and rep2 for variant and WT. Baseline from replicate r_in, target from r_tg (both orders averaged).
Exon-family bootstrap B = 1000 (seed 20261105).
Output: publication_extension/results/E13_compass/context_decomposition.json
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension" / "results" / "E13_compass"
CELLS = ["HEK", "HeLa", "K562", "MCF7", "HMC3"]
MIN_READS, SEED, NB = 20, 20261105, 1000


def clip(x):
    return np.clip(x, 0.01, 0.99)


def pear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10:
        return np.nan
    x, y = x[ok], y[ok]
    return float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / "data" / "external" / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt = d[d.snp == "none"].drop_duplicates("event_id").set_index("event_id")
    s = d[(d.seq_type == "single") & d.event_id.isin(wt.index)].reset_index(drop=True)
    n, K = len(s), len(CELLS)
    W = wt.loc[s.event_id]

    def arr(df, col):
        return df[col].to_numpy(float)
    psi = {r: np.stack([arr(s, f"{c}_rep{r}_psi_raw") for c in CELLS], 1) for r in (1, 2)}
    wps = {r: np.stack([arr(W, f"{c}_rep{r}_psi_raw") for c in CELLS], 1) for r in (1, 2)}
    q = np.ones((n, K), bool)
    for j, c in enumerate(CELLS):
        for r in (1, 2):
            q[:, j] &= (arr(s, f"{c}_rep{r}_included") + arr(s, f"{c}_rep{r}_excluded")) >= MIN_READS
            q[:, j] &= (arr(W, f"{c}_rep{r}_included") + arr(W, f"{c}_rep{r}_excluded")) >= MIN_READS
    pooled = np.stack([arr(s, f"{c}_pooled_psi_raw") for c in CELLS], 1)
    wpool = np.stack([arr(W, f"{c}_pooled_psi_raw") for c in CELLS], 1)
    dpsi_pool = np.where(q, pooled - wpool, np.nan)
    dlog_pool = np.where(q, logit(clip(pooled)) - logit(clip(wpool)), np.nan)
    wpool_q = np.where(q, wpool, np.nan)
    dps = {r: np.where(q, psi[r] - wps[r], np.nan) for r in (1, 2)}
    dlg = {r: np.where(q, logit(clip(psi[r])) - logit(clip(wps[r])), np.nan) for r in (1, 2)}
    ev = s.event_id.to_numpy()

    audit = {"singles_with_wt": n, "events": int(pd.Series(ev).nunique()),
             "qc_pass_per_cell": {c: int(q[:, j].sum()) for j, c in enumerate(CELLS)},
             "qc_pass_all5": int(q.all(1).sum()), "qc_pass_ge2": int((q.sum(1) >= 2).sum()),
             "events_with_any_qc_pass": int(pd.Series(ev[q.any(1)]).nunique()),
             "wt_psi_median_qc": {c: round(float(np.nanmedian(wpool_q[:, j])), 3) for j, c in enumerate(CELLS)}}

    # ---- leave-one-cell-out predictors (row = variant x cell x split)
    rows = []
    for j in range(K):
        others = [k for k in range(K) if k != j]
        nO = q[:, others].sum(1)
        use = q[:, j] & (nO >= 2)
        dbar = np.nanmean(dlog_pool[:, others], 1)
        tbar = np.nanmean(dpsi_pool[:, others], 1)
        pother = np.nanmean(wpool_q[:, others], 1)
        for rin, rtg in [(1, 2), (2, 1)]:
            p = clip(wps[rin][:, j])
            L = expit(logit(p) + dbar) - p
            Lo = expit(logit(clip(pother)) + dbar) - clip(pother)
            rows.append(pd.DataFrame({"i": np.flatnonzero(use), "cell": CELLS[j], "split": f"{rin}->{rtg}",
                                      "y": dps[rtg][use, j], "own_rep": dps[rin][use, j], "T": tbar[use], "L": L[use], "L_other": Lo[use],
                                      "wt_in": p[use]}))
    P = pd.concat(rows, ignore_index=True)
    P["event"] = ev[P.i]

    # ---- pairwise between-context decomposition
    prow = []
    for a, b in itertools.combinations(range(K), 2):
        rest = [k for k in range(K) if k not in (a, b)]
        use = q[:, a] & q[:, b] & (q[:, rest].sum(1) >= 1)
        dbar = np.nanmean(dlog_pool[:, rest], 1)
        for rin, rtg in [(1, 2), (2, 1)]:
            pa, pb = clip(wps[rin][:, a]), clip(wps[rin][:, b])
            La, Lb = expit(logit(pa) + dbar) - pa, expit(logit(pb) + dbar) - pb
            prow.append(pd.DataFrame({"i": np.flatnonzero(use), "pair": f"{CELLS[a]}-{CELLS[b]}", "split": f"{rin}->{rtg}",
                                      "D_obs": (dps[rtg][:, a] - dps[rtg][:, b])[use], "D_in": (dps[rin][:, a] - dps[rin][:, b])[use],
                                      "D_L": (La - Lb)[use], "wt_gap": np.abs(wpool[:, a] - wpool[:, b])[use]}))
    Q = pd.concat(prow, ignore_index=True)
    Q["event"] = ev[Q.i]

    # ---- context-invariance of the logit effect
    inv = {"within_cell_replicate_r_dlogit": {c: round(pear(dlg[1][:, j], dlg[2][:, j]), 3) for j, c in enumerate(CELLS)},
           "between_cell_r_dlogit_pooled": {f"{CELLS[a]}-{CELLS[b]}": round(pear(dlog_pool[:, a], dlog_pool[:, b]), 3) for a, b in itertools.combinations(range(K), 2)},
           "between_cell_r_dlogit_cross_replicate": {f"{CELLS[a]}-{CELLS[b]}": round(np.mean([pear(dlg[1][:, a], dlg[2][:, b]), pear(dlg[2][:, a], dlg[1][:, b])]), 3)
                                                     for a, b in itertools.combinations(range(K), 2)}}
    w = np.array(list(inv["within_cell_replicate_r_dlogit"].values()))
    inv["mean_within"] = round(float(np.nanmean(w)), 3)
    bc = np.array(list(inv["between_cell_r_dlogit_cross_replicate"].values()))
    inv["mean_between_cross_replicate"] = round(float(np.nanmean(bc)), 3)
    # disattenuated between-cell correlation of the logit effect (1 = perfectly context-invariant)
    inv["disattenuated_between_cell_r"] = {f"{CELLS[a]}-{CELLS[b]}": round(inv["between_cell_r_dlogit_cross_replicate"][f"{CELLS[a]}-{CELLS[b]}"] / np.sqrt(w[a] * w[b]), 3)
                                           for a, b in itertools.combinations(range(K), 2)}

    def stats(Pi, Qi):
        out = {k: pear(Pi[k].to_numpy(), Pi.y.to_numpy()) for k in ["T", "L", "L_other", "own_rep"]}
        out |= {f"mae_{k}": float(np.nanmean(np.abs(Pi[k] - Pi.y))) for k in ["T", "L", "L_other"]}
        r_dl = pear(Qi.D_L.to_numpy(), Qi.D_obs.to_numpy())
        rel = pear(Qi.D_in.to_numpy(), Qi.D_obs.to_numpy())
        out["pair_r_DL_Dobs"], out["pair_reliability"] = r_dl, rel
        out["S_share_explained"] = r_dl ** 2 / rel if rel and rel > 0 else np.nan
        out["H_headroom"] = 1 - out["S_share_explained"]
        g = Qi[Qi.wt_gap >= 0.1]
        r2, rl2 = pear(g.D_L.to_numpy(), g.D_obs.to_numpy()), pear(g.D_in.to_numpy(), g.D_obs.to_numpy())
        out["gap>=0.1_pair_r_DL_Dobs"], out["gap>=0.1_reliability"] = r2, rl2
        out["gap>=0.1_S"] = r2 ** 2 / rl2 if rl2 and rl2 > 0 else np.nan
        out["gap>=0.1_H"] = 1 - out["gap>=0.1_S"]
        return out
    pt = stats(P, Q)
    ue = np.unique(ev[P.i])
    pix = {e: g for e, g in P.groupby("event").indices.items()}
    qix = {e: g for e, g in Q.groupby("event").indices.items()}
    rng = np.random.default_rng(SEED)
    draws = []
    for _ in range(NB):
        e = rng.choice(ue, len(ue))
        pi = np.concatenate([pix[x] for x in e if x in pix])
        qi = np.concatenate([qix[x] for x in e if x in qix])
        draws.append(stats(P.iloc[pi], Q.iloc[qi]))
    B = pd.DataFrame(draws)

    def ci(k):
        return [round(float(B[k].quantile(.025)), 4), round(float(B[k].quantile(.975)), 4)]

    def con(a, b):
        dd = (B[a] - B[b]).dropna()
        return {"estimate": round(pt[a] - pt[b], 4), "ci": [round(float(dd.quantile(.025)), 4), round(float(dd.quantile(.975)), 4)],
                "p_boot": float(min(1, 2 * min((dd <= 0).mean(), (dd >= 0).mean())))}
    rep = {"audit": audit, "context_invariance_of_dlogit": inv,
           "n_rows_loco": len(P), "n_rows_pairs": len(Q), "n_events_bootstrap": len(ue),
           "point": {k: round(v, 4) for k, v in pt.items()}, "ci": {k: ci(k) for k in pt},
           "primary": {"L - T": con("L", "T"), "L - L_other": con("L", "L_other")},
           "gate_G-C": {"H_lower95": ci("H_headroom")[0], "proceed_with_C_if_H_lower95>0.25": bool(ci("H_headroom")[0] > 0.25)},
           "per_cell_point": {c: {k: round(pear(P[(P.cell == c)][k].to_numpy(), P[P.cell == c].y.to_numpy()), 4) for k in ["T", "L", "L_other", "own_rep"]} for c in CELLS}}
    (OUT / "context_decomposition.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()

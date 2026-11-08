"""E11 (pre-registered): does a context-matched measured baseline improve quantitative ΔPSI prediction?
Vex-seq, K562 and HepG2 (not brain, not developmental). Two-isoform propagation as in E5.

Predictors per cell line c: M matched p_c; S swapped (other cell's p); C constant (mean p_c); R raw
SpliceAI site change; N null (p_c with permuted factor a). Primary: M-S and M-C (Pearson, pooled over
cells) and the between-context test cor(pred_K - pred_H, meas_K - meas_H) vs null N.
Exon-cluster bootstrap B = 2000. Bonferroni alpha = 0.0167 (3 primary tests).
Output: publication_extension/results/E11_vexseq/
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import expit, logit
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "publication_extension" / "results" / "E10_vexseq"
OUT = ROOT / "publication_extension" / "results" / "E11_vexseq"
EPS, SEED, N_BOOT = 0.01, 20261104, 2000
KEY = ["chrom", "pos", "ref", "alt", "strand"]
CELLS = ["K562", "HepG2"]


def site(v, tr, W, col, kind):
    o = (v[col] - v.pos).to_numpy()
    ok = np.abs(o) <= W
    rr, ra = np.full(len(v), np.nan), np.full(len(v), np.nan)
    i = np.flatnonzero(ok)
    rows = v.track_row.to_numpy()[i]
    rr[i], ra[i] = tr[f"{kind}_ref"][rows, W + o[i]], tr[f"{kind}_alt"][rows, W + o[i]]
    return rr, ra


def modal_offset(v, tr, W, col, kind):
    offs = []
    for r in v.head(150).itertuples():
        c = W + getattr(r, col) - r.pos
        if 3 <= c < 2 * W - 2:
            offs.append(int(np.argmax(tr[f"{kind}_ref"][r.track_row, c - 3:c + 4])) - 3)
    vals, cnt = np.unique(offs, return_counts=True)
    return int(vals[np.argmax(cnt)]), round(float(np.mean(np.array(offs) == vals[np.argmax(cnt)])), 3)


def prop(p, a):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return p * a / (p * a + 1 - p) - p


def pear(x, y):
    return pearsonr(x, y).statistic if np.std(x) > 0 and np.std(y) > 0 else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    v = pd.read_csv(D / "vexseq_measured.tsv", sep="\t")
    v["strand"] = v.strand.astype(str)
    sc = pd.read_csv(D / "spliceai" / "scores.tsv", sep="\t")
    sc["strand"] = sc.strand.astype(str)
    sc["track_row"] = np.arange(len(sc))
    v = v.merge(sc[KEY + ["track_row"]], on=KEY, how="inner")
    t = np.load(D / "spliceai" / "tracks.npz")
    W = int(t["window"])
    tr = {k: t[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    v["acceptor"] = np.where(v.strand == "+", v["Exon up"], v["Exon down"]).astype(int)
    v["donor"] = np.where(v.strand == "+", v["Exon down"], v["Exon up"]).astype(int)
    conv = {}
    for col, kind in [("acceptor", "acc"), ("donor", "don")]:
        off, frac = modal_offset(v, tr, W, col, kind)
        conv[col] = {"modal_offset": off, "frac_at_modal": frac}
        v[col] = v[col] + off
    ar, aa = site(v, tr, W, "acceptor", "acc")
    dr, da = site(v, tr, W, "donor", "don")
    ok = np.isfinite(ar) & np.isfinite(dr)
    v, ar, aa, dr, da = v[ok].reset_index(drop=True), ar[ok], aa[ok], dr[ok], da[ok]
    a = (aa + EPS) / (ar + EPS) * (da + EPS) / (dr + EPS)
    v["a"] = a
    v["R"] = (aa - ar) + (da - dr)
    ex_id = v[["chrom", "Exon up", "Exon down"]].astype(str).agg(":".join, axis=1)
    v["exon"] = ex_id
    perm_a = rng.permutation(a)
    for c in CELLS:
        o = [x for x in CELLS if x != c][0]
        v[f"M_{c}"] = prop(v[f"{c}_base"].to_numpy(), a)
        v[f"S_{c}"] = prop(v[f"{o}_base"].to_numpy(), a)
        v[f"C_{c}"] = prop(np.full(len(v), v.drop_duplicates("exon")[f"{c}_base"].mean()), a)
        v[f"N_{c}"] = prop(v[f"{c}_base"].to_numpy(), perm_a)
    # calibrated k (train exons) on the matched predictor
    exons = v.exon.unique()
    train = set(rng.choice(exons, len(exons) // 2, replace=False))
    v["split"] = np.where(v.exon.isin(train), "train", "test")
    tm = (v.split == "train").to_numpy()
    la = np.log(a)

    def loss(k):
        return sum(np.mean((expit(logit(np.clip(v[f"{c}_base"].to_numpy()[tm], 1e-4, 1 - 1e-4)) + k * la[tm])
                            - np.clip(v[f"{c}_base"].to_numpy()[tm], 1e-4, 1 - 1e-4) - v[f"{c}_dpsi"].to_numpy()[tm]) ** 2) for c in CELLS)
    k = minimize_scalar(loss, bounds=(0, 10), method="bounded").x
    for c in CELLS:
        p = np.clip(v[f"{c}_base"].to_numpy(), 1e-4, 1 - 1e-4)
        v[f"K_{c}"] = expit(logit(p) + k * la) - p
    v.to_csv(OUT / "predictions.tsv", sep="\t", index=False)

    long = pd.concat([pd.DataFrame({"exon": v.exon, "cell": c, "y": v[f"{c}_dpsi"], "base": v[f"{c}_base"],
                                    **{pfx: v[f"{pfx}_{c}"] for pfx in ["M", "S", "C", "N", "K"]}, "R": v.R, "split": v.split})
                      for c in CELLS], ignore_index=True)
    dK = v[["M_K562", "N_K562"]].to_numpy() - v[["M_HepG2", "N_HepG2"]].to_numpy()
    meas_diff = (v.K562_dpsi - v.HepG2_dpsi).to_numpy()

    def stats(ix_long, ix_wide, test_only=False):
        L = long.iloc[ix_long]
        if test_only:
            L = L[L.split == "test"]
        out = {pfx: pear(L[pfx].to_numpy(), L.y.to_numpy()) for pfx in ["M", "S", "C", "R", "N", "K"]}
        out |= {f"mae_{pfx}": float(np.mean(np.abs(L[pfx] - L.y))) for pfx in ["M", "S", "C", "N", "K"]}
        out["between_M"] = pear(dK[ix_wide, 0], meas_diff[ix_wide])
        out["between_N"] = pear(dK[ix_wide, 1], meas_diff[ix_wide])
        return out
    ue = np.unique(v.exon)
    li = {e: np.flatnonzero(long.exon.to_numpy() == e) for e in ue}
    wi = {e: np.flatnonzero(v.exon.to_numpy() == e) for e in ue}
    pt = stats(np.arange(len(long)), np.arange(len(v)))
    pt_test = stats(np.arange(len(long)), np.arange(len(v)), test_only=True)
    big_ex = set(v.loc[(v.K562_base - v.HepG2_base).abs() >= 0.1, "exon"])
    draws, draws_big = [], []
    for _ in range(N_BOOT):
        e = rng.choice(ue, len(ue))
        draws.append(stats(np.concatenate([li[x] for x in e]), np.concatenate([wi[x] for x in e])))
    be = np.array(sorted(big_ex))
    bli = {x: li[x] for x in be}
    for _ in range(N_BOOT):
        e = rng.choice(be, len(be))
        draws_big.append(stats(np.concatenate([bli[x] for x in e]), np.concatenate([wi[x] for x in e])))
    B, BB = pd.DataFrame(draws), pd.DataFrame(draws_big)
    pt_big = stats(np.concatenate([li[x] for x in be]), np.concatenate([wi[x] for x in be]))

    def contrast(Bx, ptx, a_, b_):
        d = (Bx[a_] - Bx[b_]).dropna()
        return {"estimate": ptx[a_] - ptx[b_], "ci": [float(d.quantile(.025)), float(d.quantile(.975))],
                "p_boot": float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean())))}
    rep = {"convention_check_ref_only": conv, "n_variants": len(v), "n_exons": len(ue), "n_exons_baseline_diff_ge_0.1": len(be),
           "calibrated_k": float(k),
           "point_all": pt, "point_heldout_exons_for_K": pt_test, "point_baseline_diff_exons": pt_big,
           "ci_all": {m: [float(B[m].quantile(.025)), float(B[m].quantile(.975))] for m in B.columns},
           "primary": {"M - S (pearson)": contrast(B, pt, "M", "S"), "M - C (pearson)": contrast(B, pt, "M", "C"),
                       "between-context: M - null N": contrast(B, pt, "between_M", "between_N")},
           "secondary": {"M - R (pearson)": contrast(B, pt, "M", "R"), "mae M - S": contrast(B, pt, "mae_M", "mae_S"),
                         "baseline-diff exons: M - S": contrast(BB, pt_big, "M", "S"),
                         "baseline-diff exons: between M - N": contrast(BB, pt_big, "between_M", "between_N")}}
    (OUT / "report.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()

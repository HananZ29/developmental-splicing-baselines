"""E13d (pre-registered): sequence-feature models A (multiplicative, no fit), K (1-parameter logit scale),
B1 (3-parameter logit-additive recalibration of SpliceAI site ratios) and R (raw SpliceAI), on
  1. COMPASS test families (unseen exons AND genes; 5 cell lines),
  2. Vex-seq (independent assay; K562, HepG2),
  3. MFASS (independent assay; HEK293T).
K and B1 are fitted on COMPASS training families only. Every target is replicate-split from its baseline.
Exon-family bootstrap B = 1000, seed 20261105. Primary: B1 - A Pearson per test set (Bonferroni 0.0167).
Output: publication_extension/results/E13_compass/generalization.json, predictions_*.tsv
"""
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import expit, logit
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[2]
PE = ROOT / "publication_extension"
OUT = PE / "results" / "E13_compass"
EPS, SEED, LAM = 0.01, 20261105, 1e-3
NB = int(os.environ.get("E13_NB", "1000"))
TEST_MODE = os.environ.get("E13_TEST_MODE") == "1"   # code-path test only: skips COMPASS, fixed params, writes nothing
CELLS = ["HEK", "HeLa", "K562", "MCF7", "HMC3"]
KEY = ["chrom", "pos", "ref", "alt", "strand"]


def load(name):
    spec = importlib.util.spec_from_file_location(name, PE / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def clipp(x):
    return np.clip(x, 0.01, 0.99)


def pear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 10 and np.std(x[ok]) > 0 else np.nan


def predict(p, la, ld, model, th):
    p = clipp(p)
    if model == "A":
        a = np.exp(la + ld)
        return p * a / (p * a + 1 - p) - p
    if model == "K":
        return expit(logit(p) + th[0] * (la + ld)) - p
    if model == "B1":
        return expit(logit(p) + th[0] + th[1] * la + th[2] * ld) - p
    raise ValueError(model)


def site_ratios(v, tr, W, e11):
    """r_acc, r_don (log) and raw delta at v.acceptor / v.donor after the REF-only modal offset."""
    conv = {}
    for col, kind in [("acceptor", "acc"), ("donor", "don")]:
        off, frac = e11.modal_offset(v, tr, W, col, kind)
        conv[col] = {"modal_offset": off, "frac_at_modal": frac}
        v[col] = v[col] + off
    ar, aa = e11.site(v, tr, W, "acceptor", "acc")
    dr, da = e11.site(v, tr, W, "donor", "don")
    la, ld = np.log((aa + EPS) / (ar + EPS)), np.log((da + EPS) / (dr + EPS))
    return la, ld, (aa - ar) + (da - dr), conv


def metrics(y, pr):
    big = np.abs(y) >= 0.1
    ok = np.isfinite(y) & np.isfinite(pr)
    return {"pearson": pear(pr, y), "mae": float(np.nanmean(np.abs(pr - y))) if pr is not None else np.nan,
            "auprc_abs_dpsi_ge_0.1": float(average_precision_score(big[ok], np.abs(pr[ok]))) if big[ok].any() else np.nan}


def evaluate(df, models, label):
    """df: rows with columns y, exon, and pred_<m>. Returns point, CI, and contrasts."""
    rng = np.random.default_rng(SEED)
    y = df.y.to_numpy()
    pt = {m: metrics(y, df[f"pred_{m}"].to_numpy()) for m in models}
    pt["R"] = {"pearson": pear(df.pred_R.to_numpy(), y)}
    ex = df.exon.to_numpy()
    ue = np.unique(ex)
    ix = {e: np.flatnonzero(ex == e) for e in ue}
    draws = []
    for _ in range(NB):
        s = np.concatenate([ix[e] for e in rng.choice(ue, len(ue))])
        yy = y[s]
        dr = {m: pear(df[f"pred_{m}"].to_numpy()[s], yy) for m in models + ["R"]}
        dr |= {f"mae_{m}": float(np.mean(np.abs(df[f"pred_{m}"].to_numpy()[s] - yy))) for m in models}
        draws.append(dr)
    Bd = pd.DataFrame(draws)

    def con(a, b):
        d = (Bd[a] - Bd[b]).dropna()
        return {"estimate": (pt[a]["pearson"] if not a.startswith("mae") else pt[a[4:]]["mae"]) - (pt[b]["pearson"] if not b.startswith("mae") else pt[b[4:]]["mae"]),
                "ci": [float(d.quantile(.025)), float(d.quantile(.975))], "p_boot": float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean())))}
    return {"label": label, "n_obs": len(df), "n_exons": len(ue),
            "point": pt, "ci_pearson": {m: [float(Bd[m].quantile(.025)), float(Bd[m].quantile(.975))] for m in models + ["R"]},
            "contrasts": {"B1 - A": con("B1", "A"), "K - A": con("K", "A"), "B1 - R": con("B1", "R"), "A - R": con("A", "R"),
                          "mae B1 - A": con("mae_B1", "mae_A")}}


def main():
    e11 = load("E11_vexseq_evaluate")
    rep = {}
    if TEST_MODE:
        thK, thB = np.array([1.0]), np.array([0.0, 1.0, 1.0])
    else:
        thK, thB = compass_part(rep, e11)

    def addp(df):
        for m, th in [("A", None), ("K", thK), ("B1", thB)]:
            df[f"pred_{m}"] = predict(df.p.to_numpy(), df.la.to_numpy(), df.ld.to_numpy(), m, th)
        df["pred_R"] = df.R
        return df
    if not TEST_MODE:
        tst = addp(rep.pop("_compass_test"))
        tst.to_csv(OUT / "predictions_compass_test.tsv", sep="\t", index=False)
        rep["test1_compass_unseen_exons_genes"] = evaluate(tst, ["A", "K", "B1"], "COMPASS test families (gene-disjoint), 5 cells, replicate split")
        rep["test1_per_cell_pearson"] = {c: {m: pear(g[f"pred_{m}"].to_numpy(), g.y.to_numpy()) for m in ["A", "K", "B1", "R"]} for c, g in tst.groupby("cell")}
    vex_mfass(rep, e11, addp)
    if TEST_MODE:
        print(json.dumps({k: rep[k]["point"] for k in ["test2_vexseq", "test3_mfass"]}, indent=1, default=float))
        print(json.dumps({k: rep[k]["contrasts"]["B1 - A"] for k in ["test2_vexseq", "test3_mfass"]}, indent=1, default=float))
        return
    finish(rep)


def compass_part(rep, e11):
    # ---------------- COMPASS subset: features
    sub = pd.read_csv(OUT / "subset_variants.tsv", sep="\t")
    sc = pd.read_csv(OUT / "spliceai" / "scores.tsv", sep="\t")
    sc["strand"] = sc.strand.astype(str)
    sc["track_row"] = np.arange(len(sc))
    sub = sub.merge(sc[KEY + ["track_row"]], on=KEY, how="inner")
    t = np.load(OUT / "spliceai" / "tracks.npz")
    W = int(t["window"])
    tr = {k: t[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    plus = sub.strand == "+"
    sub["acceptor"] = np.where(plus, sub.exon_start_hg19, sub.exon_end_hg19).astype(int)
    sub["donor"] = np.where(plus, sub.exon_end_hg19, sub.exon_start_hg19).astype(int)
    la, ld, R, conv = site_ratios(sub, tr, W, e11)
    ok = np.isfinite(la) & np.isfinite(ld)
    sub, la, ld, R = sub[ok].reset_index(drop=True), la[ok], ld[ok], R[ok]
    sub["la"], sub["ld"], sub["R"] = la, ld, R
    rep["compass_convention"] = conv
    d = pd.read_csv(ROOT / "data" / "external" / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt = d[d.snp == "none"].drop_duplicates("event_id").set_index("event_id")
    W_ = wt.loc[sub.event_id]
    long = []
    for c in CELLS:
        for rin, rtg in [(1, 2), (2, 1)]:
            nv = [sub[f"{c}_rep{r}_included"].to_numpy() + sub[f"{c}_rep{r}_excluded"].to_numpy() for r in (1, 2)]
            nw = [W_[f"{c}_rep{r}_included"].to_numpy() + W_[f"{c}_rep{r}_excluded"].to_numpy() for r in (1, 2)]
            q = (nv[0] >= 20) & (nv[1] >= 20) & (nw[0] >= 20) & (nw[1] >= 20)
            p_in = W_[f"{c}_rep{rin}_psi_raw"].to_numpy()
            y = sub[f"{c}_rep{rtg}_psi_raw"].to_numpy() - W_[f"{c}_rep{rtg}_psi_raw"].to_numpy()
            long.append(pd.DataFrame({"exon": sub.event_id, "gene": sub.gene, "split": sub.split, "cell": c, "p": p_in, "y": y,
                                      "la": sub.la, "ld": sub.ld, "R": sub.R})[q])
    L = pd.concat(long, ignore_index=True).dropna(subset=["p", "y"])
    trn = L[L.split == "train"]

    def fit(model, n):
        def res(th):
            r = predict(trn.p.to_numpy(), trn.la.to_numpy(), trn.ld.to_numpy(), model, th) - trn.y.to_numpy()
            return np.concatenate([r, np.sqrt(LAM * len(r)) * th / np.sqrt(len(r))])
        return least_squares(res, np.ones(n) * (0.5 if model == "K" else 0.0) + (np.array([0, 1, 1]) if model == "B1" else 0)).x
    thK, thB = fit("K", 1), fit("B1", 3)
    rep["fitted_on_compass_train"] = {"K": thK.tolist(), "B1 [b, k_acc, k_don]": thB.tolist(), "n_train_obs": len(trn),
                                      "n_train_exons": int(trn.exon.nunique()), "n_train_genes": int(trn.gene.nunique())}
    rep["_compass_test"] = L[L.split == "test"].copy()
    return thK, thB


def vex_mfass(rep, e11, addp):
    # ---------------- Vex-seq (independent; replicate split, 3 rotations as E11b)
    e11b = load("E11b_vexseq_replicate_split")
    v = pd.read_csv(PE / "results" / "E10_vexseq" / "vexseq_measured.tsv", sep="\t")
    v["strand"] = v.strand.astype(str)
    s2 = pd.read_csv(PE / "results" / "E10_vexseq" / "spliceai" / "scores.tsv", sep="\t")
    s2["strand"] = s2.strand.astype(str)
    s2["track_row"] = np.arange(len(s2))
    v = v.merge(s2[KEY + ["track_row"]], on=KEY, how="inner")
    t2 = np.load(PE / "results" / "E10_vexseq" / "spliceai" / "tracks.npz")
    W2 = int(t2["window"])
    tr2 = {k: t2[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    v["acceptor"] = np.where(v.strand == "+", v["Exon up"], v["Exon down"]).astype(int)
    v["donor"] = np.where(v.strand == "+", v["Exon down"], v["Exon up"]).astype(int)
    la2, ld2, R2, conv2 = site_ratios(v, tr2, W2, e11)
    ok = np.isfinite(la2) & np.isfinite(ld2)
    v, la2, ld2, R2 = v[ok].reset_index(drop=True), la2[ok], ld2[ok], R2[ok]
    rep["vexseq_convention"] = conv2
    ref_id = v.chrom + "_" + v.ctrl_pos.astype(int).astype(str) + "_Control_Original"
    exon = (v.chrom + ":" + v["Exon up"].astype(str)).to_numpy()
    VR = {c: e11b.reps(c) for c in ["K562", "HepG2"]}
    vl = []
    for r in range(3):
        o = [i for i in range(3) if i != r]
        for c in ["K562", "HepG2"]:
            base = VR[c][r].reindex(ref_id).to_numpy()
            y = np.nanmean([VR[c][i].reindex(v.variant).to_numpy() for i in o], 0) - np.nanmean([VR[c][i].reindex(ref_id).to_numpy() for i in o], 0)
            vl.append(pd.DataFrame({"exon": exon, "cell": c, "rot": r, "p": base, "y": y, "la": la2, "ld": ld2, "R": R2}))
    VL = addp(pd.concat(vl, ignore_index=True).dropna(subset=["p", "y"]))
    if not TEST_MODE:
        VL.to_csv(OUT / "predictions_vexseq.tsv", sep="\t", index=False)
    rep["test2_vexseq"] = evaluate(VL, ["A", "K", "B1"], "Vex-seq, K562+HepG2, 3 replicate rotations pooled")

    # ---------------- MFASS (independent; HEK293T; R1 baseline -> R2 target and reverse)
    e5 = load("E5_mfass_evaluate")
    mv = pd.read_csv(PE / "results" / "E5_mfass" / "mfass_variants.tsv", sep="\t")
    raw = pd.read_csv(ROOT / "data" / "external" / "mfass" / "snv_data_clean.txt", sep="\t",
                      usecols=["id", "v2_index_R1", "v2_index_R2", "nat_v2_index_R1", "nat_v2_index_R2"])
    mv = mv.merge(raw, on="id", how="left")
    s3 = pd.read_csv(PE / "results" / "E5_mfass" / "spliceai" / "scores.tsv", sep="\t")
    s3["track_row"] = np.arange(len(s3))
    mv = mv.merge(s3[KEY + ["track_row"]], on=KEY, how="inner")
    t3 = np.load(PE / "results" / "E5_mfass" / "spliceai" / "tracks.npz")
    W3 = int(t3["window"])
    tr3 = {k: t3[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    conv3 = e5.check_convention(mv, tr3, W3)
    for col in ["acceptor", "donor"]:
        mv[col] = mv[col] + conv3[col]["modal_offset"]
    ar, aa = e5.site_values(mv, tr3, W3, "acceptor", "acc")
    dr, da = e5.site_values(mv, tr3, W3, "donor", "don")
    ok = np.isfinite(ar) & np.isfinite(dr)
    mv = mv[ok].reset_index(drop=True)
    ar, aa, dr, da = ar[ok], aa[ok], dr[ok], da[ok]
    rep["mfass_convention"] = conv3
    ml = []
    for rin, rtg in [("R1", "R2"), ("R2", "R1")]:
        ml.append(pd.DataFrame({"exon": mv.exon_id, "p": mv[f"nat_v2_index_{rin}"], "y": mv[f"v2_index_{rtg}"] - mv[f"nat_v2_index_{rtg}"],
                                "la": np.log((aa + EPS) / (ar + EPS)), "ld": np.log((da + EPS) / (dr + EPS)), "R": (aa - ar) + (da - dr)}))
    ML = addp(pd.concat(ml, ignore_index=True).dropna(subset=["p", "y"]))
    if not TEST_MODE:
        ML.to_csv(OUT / "predictions_mfass.tsv", sep="\t", index=False)
    rep["test3_mfass"] = evaluate(ML, ["A", "K", "B1"], "MFASS (E5 sample of 3,000), HEK293T, replicate split both directions")


def finish(rep):
    # ---------------- gate G-B
    wins = [rep[k]["contrasts"]["B1 - A"]["ci"][0] > 0 for k in ["test1_compass_unseen_exons_genes", "test2_vexseq", "test3_mfass"]]
    not_worse_R = [rep[k]["contrasts"]["B1 - R"]["ci"][1] >= 0 for k in ["test1_compass_unseen_exons_genes", "test2_vexseq", "test3_mfass"]]
    rep["gate_G-B"] = {"B1_gt_A_CI_excludes_0": wins, "B1_not_worse_than_R": not_worse_R,
                       "validated": bool(sum(wins) >= 2 and all(not_worse_R))}
    (OUT / "generalization.json").write_text(json.dumps(rep, indent=2, default=float))
    print(json.dumps({k: v for k, v in rep.items()}, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()

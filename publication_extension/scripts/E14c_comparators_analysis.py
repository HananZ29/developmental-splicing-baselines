"""E14c (PREREGISTRATION_E14.md, frozen 2026-11-05): external comparators on identical observations.

Commands:
  check     N_DS interpolation accuracy test on E13c TRAINING variants (ds_native_check.tsv).
            Writes ds_check.json with "pass". Uses no test data.
  analyse   The single pre-registered evaluation (T-fresh primary; T-E13, T-Vex, T-MFASS secondary).
            Run once, after all scoring (steps 2-4) is complete.

Env E14_TEST_MODE=1 (code-path test only): labels are randomly permuted within each set, B = 20,
and nothing is written. The output therefore carries no information about outcomes.

Predictors (§E14b):
- R_X:  Delta(acc) + Delta(don) at the exon sites, for X in SA (SpliceAI), PG (Pangolin), DS (DeltaSplice);
- A_X:  rule A with that tool's site ratios (eps = 0.01);
- N_MM: expit(logit p + MMSplice dlogit_psi) - p;
- N_DS: native DeltaSplice, PCHIP-interpolated on logit(PSI) over the 5-point grid (subsample only);
- K, B1: frozen E13d parameters.
Baselines p are clipped to [0.01, 0.99] for every method.

Interpretation of gate G-X, fixed here before any analysis: for each of N_MM, N_DS and R_PG, the
point estimate of A_SA - comparator must be >= 0.02 AND the Holm-adjusted bootstrap CI must exclude 0.
Holm-adjusted CI: the contrast with the k-th smallest bootstrap p-value (k = 1..m) uses a two-sided
percentile interval at level 1 - 0.05/(m - k + 1).
"""
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.special import expit, logit
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[2]
PE = ROOT / "publication_extension"
OUT = PE / "results" / "E14_comparators"
E13 = PE / "results" / "E13_compass"
SEED, EPS = 20261106, 0.01
TEST = os.environ.get("E14_TEST_MODE") == "1"
NB = 20 if TEST else 1000
CELLS = ["HEK", "HeLa", "K562", "MCF7", "HMC3"]
KEY = ["chrom", "pos", "ref", "alt", "strand"]
GRID = np.array([0.05, 0.2, 0.5, 0.8, 0.95])
THK = np.array([2.164415673688813])
THB = np.array([-0.2521608625197036, 1.650904408290922, 2.5275112344505155])
TOOLS = ["SA", "PG", "DS"]


def load(name):
    spec = importlib.util.spec_from_file_location(name, PE / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


E13d = load("E13d_models_A_B_generalization")
E11 = load("E11_vexseq_evaluate")
clipp, pear, predict = E13d.clipp, E13d.pear, E13d.predict


# ------------------------------------------------------------------ N_DS interpolation
def ds_interp(grid_df):
    """grid_df: rows (KEY, psi, dpsi) -> {variant key: PchipInterpolator on logit psi}."""
    f = {}
    for k, g in grid_df.groupby(KEY):
        g = g.set_index("psi").reindex(GRID)
        if g.dpsi.notna().all():
            f[tuple(k)] = PchipInterpolator(logit(GRID), g.dpsi.to_numpy())
    return f


def ds_predict(f, keys, p):
    lo, hi = logit(GRID[0]), logit(GRID[-1])
    out = np.full(len(p), np.nan)
    for i, (k, pp) in enumerate(zip(keys, p)):
        fn = f.get(k)
        if fn is not None:
            out[i] = float(fn(np.clip(logit(clipp(pp)), lo, hi)))
    return out


def check():
    d = pd.read_csv(OUT / "ds_native_check.tsv", sep="\t")
    d["strand"] = d.strand.astype(str)
    f = ds_interp(d[d.psi.isin(GRID)])
    ex = d[d.psi.isin([0.35, 0.65])]
    keys = [tuple(r) for r in ex[KEY].itertuples(index=False)]
    pred = ds_predict(f, keys, ex.psi.to_numpy())
    err = np.abs(pred - ex.dpsi.to_numpy())
    res = {"n_variants": int(ex[KEY].drop_duplicates().shape[0]), "n_checks": int(np.isfinite(err).sum()),
           "median_abs_err": float(np.nanmedian(err)), "p95_abs_err": float(np.nanquantile(err, 0.95)),
           "max_abs_err": float(np.nanmax(err)), "dpsi_abs_median_at_check": float(np.nanmedian(np.abs(ex.dpsi))),
           "rule": "pass if median <= 0.005 and p95 <= 0.02"}
    res["pass"] = bool(res["median_abs_err"] <= 0.005 and res["p95_abs_err"] <= 0.02)
    (OUT / "ds_check.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


# ------------------------------------------------------------------ tool features
def tracks(dirpath):
    sc = pd.read_csv(dirpath / "scores.tsv", sep="\t")
    sc["strand"] = sc.strand.astype(str)
    t = np.load(dirpath / "tracks.npz")
    tr = {k: t[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    row = {tuple(r): i for i, r in enumerate(sc[KEY].itertuples(index=False))}
    return row, tr, int(t["window"])


def tool_features(v, dirpath, tag, conv_rep, mfass=False):
    """Add la_tag, ld_tag, R_tag to v (needs acceptor, donor, KEY). E11 rule (E5 rule for MFASS)."""
    row, tr, W = tracks(dirpath)
    w = v.copy()
    w["track_row"] = [row.get(tuple(r), -1) for r in w[KEY].itertuples(index=False)]
    has = w.track_row >= 0
    w = w[has].copy()
    if mfass:
        e5 = load("E5_mfass_evaluate")
        conv = e5.check_convention(w, tr, W)
        for col in ["acceptor", "donor"]:
            w[col] = w[col] + conv[col]["modal_offset"]
        ar, aa = e5.site_values(w, tr, W, "acceptor", "acc")
        dr, da = e5.site_values(w, tr, W, "donor", "don")
        la, ld, R = np.log((aa + EPS) / (ar + EPS)), np.log((da + EPS) / (dr + EPS)), (aa - ar) + (da - dr)
    else:
        la, ld, R, conv = E13d.site_ratios(w, tr, W, E11)
    conv_rep[tag] = conv
    out = pd.DataFrame({f"la_{tag}": la, f"ld_{tag}": ld, f"R_{tag}": R}, index=w.index)
    return v.join(out)


def add_mmsplice(v, path, exon_cols):
    m = pd.read_csv(path, sep="\t")
    m["strand"] = m.strand.astype(str)
    m = m.rename(columns={"exon_start": exon_cols[0], "exon_end": exon_cols[1]})
    return v.merge(m[KEY + list(exon_cols) + ["mmsplice_dlogit_psi", "in_range"]], on=KEY + list(exon_cols), how="left")


# ------------------------------------------------------------------ long tables
def compass_long(v):
    d = pd.read_csv(ROOT / "data" / "external" / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt = d[d.snp == "none"].drop_duplicates("event_id").set_index("event_id")
    W_ = wt.loc[v.event_id]
    feat = [c for c in v.columns if c.startswith(("la_", "ld_", "R_"))] + ["mmsplice_dlogit_psi", "in_range"]
    long = []
    for c in CELLS:
        for rin, rtg in [(1, 2), (2, 1)]:
            nv = [v[f"{c}_rep{r}_included"].to_numpy() + v[f"{c}_rep{r}_excluded"].to_numpy() for r in (1, 2)]
            nw = [W_[f"{c}_rep{r}_included"].to_numpy() + W_[f"{c}_rep{r}_excluded"].to_numpy() for r in (1, 2)]
            q = (nv[0] >= 20) & (nv[1] >= 20) & (nw[0] >= 20) & (nw[1] >= 20)
            df = pd.DataFrame({"exon": v.event_id.to_numpy(), "cell": c, "split_dir": f"{rin}>{rtg}",
                               "p": W_[f"{c}_rep{rin}_psi_raw"].to_numpy(),
                               "y": v[f"{c}_rep{rtg}_psi_raw"].to_numpy() - W_[f"{c}_rep{rtg}_psi_raw"].to_numpy()})
            for k in KEY:
                df[k] = v[k].to_numpy()
            for k in feat:
                df[k] = v[k].to_numpy() if k in v else np.nan
            long.append(df[q])
    return pd.concat(long, ignore_index=True).dropna(subset=["p", "y"])


def build_compass(variants_path, prefix, conv_rep):
    v = pd.read_csv(variants_path, sep="\t")
    v["strand"] = v.strand.astype(str)
    plus = v.strand == "+"
    v["acceptor"] = np.where(plus, v.exon_start_hg19, v.exon_end_hg19).astype(int)
    v["donor"] = np.where(plus, v.exon_end_hg19, v.exon_start_hg19).astype(int)
    dirs = {"SA": OUT / f"spliceai_{prefix}", "PG": OUT / f"pangolin_{prefix}", "DS": OUT / f"deltasplice_{prefix}"}
    if prefix == "e13test":
        dirs["SA"] = E13 / "spliceai"
    for t in TOOLS:
        if (dirs[t] / "tracks.npz").exists():
            v = tool_features(v, dirs[t], t, conv_rep)
    v = add_mmsplice(v, OUT / f"mmsplice_{prefix}.tsv", ("exon_start_hg19", "exon_end_hg19")) \
        if (OUT / f"mmsplice_{prefix}.tsv").exists() else v.assign(mmsplice_dlogit_psi=np.nan, in_range=np.nan)
    return compass_long(v)


def build_vex(conv_rep):
    if not ((OUT / "pangolin_vexseq" / "tracks.npz").exists() and (OUT / "mmsplice_vexseq.tsv").exists()):
        return None
    e11b = load("E11b_vexseq_replicate_split")
    v = pd.read_csv(PE / "results" / "E10_vexseq" / "vexseq_measured.tsv", sep="\t")
    v["strand"] = v.strand.astype(str)
    v["acceptor"] = np.where(v.strand == "+", v["Exon up"], v["Exon down"]).astype(int)
    v["donor"] = np.where(v.strand == "+", v["Exon down"], v["Exon up"]).astype(int)
    v["exon_start"], v["exon_end"] = v[["Exon up", "Exon down"]].min(1), v[["Exon up", "Exon down"]].max(1)
    v = tool_features(v, PE / "results" / "E10_vexseq" / "spliceai", "SA", conv_rep)
    v = tool_features(v, OUT / "pangolin_vexseq", "PG", conv_rep)
    if (OUT / "deltasplice_vexseq" / "tracks.npz").exists():
        v = tool_features(v, OUT / "deltasplice_vexseq", "DS", conv_rep)
    v = add_mmsplice(v, OUT / "mmsplice_vexseq.tsv", ("exon_start", "exon_end"))
    ref_id = v.chrom + "_" + v.ctrl_pos.astype(int).astype(str) + "_Control_Original"
    exon = (v.chrom + ":" + v["Exon up"].astype(str)).to_numpy()
    VR = {c: e11b.reps(c) for c in ["K562", "HepG2"]}
    feat = [c for c in v.columns if c.startswith(("la_", "ld_", "R_"))] + ["mmsplice_dlogit_psi", "in_range"]
    vl = []
    for r in range(3):
        o = [i for i in range(3) if i != r]
        for c in ["K562", "HepG2"]:
            base = VR[c][r].reindex(ref_id).to_numpy()
            y = np.nanmean([VR[c][i].reindex(v.variant).to_numpy() for i in o], 0) - \
                np.nanmean([VR[c][i].reindex(ref_id).to_numpy() for i in o], 0)
            df = pd.DataFrame({"exon": exon, "cell": c, "rot": r, "p": base, "y": y})
            for k in KEY + feat:
                df[k] = v[k].to_numpy()
            vl.append(df)
    return pd.concat(vl, ignore_index=True).dropna(subset=["p", "y"])


def build_mfass(conv_rep):
    if not ((OUT / "pangolin_mfass" / "tracks.npz").exists() and (OUT / "mmsplice_mfass.tsv").exists()):
        return None
    mv = pd.read_csv(PE / "results" / "E5_mfass" / "mfass_variants.tsv", sep="\t")
    raw = pd.read_csv(ROOT / "data" / "external" / "mfass" / "snv_data_clean.txt", sep="\t",
                      usecols=["id", "v2_index_R1", "v2_index_R2", "nat_v2_index_R1", "nat_v2_index_R2"])
    mv = mv.merge(raw, on="id", how="left")
    mv["strand"] = mv.strand.astype(str)
    mv["exon_start"], mv["exon_end"] = mv[["acceptor", "donor"]].min(1), mv[["acceptor", "donor"]].max(1)
    mv = tool_features(mv, PE / "results" / "E5_mfass" / "spliceai", "SA", conv_rep, mfass=True)
    mv = tool_features(mv, OUT / "pangolin_mfass", "PG", conv_rep, mfass=True)
    if (OUT / "deltasplice_mfass" / "tracks.npz").exists():
        mv = tool_features(mv, OUT / "deltasplice_mfass", "DS", conv_rep, mfass=True)
    mv = add_mmsplice(mv, OUT / "mmsplice_mfass.tsv", ("exon_start", "exon_end"))
    feat = [c for c in mv.columns if c.startswith(("la_", "ld_", "R_"))] + ["mmsplice_dlogit_psi", "in_range"]
    ml = []
    for rin, rtg in [("R1", "R2"), ("R2", "R1")]:
        df = pd.DataFrame({"exon": mv.exon_id, "p": mv[f"nat_v2_index_{rin}"],
                           "y": mv[f"v2_index_{rtg}"] - mv[f"nat_v2_index_{rtg}"]})
        for k in KEY + feat:
            df[k] = mv[k].to_numpy()
        ml.append(df)
    return pd.concat(ml, ignore_index=True).dropna(subset=["p", "y"])


# ------------------------------------------------------------------ predictions
def add_predictions(L, nds=None):
    p = L.p.to_numpy()
    pc = clipp(p)
    for t in TOOLS:
        if f"la_{t}" in L:
            L[f"A_{t}"] = predict(p, L[f"la_{t}"].to_numpy(), L[f"ld_{t}"].to_numpy(), "A", None)
            L[f"Rr_{t}"] = L[f"R_{t}"]
    L["K"] = predict(p, L.la_SA.to_numpy(), L.ld_SA.to_numpy(), "K", THK)
    L["B1"] = predict(p, L.la_SA.to_numpy(), L.ld_SA.to_numpy(), "B1", THB)
    L["N_MM"] = expit(logit(pc) + L.mmsplice_dlogit_psi.to_numpy()) - pc
    if nds is not None:
        keys = [tuple(r) for r in L[KEY].itertuples(index=False)]
        L["N_DS"] = ds_predict(nds, keys, p)
    return L


# ------------------------------------------------------------------ evaluation
def evaluate(df, preds, contrasts, label):
    rng = np.random.default_rng(SEED)
    y = df.y.to_numpy()
    if TEST:
        y = rng.permutation(y)
    P = {m: df[m].to_numpy() for m in preds}
    raw = {m for m in preds if m.startswith("Rr_")}
    big = np.abs(y) >= 0.1
    pt = {}
    for m in preds:
        ok = np.isfinite(P[m]) & np.isfinite(y)
        pt[m] = {"pearson": pear(P[m], y), "n": int(ok.sum()),
                 "auprc_abs_dpsi_ge_0.1": float(average_precision_score(big[ok], np.abs(P[m][ok]))) if big[ok].any() else None}
        if m not in raw:
            pt[m]["mae"] = float(np.mean(np.abs(P[m][ok] - y[ok])))
    ex = df.exon.to_numpy()
    ue = np.unique(ex)
    ix = {e: np.flatnonzero(ex == e) for e in ue}
    draws = []
    for _ in range(NB):
        s = np.concatenate([ix[e] for e in rng.choice(ue, len(ue))])
        draws.append({m: pear(P[m][s], y[s]) for m in preds})
    Bd = pd.DataFrame(draws)
    cons = {}
    for a, b in contrasts:
        d = (Bd[a] - Bd[b]).dropna()
        cons[f"{a} - {b}"] = {"estimate": pt[a]["pearson"] - pt[b]["pearson"], "draws": d.to_numpy(),
                              "ci95": [float(d.quantile(.025)), float(d.quantile(.975))],
                              "p_boot": float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean())))}
    return {"label": label, "n_obs": len(df), "n_families": len(ue), "point": pt,
            "ci_pearson": {m: [float(Bd[m].quantile(.025)), float(Bd[m].quantile(.975))] for m in preds},
            "contrasts": cons}


def holm(cons, names, alpha=0.05):
    order = sorted(names, key=lambda n: cons[n]["p_boot"])
    m, prev = len(order), 0.0
    for k, n in enumerate(order):
        lev = alpha / (m - k)
        d = cons[n]["draws"]
        cons[n]["holm_level"] = lev
        cons[n]["ci_holm"] = [float(np.quantile(d, lev / 2)), float(np.quantile(d, 1 - lev / 2))]
        prev = max(prev, min(1.0, cons[n]["p_boot"] * (m - k)))
        cons[n]["p_holm"] = prev
    return cons


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "draws"}
    return o


# ------------------------------------------------------------------ main analysis
def analyse():
    rep = {"test_mode": TEST, "conventions": {}}
    ds_ok = (OUT / "ds_check.json").exists() and json.loads((OUT / "ds_check.json").read_text())["pass"]
    nds = None
    if ds_ok and (OUT / "ds_native_fresh.tsv").exists():
        g = pd.read_csv(OUT / "ds_native_fresh.tsv", sep="\t")
        g["strand"] = g.strand.astype(str)
        nds = ds_interp(g)
    rep["N_DS_available"] = nds is not None

    # ---- T-fresh (primary)
    conv = {}
    L = add_predictions(build_compass(OUT / "fresh_variants.tsv", "fresh", conv), nds)
    rep["conventions"]["T-fresh"] = conv
    need = [f"{k}_{t}" for t in TOOLS for k in ("la", "ld")] + ["mmsplice_dlogit_psi"]
    elig = np.all(np.isfinite(L[need].to_numpy(dtype=float)), axis=1)
    rep["T-fresh_rows"] = {"all": len(L), "eligible_all_tools": int(elig.sum()),
                           "dropped": {c: int((~np.isfinite(L[c].to_numpy(dtype=float))).sum()) for c in need}}
    F = L[elig].copy()
    preds = ["A_SA", "A_PG", "A_DS", "Rr_SA", "Rr_PG", "Rr_DS", "N_MM", "K", "B1"]
    p1 = [("A_SA", "N_MM"), ("A_SA", "Rr_PG")]
    p2 = [("A_SA", "Rr_SA"), ("A_PG", "Rr_PG"), ("A_DS", "Rr_DS")]
    res = evaluate(F, preds, p1 + p2 + [("A_SA", "A_PG"), ("A_SA", "A_DS"), ("A_SA", "B1")], "T-fresh (all rows)")
    res["mmsplice_in_range_frac"] = float(F.in_range.astype(float).mean())
    if nds is not None:
        S = F[np.isfinite(F.N_DS)].copy()
        sub = evaluate(S, preds + ["N_DS"], [("A_SA", "N_DS")] + p1 + p2, "T-fresh N_DS subsample")
        res["N_DS_subsample"] = sub
        res["contrasts"]["A_SA - N_DS"] = sub["contrasts"]["A_SA - N_DS"]
    P1 = [f"{a} - {b}" for a, b in p1] + (["A_SA - N_DS"] if nds is not None else [])
    holm(res["contrasts"], P1)
    holm(res["contrasts"], [f"{a} - {b}" for a, b in p2])
    gx = {n: bool(res["contrasts"][n]["estimate"] >= 0.02 and res["contrasts"][n]["ci_holm"][0] > 0) for n in P1}
    better = {n: bool(res["contrasts"][n]["ci_holm"][1] < 0) for n in P1}
    rep["gate_G-X"] = {"per_comparator_pass": gx, "comparator_significantly_better": better,
                       "N_DS_included": nds is not None,
                       "pass": bool(nds is not None and all(gx.values())),
                       "note": "pass requires all three comparators, including N_DS"}
    res["per_cell_pearson"] = {c: {m: pear(g[m].to_numpy(), g.y.to_numpy()) for m in preds} for c, g in F.groupby("cell")}
    io = F.in_range.astype(bool)
    res["mmsplice_in_range_subset"] = {m: pear(F[m].to_numpy()[io], F.y.to_numpy()[io]) for m in preds}
    rep["T-fresh"] = res
    if not TEST:
        F.to_csv(OUT / "predictions_fresh.tsv", sep="\t", index=False)

    # ---- secondary: T-E13 (previously seen)
    conv = {}
    if (OUT / "mmsplice_e13test.tsv").exists():
        e = pd.read_csv(E13 / "subset_variants.tsv", sep="\t")
        tmp = OUT / "_e13test_variants.tsv"
        e[e.split == "test"].to_csv(tmp, sep="\t", index=False)
        L2 = add_predictions(build_compass(tmp, "e13test", conv))
        tmp.unlink()
        have = [t for t in TOOLS if f"la_{t}" in L2]
        need2 = [f"{k}_{t}" for t in have for k in ("la", "ld")] + ["mmsplice_dlogit_psi"]
        L2 = L2[np.all(np.isfinite(L2[need2].to_numpy(dtype=float)), axis=1)]
        pr2 = [f"A_{t}" for t in have] + [f"Rr_{t}" for t in have] + ["N_MM", "K", "B1"]
        c2 = [(a, b) for a, b in p1 + p2 if a in pr2 and b in pr2]
        rep["T-E13 (previously seen)"] = evaluate(L2, pr2, c2, "E13c test families, previously seen")
        rep["conventions"]["T-E13"] = conv
    # ---- secondary: T-Vex and T-MFASS (previously seen), replicate split exactly as E13d
    for name, builder in [("T-Vex (previously seen; MMSplice fitted on Vex-seq: not interpreted)", build_vex),
                          ("T-MFASS (previously seen)", build_mfass)]:
        conv = {}
        L3 = builder(conv)
        if L3 is None:
            rep[name] = "not run: Pangolin/MMSplice scores missing"
            continue
        L3 = add_predictions(L3)
        have = [t for t in TOOLS if f"la_{t}" in L3]
        need3 = [f"{k}_{t}" for t in have for k in ("la", "ld")] + ["mmsplice_dlogit_psi"]
        L3 = L3[np.all(np.isfinite(L3[need3].to_numpy(dtype=float)), axis=1)]
        pr3 = [f"A_{t}" for t in have] + [f"Rr_{t}" for t in have] + ["N_MM", "K", "B1"]
        c3 = [(a, b) for a, b in p1 + p2 if a in pr3 and b in pr3]
        rep[name] = evaluate(L3, pr3, c3, name)
        rep["conventions"][name] = conv
    out = strip(rep)
    if TEST:
        print(json.dumps({"rows": rep["T-fresh_rows"], "n_fam": res["n_families"],
                          "point_r_PERMUTED": {m: round(v["pearson"], 3) for m, v in res["point"].items()},
                          "contrasts": list(res["contrasts"]), "gate_keys": list(gx)}, indent=1, default=float))
        return
    (OUT / "comparators.json").write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=1, default=float)[:8000])


if __name__ == "__main__":
    {"check": check, "analyse": analyse}[sys.argv[1]]()

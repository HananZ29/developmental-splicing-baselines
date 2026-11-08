"""Baseline ablation: propagate SpliceAI predictions through fetal / adult / generic isoform baselines
and evaluate against developing-brain QTL ground truth.

Conditions (only the baseline pi differs; splice scorer, window, eps, alpha, isoform universe and
evaluation set are identical):
  fetal       Leung et al. fetal cortex Iso-Seq, pooled FL proportions (FSM, ISM excluded)
  adult       Leung et al. adult cortex Iso-Seq, same definition
  generic     GTEx v8 mean over 54 tissues of per-tissue mean transcript proportions
  gtex_cortex GTEx v8 Brain - Cortex (independent short-read adult reference; sensitivity)
  uniform     equal weights over annotated isoforms (no-baseline control)

Evaluations (pre-specified; see docs/methods.md for what each one does and does not test):
  sQTL arm  (primary)  : predicted change in lead-intron usage vs observed sQTL z (ALT allele)
  isoQTL arm (secondary): predicted change in lead-isoform proportion vs observed isoQTL z, where the
                          observed phenotype is isoform ABUNDANCE (VST), not proportion
  Metrics: direction accuracy, Spearman rho (signed), AUPRC/AUROC of lead variants vs matched controls.
  Uncertainty: cluster bootstrap over 1-Mb genomic blocks (B=2000), paired across conditions.

Usage: python experiments/08_propagate_evaluate.py
Outputs: results/experiment/*
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from propagation import (ALPHA, EPS, intron_membership, intron_usage, propagate,  # noqa: E402
                         smooth_baseline, transcript_factors)
from splice_graph import load_lift37  # noqa: E402

OUT = Path(os.environ.get("EXPERIMENT_OUT", ROOT / "results" / "experiment_v2"))
SPLICEAI = Path(os.environ.get("SPLICEAI_DIR", ROOT / "results" / "spliceai"))
CONDITIONS = ["fetal", "adult", "generic", "gtex_cortex", "uniform"]
PRIMARY = ["fetal", "adult", "generic"]
PAIRS = [("fetal", "adult"), ("fetal", "generic"), ("adult", "generic"), ("fetal", "uniform")]
INFORMATIVE = 0.05   # baseline-independent subset: max |SpliceAI delta| at an annotated site >= 0.05
# v2 (primary): LeafCutter-style connected-component intron clusters; v1: CLUSTER_MODE=direct
CLUSTER_MODE = os.environ.get("CLUSTER_MODE", "component")
REAL_WEIGHT = 0.01   # post-hoc audit subset: target weight >= 0.01 under both fetal and adult
N_BOOT, SEED = 2000, 11
KEY = ["chrom", "pos", "ref", "alt", "strand"]


def load_baselines():
    raw = {}
    for stage in ["fetal", "adult"]:
        t = pd.read_csv(ROOT / "results" / "proportions" / f"{stage}_props_ism-exclude.tsv.gz", sep="\t")
        t = t[t["feature"] != "NOVEL"]
        raw[stage] = dict(zip(t["feature"].str.split(".").str[0], t["prop_pooled"]))
    g = pd.read_csv(ROOT / "results" / "gtex" / "baselines.tsv.gz", sep="\t")
    raw["generic"] = dict(zip(g["enst"], g["generic_all"].fillna(0)))
    raw["gtex_cortex"] = dict(zip(g["enst"], g["gtex_cortex"].fillna(0)))
    raw["uniform"] = {}
    return raw


def predict_all(v, by_gene, tracks, window, raw, eps=EPS, alpha=ALPHA, max_dist=None):
    rows = []
    for r in v.itertuples(index=False):
        ts = sorted(by_gene[r.ensg], key=lambda t: t.enst)
        names = [t.enst for t in ts]
        tr = {k: tracks[k][r.track_row] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
        f = transcript_factors(ts, r.pos, tr, window, eps, max_dist)
        site_delta, annotated = 0.0, set()
        for t in ts:
            for pos, kind in t.sites():
                o = pos - r.pos
                annotated.add((o, kind))
                if abs(o) <= window:
                    key = "acc" if kind == "acceptor" else "don"
                    site_delta = max(site_delta, abs(tr[f"{key}_alt"][window + o] - tr[f"{key}_ref"][window + o]))
        # strongest SpliceAI event within +-50 nt: is it at an annotated site the model can use?
        ev = max([("acceptor", r.DS_AG, r.DP_AG), ("acceptor", r.DS_AL, r.DP_AL),
                  ("donor", r.DS_DG, r.DP_DG), ("donor", r.DS_DL, r.DP_DL)], key=lambda e: e[1])
        out = dict(arm=r.arm, role=r.role, ensg=r.ensg, variant_id=r.variant_id, chrom=r.chrom, pos=r.pos,
                   z_obs=r.z_obs, site_maxdelta=site_delta, DS_max=r.DS_max, n_tx=len(ts),
                   n_tx_with_factor_ne_1=int((np.abs(np.log(f)) > 1e-6).sum()),
                   top_event_at_annotated_site=(ev[2], ev[0]) in annotated, top_event_ds=ev[1],
                   strict_control=bool(getattr(r, "strict_control", True)))
        if r.arm == "sQTL":
            introns, index, m, comp = intron_membership(ts, CLUSTER_MODE)
            ok = np.array([len(c) >= 2 for c in comp])
            tgt = None
            if r.role == "positive":
                s, e = map(int, str(r.target).split("-"))
                tgt = index[(s, e)]
        else:
            tgt = names.index(r.target) if r.role == "positive" else None
        for c in CONDITIONS:
            pi = smooth_baseline(raw[c], names, alpha)
            new = propagate(pi, f)
            if r.arm == "sQTL":
                du = intron_usage(new, m, comp) - intron_usage(pi, m, comp)
                out[f"score_{c}"] = float(np.nanmax(np.abs(du[ok]))) if ok.any() else 0.0
                out[f"pred_{c}"] = float(du[tgt]) if tgt is not None else np.nan
                if tgt is not None:
                    out[f"base_{c}"] = float(intron_usage(pi, m, comp)[tgt])
            else:
                d = new - pi
                out[f"score_{c}"] = float(np.abs(d).max())
                out[f"pred_{c}"] = float(d[tgt]) if tgt is not None else np.nan
                if tgt is not None:
                    out[f"base_{c}"] = float(pi[tgt])
        rows.append(out)
    return pd.DataFrame(rows)


def direction_acc(pred, z):
    s = np.sign(pred)
    return float(np.mean(np.where(s == 0, 0.5, s == np.sign(z))))


def metrics(df, conds):
    pos = df[df["role"] == "positive"]
    inf = pos[pos["site_maxdelta"] >= INFORMATIVE]
    res = {}
    for c in conds:
        m = {"dir_acc_all": direction_acc(pos[f"pred_{c}"], pos["z_obs"]),
             "spearman_all": spearmanr(pos[f"pred_{c}"], pos["z_obs"]).statistic,
             # post-hoc (added in the audit, not pre-specified):
             "posthoc_spearman_abs": spearmanr(pos[f"pred_{c}"].abs(), pos["z_obs"].abs()).statistic}
        real = pos[(pos["base_fetal"] >= REAL_WEIGHT) & (pos["base_adult"] >= REAL_WEIGHT)]
        if len(real) >= 5:
            m["posthoc_dir_acc_realweight"] = direction_acc(real[f"pred_{c}"], real["z_obs"])
            m["posthoc_spearman_realweight"] = spearmanr(real[f"pred_{c}"], real["z_obs"]).statistic
        if len(inf) >= 5:
            m["dir_acc_informative"] = direction_acc(inf[f"pred_{c}"], inf["z_obs"])
            m["spearman_informative"] = spearmanr(inf[f"pred_{c}"], inf["z_obs"]).statistic
        y = (df["role"] == "positive").astype(int)
        m["auprc"] = average_precision_score(y, df[f"score_{c}"])
        m["auroc"] = roc_auc_score(y, df[f"score_{c}"])
        strict = df[(df["role"] == "positive") | df["strict_control"]]
        ys = (strict["role"] == "positive").astype(int)
        m["auprc_strict_controls"] = average_precision_score(ys, strict[f"score_{c}"])
        m["auroc_strict_controls"] = roc_auc_score(ys, strict[f"score_{c}"])
        res[c] = m
    return res


def reference_scores(df):
    y = (df["role"] == "positive").astype(int)
    return {k: {"auprc": average_precision_score(y, df[k]), "auroc": roc_auc_score(y, df[k])}
            for k in ["DS_max", "site_maxdelta"]}


def cluster_bootstrap(df, conds, rng):
    block = df["chrom"] + ":" + (df["pos"] // 1_000_000).astype(str)
    groups = {b: idx for b, idx in df.groupby(block).groups.items()}
    keys = list(groups)
    draws = []
    for _ in range(N_BOOT):
        pick = rng.choice(len(keys), len(keys), replace=True)
        sub = df.loc[np.concatenate([groups[keys[i]] for i in pick])]
        if sub["role"].nunique() < 2 or (sub["role"] == "positive").sum() < 5:
            continue
        m = metrics(sub, conds)
        draws.append({f"{c}|{k}": v for c, d in m.items() for k, v in d.items()})
    return pd.DataFrame(draws), len(keys)


def summarize_boot(point, boot, conds):
    rows = []
    metric_names = sorted({k for d in point.values() for k in d})
    for k in metric_names:
        for c in conds:
            if f"{c}|{k}" not in boot:
                continue
            b = boot[f"{c}|{k}"].dropna()
            rows.append(dict(metric=k, condition=c, estimate=point[c][k], ci_low=b.quantile(.025), ci_high=b.quantile(.975)))
        for a, c in PAIRS:
            if f"{a}|{k}" not in boot:
                continue
            d = (boot[f"{a}|{k}"] - boot[f"{c}|{k}"]).dropna()
            rows.append(dict(metric=k, condition=f"{a} - {c}", estimate=point[a][k] - point[c][k],
                             ci_low=d.quantile(.025), ci_high=d.quantile(.975),
                             p_boot=min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    v = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    sc = pd.read_csv(SPLICEAI / "scores.tsv", sep="\t")
    sc["track_row"] = np.arange(len(sc))
    n_eligible = len(v)
    ds_cols = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL", "DS_max"]
    v = v.merge(sc[KEY + ["track_row"] + ds_cols], on=KEY, how="inner")
    # strict controls: gene not significant in the other QTL arm (label audit, 09_validate_labels.py)
    flag = pd.read_csv(ROOT / "results" / "variants" / "controls_flagged.tsv", sep="\t")
    flag = flag.set_index(["arm", "variant_id", "ensg"])["gene_sig_in_other_arm"]
    flag = flag[~flag.index.duplicated()]
    key = list(zip(v["arm"], v["variant_id"], v["ensg"]))
    v["strict_control"] = (v["role"] == "control") & ~pd.Series(flag.reindex(key).fillna(True).to_numpy(dtype=bool))
    tr = np.load(SPLICEAI / "tracks.npz")
    tracks = {k: tr[k] for k in ["acc_ref", "acc_alt", "don_ref", "don_alt"]}
    window = int(tr["window"])
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", set(v["ensg"]))
    raw = load_baselines()

    preds = predict_all(v, by_gene, tracks, window, raw)
    preds.to_csv(OUT / "predictions.tsv", sep="\t", index=False)

    report, tables = {}, []
    for arm in ["sQTL", "isoQTL"]:
        df = preds[preds["arm"] == arm].reset_index(drop=True)
        pos = df[df["role"] == "positive"]
        point = metrics(df, CONDITIONS)
        boot, n_blocks = cluster_bootstrap(df, CONDITIONS, rng)
        tab = summarize_boot(point, boot, CONDITIONS).assign(arm=arm)
        tables.append(tab)
        inf = pos[pos["site_maxdelta"] >= INFORMATIVE]
        binom = {c: binomtest(int(((np.sign(inf[f"pred_{c}"]) == np.sign(inf["z_obs"])) & (inf[f"pred_{c}"] != 0)).sum()),
                              int((inf[f"pred_{c}"] != 0).sum()), 0.5).pvalue for c in CONDITIONS} if len(inf) else {}
        flips = {f"{a} vs {c}": int((np.sign(pos[f"pred_{a}"]) != np.sign(pos[f"pred_{c}"])).sum()) for a, c in PAIRS}
        report[arm] = {
            "n_positive": int(len(pos)), "n_control": int((df["role"] == "control").sum()),
            "n_strict_control": int(df["strict_control"].sum()),
            "positives_with_any_tx_factor_ne_1": int((pos["n_tx_with_factor_ne_1"] > 0).sum()),
            "positives_abs_pred_ge_0.01": {c: int((pos[f"pred_{c}"].abs() >= 0.01).sum()) for c in CONDITIONS},
            "positives_top_event_at_annotated_site": int(pos["top_event_at_annotated_site"].sum()),
            "positives_top_event_ds_ge_0.2_unannotated": int(((pos["top_event_ds"] >= 0.2)
                                                              & ~pos["top_event_at_annotated_site"]).sum()),
            "n_genes_positive": int(pos["ensg"].nunique()), "n_genes_total": int(df["ensg"].nunique()),
            "n_isoforms_in_universe_total": int(df.drop_duplicates("ensg")["n_tx"].sum()),
            "n_positive_informative": int(len(inf)), "n_blocks": n_blocks, "n_boot_used": int(len(boot)),
            "baseline_sign_flips_among_positives": flips,
            "median_abs_pred_positive": {c: float(pos[f"pred_{c}"].abs().median()) for c in CONDITIONS},
            "binomial_p_direction_informative": binom,
            "reference_scores": reference_scores(df),
            "point": point,
        }
    table = pd.concat(tables)
    table.to_csv(OUT / "metrics_with_ci.tsv", sep="\t", index=False, float_format="%.4g")

    # failed examples: informative positives where every primary baseline predicts the wrong sign
    inf = preds[(preds["role"] == "positive") & (preds["site_maxdelta"] >= INFORMATIVE)]
    wrong = inf[np.all([np.sign(inf[f"pred_{c}"]) != np.sign(inf["z_obs"]) for c in PRIMARY], axis=0)]
    wrong.sort_values("site_maxdelta", ascending=False).to_csv(OUT / "failed_examples.tsv", sep="\t", index=False)

    # sensitivity grid, declared before results were seen: eps x alpha x site distance
    sens = []
    for eps in [0.001, 0.01, 0.1]:
        for alpha in [1e-4, 1e-3, 1e-2]:
            for max_dist in [50, window]:
                p = predict_all(v, by_gene, tracks, window, raw, eps=eps, alpha=alpha, max_dist=max_dist)
                for arm in ["sQTL", "isoQTL"]:
                    for c, m in metrics(p[p["arm"] == arm].reset_index(drop=True), CONDITIONS).items():
                        sens.append(dict(eps=eps, alpha=alpha, max_dist=max_dist, arm=arm, condition=c, **m))
    pd.DataFrame(sens).to_csv(OUT / "sensitivity_grid.tsv", sep="\t", index=False, float_format="%.4g")

    report["cluster_mode"] = CLUSTER_MODE
    report["counts"] = {"eligible_rows": n_eligible, "scored_rows": int(len(v)),
                        "unique_variant_strands_scored": int(len(sc)), "failed_examples_all_wrong": int(len(wrong))}
    (OUT / "report.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=float))
    print(table.to_string())


if __name__ == "__main__":
    main()

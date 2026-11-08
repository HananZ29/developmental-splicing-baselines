"""Figures for the publication extension (static PNG; same palette as results/figures).

FigE1 root cause: structural sign invariance and the metric ceiling (E2)
FigE2 stage-swap: fetal - adult within prenatal (Wen) and adult (GTEx) sQTL, plus interaction (E3)
FigE3 fetal baselines: coverage and performance of Leung vs Patowary fetal baselines (E4)
FigE4 MFASS: predicted vs measured dPSI and baseline ablation (E5), if results exist
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import json  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
R, F = ROOT / "publication_extension" / "results", ROOT / "publication_extension" / "figures"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = {"fetal": "#2a78d6", "adult": "#eb6834", "generic": "#1baf7a", "gtex_cortex": "#eda100", "uniform": "#8a8986",
     "patowary_fetal": "#6250d6", "patowary_gz": "#e87ba4", "patowary_cp": "#008300"}


def style(ax, grid="x"):
    ax.set_facecolor(SURFACE)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def fig_e1():
    s = json.loads((R / "E2_root_cause" / "summary.json").read_text())
    pv = pd.read_csv(R / "E2_root_cause" / "per_variant.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), facecolor=SURFACE)
    ax = axes[0]
    style(ax)
    arms = ["sQTL", "isoQTL"]
    inv = [s[a]["sign_invariant_by_structure"] / s[a]["n_positive"] for a in arms]
    ax.barh(arms, inv, color=C["fetal"], height=0.5)
    ax.barh(arms, [1 - x for x in inv], left=inv, color=GRID, height=0.5)
    for i, x in enumerate(inv):
        ax.text(x / 2, i, f"{x:.0%} sign fixed\nfor any baseline", ha="center", va="center", fontsize=8, color="white")
    ax.set_xlim(0, 1)
    ax.set_title("Predicted sign fixed by transcript structure", fontsize=9, color=INK, loc="left")
    ax = axes[1]
    style(ax)
    ceil = [s[a]["max_possible_abs_diff_dir_acc_fetal_vs_adult"] for a in arms]
    ax.barh(arms, ceil, color=C["adult"], height=0.5)
    for i, x in enumerate(ceil):
        ax.text(x + 0.001, i, f"±{x:.3f} ({s[arms[i]]['observed_sign_flips_fetal_vs_adult']} sign flips)", va="center", fontsize=8, color=INK2)
    ax.set_xlim(0, 0.06)
    ax.set_title("Largest possible |fetal − adult| direction accuracy", fontsize=9, color=INK, loc="left")
    ax = axes[2]
    style(ax, "both")
    a = pv[pv.arm == "sQTL"]
    ax.scatter(np.abs(a.pred_fetal) + 1e-7, a.pred_absmax + 1e-7, s=12, color=C["fetal"], ec=SURFACE, lw=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    lim = [1e-7, 1]
    ax.plot(lim, lim, color=INK2, lw=1, ls=(0, (3, 3)))
    ax.set_xlabel("|Δu| under fetal baseline", fontsize=8, color=INK2)
    ax.set_ylabel("max |Δu| over 1,000 random baselines", fontsize=8, color=INK2)
    ax.set_title("…but magnitude is highly baseline-dependent (sQTL)", fontsize=9, color=INK, loc="left")
    fig.suptitle("Root cause: the propagation model fixes direction and rank; the baseline changes only magnitude",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(F / "figE1_root_cause.png", dpi=160, facecolor=SURFACE)


def fig_e2():
    t = pd.read_csv(R / "E34" / "metrics_with_ci.tsv", sep="\t")
    it = pd.read_csv(R / "E34" / "stage_swap_interaction.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), facecolor=SURFACE, sharey=True)
    for ax, m, lab in zip(axes, ["dir_acc_all", "spearman_all"], ["Direction accuracy", "Spearman ρ"]):
        style(ax)
        rows = []
        for arm, nm in [("Wen_sQTL", "Prenatal sQTL (Wen; n=147)"), ("GTEx_sQTL", "Adult sQTL (GTEx; n=82)")]:
            r = t[(t.arm == arm) & (t.metric == m) & (t.condition == "fetal - adult")].iloc[0]
            rows.append((nm + ": fetal − adult", r.estimate, r.ci_low, r.ci_high, C["fetal"]))
        r = it[(it.fetal_baseline == "fetal") & (it.metric == m)].iloc[0]
        rows.append(("Interaction (prenatal − adult)", r.interaction, r.ci_low, r.ci_high, INK))
        for i, (nm, e, lo, hi, col) in enumerate(rows):
            ax.plot([lo, hi], [i, i], color=col, lw=2, solid_capstyle="round")
            ax.plot(e, i, "o", ms=7, color=col, mec=SURFACE, mew=1.5)
            ax.text(e, i - 0.25, f"{e:+.3f}", fontsize=7, color=INK2, ha="center", va="bottom")
        ax.axvline(0, color=INK2, lw=1)
        ax.set_yticks(range(len(rows)), [x[0] for x in rows], fontsize=8, color=INK)
        ax.set_ylim(len(rows) - 0.5, -0.8)
        ax.set_title(lab, fontsize=9, color=INK, loc="left")
        ax.set_xlabel("Difference (95% block-bootstrap CI)", fontsize=8, color=INK2)
    fig.suptitle("Stage-swap: a stage-matched baseline should help more on prenatal than adult associations (interaction > 0)",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE2_stage_swap.png", dpi=160, facecolor=SURFACE)


def fig_e3():
    rep = json.loads((R / "E34" / "report.json").read_text())
    t = pd.read_csv(R / "E34" / "metrics_with_ci.tsv", sep="\t")
    conds = ["fetal", "patowary_fetal", "patowary_gz", "patowary_cp", "adult", "generic", "uniform"]
    lab = {"fetal": "Leung fetal", "patowary_fetal": "Patowary fetal", "patowary_gz": "Patowary GZ", "patowary_cp": "Patowary CP",
           "adult": "Leung adult", "generic": "GTEx generic*", "uniform": "Uniform"}
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9), facecolor=SURFACE)
    ax = axes[0]
    style(ax)
    cov = rep["coverage_agreement"]["nonzero_weight_fraction"]
    cs = [c for c in conds if c in cov]
    ax.barh([lab[c] for c in cs], [cov[c] for c in cs], color=[C[c] for c in cs], height=0.6)
    for i, c in enumerate(cs):
        ax.text(cov[c] + 0.01, i, f"{cov[c]:.0%}", va="center", fontsize=7, color=INK2)
    ax.set_xlim(0, 1.1)
    ax.invert_yaxis()
    ax.set_title("Transcripts with nonzero weight\n(863 evaluation genes)", fontsize=9, color=INK, loc="left")
    for ax, arm in zip(axes[1:], ["Wen_sQTL", "Wen_isoQTL", "GTEx_sQTL"]):
        style(ax)
        for i, c in enumerate(conds):
            r = t[(t.arm == arm) & (t.metric == "spearman_all") & (t.condition == c)].iloc[0]
            ax.plot([r.ci_low, r.ci_high], [i, i], color=C[c], lw=2, solid_capstyle="round")
            ax.plot(r.estimate, i, "o", ms=6, color=C[c], mec=SURFACE, mew=1.5)
        ax.set_yticks(range(len(conds)), [lab[c] for c in conds], fontsize=8, color=INK)
        ax.set_ylim(len(conds) - 0.5, -0.5)
        ax.set_title(f"{arm}: Spearman ρ(pred, z)", fontsize=9, color=INK, loc="left")
    fig.suptitle("A 3× better-covered, independent fetal baseline (Patowary) does not change performance",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE3_fetal_baselines.png", dpi=160, facecolor=SURFACE)


def fig_e4():
    p = R / "E5_mfass" / "predictions.tsv"
    if not p.exists():
        return
    d = pd.read_csv(p, sep="\t")
    test = d[d.split == "test"]
    met = pd.read_csv(R / "E5_mfass" / "metrics.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), facecolor=SURFACE)
    for ax, col, nm in [(axes[0], "P0_raw", "SpliceAI site change (P0)"), (axes[1], "P1_prop_measured", "Propagated, measured baseline (P1)")]:
        style(ax, "both")
        ax.scatter(test[col], test.v2_dpsi, s=6, alpha=0.5, color=C["fetal"], lw=0)
        r = met[(met.predictor == col) & (met.metric == "pearson")].iloc[0]
        ax.set_title(f"{nm}\nPearson r = {r.estimate:.3f} [{r.ci_low:.3f}, {r.ci_high:.3f}]", fontsize=9, color=INK, loc="left")
        ax.set_xlabel("Predicted", fontsize=8, color=INK2)
        ax.set_ylabel("Measured ΔPSI (MFASS)", fontsize=8, color=INK2)
    ax = axes[2]
    style(ax)
    order = ["P1_prop_measured", "P5_calibrated", "P2_prop_constant", "P3_prop_shuffled", "P4_prop_half", "P6_null_permuted_factor", "P0_raw"]
    nm = {"P1_prop_measured": "Measured baseline", "P5_calibrated": "Measured + calibrated k", "P2_prop_constant": "Constant baseline",
          "P3_prop_shuffled": "Shuffled baseline", "P4_prop_half": "Baseline = 0.5", "P6_null_permuted_factor": "Null: permuted SpliceAI",
          "P0_raw": "SpliceAI raw"}
    for i, c in enumerate(order):
        r = met[(met.predictor == c) & (met.metric == "pearson")].iloc[0]
        col = C["fetal"] if c in ("P1_prop_measured", "P5_calibrated") else INK2
        ax.plot([r.ci_low, r.ci_high], [i, i], color=col, lw=2, solid_capstyle="round")
        ax.plot(r.estimate, i, "o", ms=6, color=col, mec=SURFACE, mew=1.5)
        ax.text(r.estimate, i - 0.28, f"{r.estimate:.3f}", fontsize=7, color=INK2, ha="center", va="bottom")
    ax.set_yticks(range(len(order)), [nm[c] for c in order], fontsize=8, color=INK)
    ax.set_ylim(len(order) - 0.5, -0.8)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_title("Pearson r with measured ΔPSI (held-out exons)", fontsize=9, color=INK, loc="left")
    fig.suptitle("Measured-effect test (MFASS minigene, HEK293T): does the true baseline improve predicted ΔPSI?",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE4_mfass.png", dpi=160, facecolor=SURFACE)


def fig_e5():
    ta = pd.read_csv(R / "E7_local_coverage" / "transcript_agreement.tsv", sep="\t")
    s = json.loads((R / "E7_local_coverage" / "summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), facecolor=SURFACE)
    ax = axes[0]
    style(ax, "y")
    bins = ["1-9", "10-99", ">=100"]
    pairs = [("d209~d334", "Patowary donor vs donor", C["patowary_fetal"]), ("leung_FetalCTX1~leung_FetalCTX2", "Leung fetal replicate 1 vs 2", C["fetal"]),
             ("fetal~adult", "Leung fetal vs adult", C["adult"])]
    x = np.arange(len(bins))
    for j, (pr, lab, col) in enumerate(pairs):
        vals = [ta[(ta.pair == pr) & (ta.abundance_bin_mean_donor_reads == b)].spearman.iloc[0] for b in bins]
        ax.bar(x + (j - 1) * 0.26, vals, width=0.24, color=col, label=lab)
    ax.set_xticks(x, [f"{b} reads" for b in bins], fontsize=8, color=INK)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Spearman of within-gene proportions", fontsize=8, color=INK2)
    ax.set_title("Agreement at the 1,670 cluster transcripts, by abundance\n(55–84% of transcripts have 0 reads: not shown, undefined)",
                 fontsize=9, color=INK, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax = axes[1]
    style(ax)
    src = ["fetal", "adult", "patowary", "d209", "d334", "d336"]
    lab = {"fetal": "Leung fetal", "adult": "Leung adult", "patowary": "Patowary pooled", "d209": "Patowary 209",
           "d334": "Patowary 334", "d336": "Patowary 336"}
    vals = [s["target_intron_unsupported_variants"][k] for k in src]
    ax.barh([lab[k] for k in src], vals, color=[C["fetal"], C["adult"], C["patowary_fetal"]] + [C["patowary_gz"]] * 3, height=0.6)
    for i, v in enumerate(vals):
        ax.text(v + 1, i, f"{v}/147", va="center", fontsize=7, color=INK2)
    ax.invert_yaxis()
    ax.set_xlim(0, 85)
    ax.set_title("Evaluated lead introns with no supporting read", fontsize=9, color=INK, loc="left")
    fig.suptitle("Baseline quality at the evaluated loci: genome-wide donor correlation does not reflect low-expression transcripts",
                 fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE5_local_coverage.png", dpi=160, facecolor=SURFACE)


def fig_e6():
    s = json.loads((R / "E9_sensitivity" / "summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), facecolor=SURFACE)
    ax = axes[0]
    style(ax)
    B = s["B_direction_accuracy_with_tolerance"]
    keys = list(B)
    for i, k in enumerate(keys):
        e, (lo, hi) = B[k]["fetal_minus_adult"], B[k]["ci"]
        ax.plot([lo, hi], [i, i], color=C["fetal"], lw=2, solid_capstyle="round")
        ax.plot(e, i, "o", ms=7, color=C["fetal"], mec=SURFACE, mew=1.5)
        ax.text(e, i - 0.25, f"{e:+.4f}", fontsize=7, color=INK2, ha="center", va="bottom")
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(range(len(keys)), ["v2 as reported (τ = 0)", "τ = 1e-12 (numerical zero)", "τ = 1e-6", "τ = 1e-3"], fontsize=8, color=INK)
    ax.set_ylim(len(keys) - 0.5, -0.8)
    ax.set_xlabel("Fetal − adult direction accuracy, Wen sQTL (95% CI)", fontsize=8, color=INK2)
    ax.set_title("Treating |pred| ≤ τ as a tie removes the v2 point estimate", fontsize=9, color=INK, loc="left")
    ax = axes[1]
    style(ax)
    P = s["C_oracle_power"]
    names = ["Direction accuracy", "Spearman ρ*", "Pearson on magnitudes*"]
    vals = [P["power_direction_accuracy_mcnemar"], P["power_spearman_bootstrap"], P["power_pearson_on_magnitudes_bootstrap"]]
    ax.barh(names, vals, color=[C["adult"], INK2, INK2], height=0.55)
    for i, v in enumerate(vals):
        ax.text(v + 0.01, i, f"{v:.0%}", va="center", fontsize=8, color=INK2)
    ax.set_xlim(0, 1.12)
    ax.invert_yaxis()
    ax.set_title("Power to detect a fetal baseline that is correct by construction\n(*assumes z linear in the predicted change)",
                 fontsize=9, color=INK, loc="left")
    fig.suptitle("The direction benchmark is structurally insensitive to the baseline", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE6_tolerance_power.png", dpi=160, facecolor=SURFACE)


def fig_e7():
    s = json.loads((R / "E11_vexseq" / "report.json").read_text())
    rs = json.loads((R / "E11_vexseq" / "replicate_split.json").read_text())
    v = pd.read_csv(R / "E11_vexseq" / "predictions.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.9), facecolor=SURFACE, gridspec_kw={"width_ratios": [1, 1, 1.1]})
    ax = axes[0]
    style(ax)
    names = [("M", "Matched baseline"), ("S", "Swapped (other cell line)"), ("C", "Constant baseline"),
             ("R", "Raw SpliceAI"), ("N", "Null (permuted factor)")]
    cols = [C["fetal"], C["adult"], C["uniform"], INK2, GRID]
    for i, ((k, lab), col) in enumerate(zip(names, cols)):
        e, (lo, hi) = s["point_all"][k], s["ci_all"][k]
        ax.plot([lo, hi], [i, i], color=col if k != "N" else INK2, lw=2, solid_capstyle="round")
        ax.plot(e, i, "o", ms=7, color=col if k != "N" else INK2, mec=SURFACE, mew=1.5)
        ax.text(hi + 0.02, i, f"{e:.3f}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(range(len(names)), [n for _, n in names], fontsize=8, color=INK)
    ax.set_ylim(len(names) - 0.5, -0.6)
    ax.set_xlim(-0.1, 0.9)
    ax.set_xlabel("Pearson r with measured ΔPSI (95% CI)", fontsize=8, color=INK2)
    ax.set_title("Predictor accuracy, K562 + HepG2 pooled", fontsize=9, color=INK, loc="left")
    ax = axes[1]
    style(ax)
    rows = [("M − S", "M - S (pearson)", "M - S"), ("M − C", "M - C (pearson)", "M - C"),
            ("Between-context\nM − null", "between-context: M - null N", "between M - N")]
    for i, (lab, k1, k2) in enumerate(rows):
        for j, (src, col, nm) in enumerate([(s["primary"][k1], C["fetal"], "E11"), (rs["all_exons"][k2], C["adult"], "replicate split")]):
            y = i + (j - 0.5) * 0.3
            ax.plot(src["ci"], [y, y], color=col, lw=2, solid_capstyle="round")
            ax.plot(src["estimate"], y, "o", ms=6, color=col, mec=SURFACE, mew=1.2, label=nm if i == 0 else None)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=8, color=INK)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlabel("Difference in Pearson r (95% CI)", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_title("Primary contrasts (pre-registered) and robustness", fontsize=9, color=INK, loc="left")
    ax = axes[2]
    style(ax, grid="both")
    big = (v.K562_base - v.HepG2_base).abs() >= 0.1
    px, mx = v.M_K562 - v.M_HepG2, v.K562_dpsi - v.HepG2_dpsi
    ax.scatter(px[~big], mx[~big], s=5, color=C["uniform"], alpha=0.35, lw=0, label="other exons")
    ax.scatter(px[big], mx[big], s=7, color=C["fetal"], alpha=0.7, lw=0, label="32 exons, baseline differs ≥ 0.1")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("Predicted ΔPSI(K562) − ΔPSI(HepG2), matched baselines", fontsize=8, color=INK2)
    ax.set_ylabel("Measured ΔPSI(K562) − ΔPSI(HepG2)", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.set_title(f"Between-context effect differences (r = {s['point_all']['between_M']:.2f}; null {s['point_all']['between_N']:.2f})",
                 fontsize=9, color=INK, loc="left")
    fig.suptitle("Vex-seq (two non-neural cell lines): a context-matched baseline improves quantitative ΔPSI prediction, modestly",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE7_vexseq.png", dpi=160, facecolor=SURFACE)


def fig_e8():
    a = json.loads((R / "E12_gse310131" / "summary_E12a.json").read_text())["local_coverage"]
    b = json.loads((R / "E12_gse310131" / "summary_E12b.json").read_text())["sensitivity"]["Wen_sQTL"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.0), facecolor=SURFACE, gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    style(ax)
    src = [("fetal", "Leung fetal cortex"), ("adult", "Leung adult cortex"), ("patowary", "Patowary fetal cortex"),
           ("gse_t00", "GSE310131 t00 iPSC"), ("gse_t04", "GSE310131 t04 intermediate"), ("gse_t30", "GSE310131 t30 NGN2 neuron")]
    cols = [C["fetal"], C["adult"], C["patowary_fetal"], INK2, INK2, INK2]
    vals = [a[f"lead_intron_unsupported_transcript_level_{k}"] for k, _ in src]
    ax.barh([n for _, n in src], vals, color=cols, height=0.6)
    for i, v in enumerate(vals):
        ax.text(v + 1, i, str(v), va="center", fontsize=8, color=INK2)
    for i, st in enumerate(["t00", "t04", "t30"]):
        j = a[f"lead_intron_zero_junction_reads_{st}"]
        ax.plot(j, 3 + i, "D", ms=6, color=C["generic"], mec=SURFACE, label="junction level (all GSE isoforms)" if i == 0 else None)
    ax.invert_yaxis()
    ax.set_xlim(0, 115)
    ax.set_xlabel("Lead introns with no supporting read (of 147 Wen sQTL positives)", fontsize=8, color=INK2)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_title("Local support at the evaluated loci", fontsize=9, color=INK, loc="left")
    ax = axes[1]
    style(ax)
    pairs = [("gse_t30_L1 vs gse_t30_L2", "t30 lineage 1 vs 2 (same stage)"), ("gse_t30_L1 vs gse_t30_L3", "t30 lineage 1 vs 3"),
             ("gse_t30 vs gse_t04", "t30 vs t04 (within study)"), ("gse_t30 vs gse_t00", "t30 vs t00 (within study)"),
             ("fetal vs adult", "Leung fetal vs adult"), ("gse_t30 vs patowary_fetal", "t30 vs Patowary fetal (cross-study)"),
             ("gse_t30 vs fetal", "t30 vs Leung fetal (cross-study)"), ("gse_t00 vs fetal", "t00 vs Leung fetal (cross-study)")]
    v = [b[k]["frac_ratio_outside_[0.5,2]"] for k, _ in pairs]
    fl = [b[k]["tau=1e-06"]["sign_flips"] for k, _ in pairs]
    nb = [b[k]["tau=1e-06"]["both_above"] for k, _ in pairs]
    cc = [GRID, GRID, C["uniform"], C["uniform"], C["fetal"], INK2, INK2, INK2]
    ax.barh([n for _, n in pairs], v, color=cc, height=0.6)
    for i, (x, f_, n_) in enumerate(zip(v, fl, nb)):
        ax.text(x + 0.01, i, f"{x:.0%}   sign flips {f_}/{n_}", va="center", fontsize=7.5, color=INK2)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.75)
    ax.set_xlabel("Share of variants whose predicted |Δu| changes more than 2-fold (Wen sQTL)", fontsize=8, color=INK2)
    ax.set_title("Prediction sensitivity: study/platform outweighs differentiation stage", fontsize=9, color=INK, loc="left")
    fig.suptitle("GSE310131 (one donor, iPSC → NGN2 neuron): deeper local support; stage changes magnitudes modestly, never direction",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(F / "figE8_gse310131.png", dpi=160, facecolor=SURFACE)


if __name__ == "__main__":
    F.mkdir(parents=True, exist_ok=True)
    fig_e8()
    fig_e7()
    fig_e1()
    fig_e2()
    fig_e3()
    fig_e4()
    fig_e5()
    fig_e6()
    print("figures in", F)

"""Manuscript Figures 2, 3c, 4, 5 and 6 (research_documentation/08_FIGURES_AND_TABLES_PLAN.md).

Every plotted number is read from a results file produced by an earlier script; nothing is
re-estimated here, except Spearman correlations and simple shares computed directly from
predictions.tsv (outcome-blind). Outputs go to publication_extension/figures/paper/ as PDF and PNG,
plus figure_values.json with every number shown, for traceability.

Run: python -E publication_extension/scripts/make_paper_figures.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "publication_extension/results"
OUT = ROOT / "publication_extension/figures/paper"
OUT.mkdir(parents=True, exist_ok=True)

SRC = {
    "e15_results": RES / "E15_devstage/E15_results.json",
    "e15_predictions": RES / "E15_devstage/predictions.tsv",
    "e15_diagnostics": RES / "E15_devstage/baseline_diagnostics.json",
    "e15_gate0_mei": RES / "E15_devstage/gate0_mei.json",
    "e2_summary": RES / "E2_root_cause/summary.json",
    "e13_decomposition": RES / "E13_compass/context_decomposition.json",
    "e13_headroom": RES / "E13_compass/headroom_scale.json",
    "e13_generalization": RES / "E13_compass/generalization.json",
    "e14_comparators": RES / "E14_comparators/comparators.json",
}

# Reference categorical palette (dataviz skill, light mode), fixed order; U/R use neutral ink.
C = {"F1": "#2a78d6", "A1": "#eb6834", "G": "#1baf7a", "U": "#8a8984"}
CM = {"A": "#2a78d6", "K": "#eb6834", "B1": "#1baf7a", "R": "#8a8984"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
VALUES = {"sources": {k: str(v.relative_to(ROOT)) for k, v in SRC.items()}}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 150,
})


def load(k):
    return json.loads(SRC[k].read_text())


def panel(ax, letter, dx=-0.14):
    ax.text(dx, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def hgrid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


# ------------------------------------------------------------------ Figure 2 (outcome-blind)
def fig2():
    pr = pd.read_csv(SRC["e15_predictions"], sep="\t")
    diag = load("e15_diagnostics")["F1_stability"]
    e2 = load("e2_summary")
    cl = pr.drop_duplicates("cluster")["cluster_max_abs_usage_diff_F1_A1"].dropna().to_numpy()
    rho_pred = float(spearmanr(pr.pred_F1, pr.pred_A1).statistic)
    elig = diag["F1_vs_A1_cluster_max_abs_usage_difference_eligible"]
    dist = diag["frac_eligible_clusters_F1_A1_distinguishable(|du|>0.05_and_>2_bootstrap_SE)"]
    v = VALUES["fig2"] = {
        "a_n_analysed_clusters": int(len(cl)),
        "a_share_analysed_clusters_max_abs_du_gt_0.05": float((cl > 0.05).mean()),
        "a_median_analysed_clusters": float(np.median(cl)),
        "a_eligible_clusters_n": elig["n"], "a_eligible_median": elig["p50"],
        "a_eligible_frac_distinguishable": dist,
        "b_n_units": int(len(pr)), "b_spearman_pred_F1_A1": rho_pred,
        "c": {k: {"n": e2[k]["n_positive"], "sign_invariant": e2[k]["sign_invariant_by_structure"],
                  "flippable": e2[k]["sign_flippable_in_principle"],
                  "observed_flips_fetal_vs_adult": e2[k]["observed_sign_flips_fetal_vs_adult"]}
              for k in ("sQTL", "isoQTL")},
    }

    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"width_ratios": [1.15, 1, 0.9]})
    ax = axs[0]
    bins = np.logspace(-4, 0, 41)
    ax.hist(np.clip(cl, 1e-4, 1), bins=bins, color=C["F1"], edgecolor="white", lw=0.4)
    ax.set_xscale("log")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.75)
    ax.axvline(0.05, color=INK, lw=0.8, ls="--")
    ax.text(0.058, ax.get_ylim()[1] * 0.97, "0.05", fontsize=7, color=INK, va="top")
    ax.set_xlabel("Cluster max |u$_{F1}$ − u$_{A1}$|")
    ax.set_ylabel("Analysed clusters")
    ax.set_title("Baselines differ", loc="left")
    ax.text(0.02, 0.97, f"{v['a_share_analysed_clusters_max_abs_du_gt_0.05']:.0%} of {len(cl):,} analysed\n"
            f"clusters > 0.05\n{dist:.0%} of {elig['n']:,} eligible\nclusters distinguishable*",
            transform=ax.transAxes, va="top", fontsize=6.8, color=INK2)
    hgrid(ax)
    panel(ax, "a")

    ax = axs[1]
    ax.scatter(pr.pred_A1, pr.pred_F1, s=6, color=C["F1"], alpha=0.45, lw=0)
    lim = np.nanpercentile(np.abs(np.r_[pr.pred_F1, pr.pred_A1]), 99.5)
    ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=0.8, ls="--")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    v["b_axis_limit_abs"] = float(lim)
    v["b_units_outside_axes"] = int(((pr.pred_F1.abs() > lim) | (pr.pred_A1.abs() > lim)).sum())
    ax.set_xlabel("Predicted Δusage, adult baseline (A1)")
    ax.set_ylabel("Predicted Δusage, fetal (F1)")
    ax.set_title("Predictions barely differ", loc="left")
    ax.text(0.04, 0.96, f"Spearman ρ = {rho_pred:.3f}\nn = {len(pr):,} units", transform=ax.transAxes,
            va="top", fontsize=7, color=INK)
    panel(ax, "b")

    ax = axs[2]
    names = ["sQTL", "isoQTL"]
    inv = np.array([v["c"][k]["sign_invariant"] / v["c"][k]["n"] for k in names])
    ax.bar(names, inv, color=C["F1"], width=0.55, label="Sign fixed by structure")
    ax.bar(names, 1 - inv, bottom=inv, color=GRID, width=0.55, edgecolor="white", lw=1,
           label="Sign can depend on baseline")
    for i, k in enumerate(names):
        ax.text(i, inv[i] / 2, f"{inv[i]:.0%}", ha="center", va="center", color="white", fontsize=7.5,
                fontweight="bold")
        ax.text(i, 1.02, f"n={v['c'][k]['n']}\n{v['c'][k]['observed_flips_fetal_vs_adult']} observed\nF1/A1 flips",
                ha="center", va="bottom", fontsize=6.3, color=INK2)
    ax.set_ylim(0, 1.32)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0%", "50%", "100%"])
    ax.set_ylabel("Positive v2 variants")
    ax.set_title("Two-isoform sign (v2)", loc="left")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), frameon=False, ncol=1, fontsize=6.5)
    panel(ax, "c")
    fig.tight_layout(w_pad=1.6)
    save(fig, "Fig2_baseline_vs_prediction_differences")


# ------------------------------------------------------------------ Figure 3c
def fig3c():
    d = load("e13_decomposition")
    h = load("e13_headroom")
    S, Sci = d["point"]["S_share_explained"], d["ci"]["S_share_explained"]
    H, Hci = d["point"]["H_headroom"], d["ci"]["H_headroom"]
    expl = h["explained_by_baseline_only_as_frac_of_var_dPSI"]
    unex = h["unexplained_as_frac_of_var_dPSI"]
    VALUES["fig3c"] = {"S": S, "S_ci": Sci, "H": H, "H_ci": Hci, "explained_frac_var_dPSI": expl,
                       "unexplained_frac_var_dPSI": unex, "reliability": h["reliability"]}

    fig, axs = plt.subplots(1, 2, figsize=(5.6, 2.1), gridspec_kw={"width_ratios": [1, 1]})
    ax = axs[0]
    ys = [1, 0]
    for y, (lab, val, ci, col) in zip(ys, [("Explained by\nreference PSI (S)", S, Sci, CM["A"]),
                                          ("Residual\nheadroom (H)", H, Hci, CM["R"])]):
        ax.barh(y, val, color=col, height=0.5)
        ax.errorbar(val, y, xerr=[[val - ci[0]], [ci[1] - val]], color=INK, lw=0.9, capsize=2.5)
        ax.text(ci[1] + 0.02, y, f"{val:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]", va="center", fontsize=7, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels(["Explained by\nreference PSI (S)", "Residual\nheadroom (H)"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of reliable between-cell ΔPSI difference")
    ax.set_title("Context-specificity (COMPASS, 5 lines)", loc="left")
    hgrid(ax, "x")
    panel(ax, "c")

    ax = axs[1]
    ax.barh([1, 0], [expl * 100, unex * 100], color=[CM["A"], CM["R"]], height=0.5)
    for y, val in zip([1, 0], [expl, unex]):
        ax.text(val * 100 + 0.2, y, f"{val:.1%}", va="center", fontsize=7, color=INK)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["Explained by\nreference PSI", "Headroom for\ncontext models"])
    ax.set_xlim(0, 10)
    ax.set_xlabel("% of total ΔPSI variance")
    ax.set_title("Absolute scale", loc="left")
    hgrid(ax, "x")
    fig.tight_layout(w_pad=2.0)
    save(fig, "Fig3c_compass_context_decomposition")


# ------------------------------------------------------------------ Figure 4
def fig4():
    g = load("e13_generalization")
    sets = [("test1_compass_unseen_exons_genes", "COMPASS test\n(never fitted)"),
            ("test2_vexseq", "Vex-seq\n(seen earlier)"), ("test3_mfass", "MFASS\n(seen earlier)")]
    models = ["A", "K", "B1", "R"]
    labels = {"A": "A (propagation, no parameters)", "K": "K (global logit, 1 parameter)",
              "B1": "B1 (logit, 3 parameters)", "R": "R (raw SpliceAI)"}
    VALUES["fig4"] = {lab.split("\n")[0]: {"n_obs": g[k]["n_obs"],
                                          "pearson": {m: g[k]["point"][m]["pearson"] for m in models},
                                          "ci": {m: g[k]["ci_pearson"][m] for m in models},
                                          "B1-A": g[k]["contrasts"]["B1 - A"], "K-A": g[k]["contrasts"]["K - A"]}
                      for k, lab in sets}
    VALUES["fig4"]["gate_G-B"] = g["gate_G-B"]

    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    w = 0.19
    for j, m in enumerate(models):
        xs = np.arange(len(sets)) + (j - 1.5) * w
        pts = [g[k]["point"][m]["pearson"] for k, _ in sets]
        ci = np.array([g[k]["ci_pearson"][m] for k, _ in sets])
        ax.errorbar(xs, pts, yerr=[np.array(pts) - ci[:, 0], ci[:, 1] - np.array(pts)], fmt="o",
                    color=CM[m], ms=4.5, lw=1, capsize=2, label=labels[m])
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([lab for _, lab in sets])
    ax.set_ylabel("Pearson r (predicted vs measured ΔPSI)")
    ax.set_ylim(0.2, 0.85)
    ax.set_title("Accuracy by test assay", loc="left")
    ax.legend(loc="lower left", frameon=False, fontsize=6.5)
    hgrid(ax)
    panel(ax, "a", dx=-0.2)

    ax = axs[1]
    rows = []
    for k, lab in sets:
        for c, col in (("B1 - A", CM["B1"]), ("K - A", CM["K"])):
            rows.append((lab.split("\n")[0], c, g[k]["contrasts"][c], col))
    ys = np.arange(len(rows))[::-1]
    for y, (lab, c, r, col) in zip(ys, rows):
        ax.errorbar(r["estimate"], y, xerr=[[r["estimate"] - r["ci"][0]], [r["ci"][1] - r["estimate"]]],
                    fmt="s", color=col, ms=4, lw=1, capsize=2)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{lab}: {c.replace(' - ', '−')}" for lab, c, _, _ in rows], fontsize=6.8)
    ax.set_xlabel("Δ Pearson r (calibrated − A), 95% CI")
    ax.set_title("Calibration does not transfer", loc="left")
    hgrid(ax, "x")
    panel(ax, "b")
    fig.tight_layout(w_pad=1.5)
    save(fig, "Fig4_calibration_generalization")


# ------------------------------------------------------------------ Figure 5 (primary E15 result)
def f1a1(entry):
    """The F1-versus-A1 contrast of an E15 result block, whatever its key prefix (x, e)."""
    for k, v in entry["contrasts"].items():
        a, b = k.split("-")
        if a.lstrip("xe").startswith("F") and b.lstrip("xe").startswith("A"):
            return v
    raise KeyError("no F1-A1 contrast")


def fig5():
    d = load("e15_results")
    mei = load("e15_gate0_mei")["MEI"]
    p = d["primary"]
    bl = ["F1", "A1", "G", "U"]
    names = {"F1": "Fetal (F1)", "A1": "Adult (A1)", "G": "Generic GTEx (G)", "U": "Uniform (U)"}
    sens = d["sensitivity"]
    srows = [("Primary (rule P)", p),
             ("Rule D (closest to splice site)", sens["rule_D"]),
             ("Strong variants (DS ≥ 0.10)", sens["strong_DSmax>=0.10"]),
             ("Weak variants (DS < 0.10)", sens["weak_DSmax<0.10"]),
             ("One unit per signal group", sens["one_unit_per_signal_group"]),
             ("F1 incl. HSB155 (13 donors)", sens["F1_all"]),
             ("Standard regions only", sens["standard_regions"]),
             ("Spliced fraction ≥ 0.10", sens["spliced_fraction>=0.10"]),
             ("No thinning", sens["no_thinning"]),
             ("ε = 0.001", sens["eps=0.001"]),
             ("ε = 0.1", sens["eps=0.1"]),
             ("Truth: conditional slope", sens["secondary_truth_conditional_bwd_slope"]),
             ("Baseline difference T1 (low)", sens["by_baseline_difference_tertile"]["T1"]),
             ("Baseline difference T2", sens["by_baseline_difference_tertile"]["T2"]),
             ("Baseline difference T3 (high)", sens["by_baseline_difference_tertile"]["T3"]),
             ("Adult truth (GTEx cortex)*", d["crossover"]["adult_truth"])]
    VALUES["fig5"] = {"MEI": mei, "n_units": p["n_units"], "n_components": p["n_components"],
                      "rho": p["rho"], "rho_95CI": p["rho_95CI"], "contrasts": p["contrasts"],
                      "sensitivity_F1-A1": {lab: {"n_units": e["n_units"], **f1a1(e)} for lab, e in srows}}

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[0.8, 1.35], hspace=0.42, wspace=0.42)
    ax = fig.add_subplot(gs[0, 0])
    for i, b in enumerate(bl):
        r, ci = p["rho"][b], p["rho_95CI"][b]
        ax.errorbar(i, r, yerr=[[r - ci[0]], [ci[1] - r]], fmt="o", color=C[b], ms=5, lw=1.1, capsize=2.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["F1", "A1", "G", "U"])
    ax.set_ylabel("Spearman ρ (pred vs SuSiE z)")
    ax.set_ylim(0.55, 0.75)
    ax.set_title(f"Agreement by baseline (n = {p['n_units']:,})", loc="left")
    hgrid(ax)
    panel(ax, "a")

    ax = fig.add_subplot(gs[0, 1])
    cons = ["F1-A1", "F1-G", "F1-U"]
    ys = np.arange(len(cons))[::-1]
    ax.axvspan(-mei, mei, color=GRID, lw=0)
    ax.axvline(0, color=INK2, lw=0.8)
    for y, c in zip(ys, cons):
        r = p["contrasts"][c]
        ax.plot(r["ci_holm_level"], [y, y], color=INK2, lw=0.9)
        ax.errorbar(r["delta_rho"], y, xerr=[[r["delta_rho"] - r["ci95"][0]], [r["ci95"][1] - r["delta_rho"]]],
                    fmt="o", color=C["F1"], ms=5, lw=2.2, capsize=0)
        ax.text(0.047, y, r["decision"], va="center", fontsize=6.5, color=INK2)
    ax.set_yticks(ys)
    ax.set_yticklabels(["F1 − A1", "F1 − G", "F1 − U"])
    ax.set_xlim(-0.05, 0.10)
    ax.set_xticks([-0.04, -0.02, 0, 0.02, 0.04])
    ax.set_xlabel("Δρ (thick: 95% CI; thin: Holm-level CI)")
    ax.set_title(f"Primary contrasts (grey band: ±MEI {mei})", loc="left")
    panel(ax, "b")

    ax = fig.add_subplot(gs[1, :])
    ys = np.arange(len(srows))[::-1]
    ax.axvspan(-mei, mei, color=GRID, lw=0)
    ax.axvline(0, color=INK2, lw=0.8)
    for y, (lab, e) in zip(ys, srows):
        r = f1a1(e)
        col = C["F1"] if lab.startswith("Primary") else INK2
        ax.errorbar(r["delta_rho"], y, xerr=[[r["delta_rho"] - r["ci95"][0]], [r["ci95"][1] - r["delta_rho"]]],
                    fmt="o", color=col, ms=4, lw=1.2, capsize=1.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{lab} (n={e['n_units']:,})" for lab, e in srows], fontsize=7)
    ax.set_xlim(-0.03, 0.045)
    ax.text(mei, len(srows) - 0.35, " +MEI", fontsize=6.5, color=INK2, va="bottom")
    ax.text(-mei, len(srows) - 0.35, "−MEI ", fontsize=6.5, color=INK2, va="bottom", ha="right")
    ax.set_xlabel("Δρ F1 − A1, 95% CI")
    ax.set_title("F1 − A1 under pre-specified sensitivity analyses", loc="left")
    hgrid(ax, "x")
    panel(ax, "c")
    save(fig, "Fig5_E15_primary_result")


def fig6():
    """E14 external comparators (single frozen E14c run)."""
    d = load("e14_comparators")
    T = d["T-fresh"]
    S = T["N_DS_subsample"]
    kind = {"A": CM["A"], "R": CM["R"], "N": CM["K"], "F": CM["B1"]}
    rows = [("A_SA", "A", "SpliceAI, rule A"), ("A_PG", "A", "Pangolin, rule A"), ("A_DS", "A", "DeltaSplice, rule A"),
            ("Rr_SA", "R", "SpliceAI, raw"), ("Rr_PG", "R", "Pangolin, raw (native)"), ("Rr_DS", "R", "DeltaSplice, raw"),
            ("N_MM", "N", "MMSplice, native"), ("K", "F", "K (1 fitted param.)"), ("B1", "F", "B1 (3 fitted param.)")]
    p1 = ["A_SA - N_MM", "A_SA - Rr_PG", "A_SA - N_DS"]
    p2 = ["A_SA - Rr_SA", "A_PG - Rr_PG", "A_DS - Rr_DS"]
    VALUES["fig6"] = {"n_obs": T["n_obs"], "n_families": T["n_families"], "sub_n_obs": S["n_obs"],
                      "sub_n_families": S["n_families"], "pearson": {m: T["point"][m]["pearson"] for m, _, _ in rows},
                      "pearson_ci": {m: T["ci_pearson"][m] for m, _, _ in rows},
                      "N_DS_sub_pearson": S["point"]["N_DS"]["pearson"], "A_SA_sub_pearson": S["point"]["A_SA"]["pearson"],
                      "contrasts": {c: {k: T["contrasts"][c][k] for k in ("estimate", "ci95", "ci_holm", "p_holm")}
                                    for c in p1 + p2},
                      "gate_G-X": d["gate_G-X"]}

    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.1], height_ratios=[1.1, 1], hspace=0.55, wspace=0.75)

    ax = fig.add_subplot(gs[0, 0])
    ys = np.arange(len(rows) + 1)[::-1]
    for y, (m, k, lab) in zip(ys, rows):
        r, ci = T["point"][m]["pearson"], T["ci_pearson"][m]
        ax.errorbar(r, y, xerr=[[r - ci[0]], [ci[1] - r]], fmt="o", color=kind[k], ms=4.5, lw=1.2, capsize=0)
    r, ci = S["point"]["N_DS"]["pearson"], S["ci_pearson"]["N_DS"]
    ax.errorbar(r, ys[-1], xerr=[[r - ci[0]], [ci[1] - r]], fmt="D", color=kind["N"], ms=4, lw=1.2, capsize=0,
                mfc="white")
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, _, lab in rows] + [f"DeltaSplice, native*"], fontsize=7)
    ax.set_xlabel("Pearson r with measured ΔPSI (95% CI)")
    ax.set_xlim(0.15, 0.9)
    ax.set_title(f"T-fresh, {T['n_obs']:,} obs", loc="left")
    hgrid(ax, "x")
    panel(ax, "a", dx=-0.62)

    ax = fig.add_subplot(gs[0, 1])
    labs = {"A_SA - N_MM": "A_SA − MMSplice", "A_SA - Rr_PG": "A_SA − Pangolin",
            "A_SA - N_DS": "A_SA − DeltaSplice native*", "A_SA - Rr_SA": "A_SA − R_SA",
            "A_PG - Rr_PG": "A_PG − R_PG", "A_DS - Rr_DS": "A_DS − R_DS"}
    cons = p1 + p2
    ys = np.array([6.5, 5.5, 4.5, 2.5, 1.5, 0.5])
    ax.axvline(0, color=INK2, lw=0.8)
    ax.axvline(0.02, color=INK2, lw=0.8, ls=(0, (3, 2)))
    for y, c in zip(ys, cons):
        e = T["contrasts"][c]
        ax.plot(e["ci_holm"], [y, y], color=INK2, lw=0.9)
        ax.errorbar(e["estimate"], y, xerr=[[e["estimate"] - e["ci95"][0]], [e["ci95"][1] - e["estimate"]]],
                    fmt="o", color=CM["A"], ms=4.5, lw=2.2, capsize=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([labs[c] for c in cons], fontsize=7)
    ax.text(-0.155, 7.2, "P1 (gate G-X)", fontsize=7, color=INK2, va="bottom")
    ax.text(-0.155, 3.2, "P2 (rule A vs raw)", fontsize=7, color=INK2, va="bottom")
    ax.text(0.021, -0.35, "G-X margin 0.02", fontsize=6.5, color=INK2, va="bottom")
    ax.set_ylim(-0.4, 7.9)
    ax.set_xlim(-0.16, 0.22)
    ax.set_xlabel("Δr (thick: 95% CI; thin: Holm-level CI)")
    ax.set_title("Primary contrasts on T-fresh", loc="left")
    hgrid(ax, "x")
    panel(ax, "b", dx=-0.62)

    ax = fig.add_subplot(gs[1, :])
    sets = [("T-E13 (previously seen)", "COMPASS test (E13)"),
            ("T-Vex (previously seen; MMSplice fitted on Vex-seq: not interpreted)", "Vex-seq"),
            ("T-MFASS (previously seen)", "MFASS")]
    models = [("A_SA", "A", "SpliceAI rule A"), ("Rr_SA", "R", "SpliceAI raw"), ("Rr_PG", "R", "Pangolin raw"),
              ("N_MM", "N", "MMSplice native")]
    VALUES["fig6"]["secondary"] = {}
    for i, (key, name) in enumerate(sets):
        e = d[key]
        VALUES["fig6"]["secondary"][name] = {"n_obs": e["n_obs"], "n_families": e["n_families"],
                                            "pearson": {m: e["point"][m]["pearson"] for m, _, _ in models},
                                            "ci": {m: e["ci_pearson"][m] for m, _, _ in models}}
        for j, (m, k, lab) in enumerate(models):
            x = i + (j - 1.5) * 0.17
            r, ci = e["point"][m]["pearson"], e["ci_pearson"][m]
            hollow = name == "Vex-seq" and m == "N_MM"
            ax.errorbar(x, r, yerr=[[r - ci[0]], [ci[1] - r]], fmt="s" if m.startswith("Rr_PG") else "o",
                        color=kind[k], ms=4.5, lw=1.2, capsize=0, mfc="white" if hollow else kind[k],
                        label=lab if i == 0 else None)
    ax.set_xticks(range(3))
    ax.set_xticklabels([n + "\n" + f"{d[k]['n_obs']:,} obs, {d[k]['n_families']:,} families" for k, n in sets])
    ax.set_ylabel("Pearson r (95% CI)")
    ax.set_ylim(0.2, 0.86)
    ax.set_title("Secondary sets, previously used in earlier experiments", loc="left")
    ax.legend(ncol=4, frameon=False, loc="upper right", bbox_to_anchor=(1, 1.02))
    hgrid(ax)
    panel(ax, "c", dx=-0.07)
    save(fig, "Fig6_E14_external_comparators")


if __name__ == "__main__":
    fig2()
    fig3c()
    fig4()
    fig5()
    fig6()
    (OUT / "figure_values.json").write_text(json.dumps(VALUES, indent=1, default=float))
    print("written to", OUT)

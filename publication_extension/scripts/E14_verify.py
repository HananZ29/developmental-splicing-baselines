"""Independent verification of E14 comparators.json (added 2026-11-08, after the single frozen run).

This does not change the frozen E14 result. It checks:
  1. the rule-A and N_MM predictions recompute from the raw tool features in predictions_fresh.tsv;
  2. every T-fresh point Pearson r (full set and N_DS subsample) recomputes with scipy;
  3. the primary contrasts, re-bootstrapped with separately written code (exon-family resampling,
     B = 1,000, seed 20261108, not the frozen 20261106), give the same gate G-X and P2 conclusions;
     intervals are expected to agree only within Monte Carlo error.
Writes publication_extension/results/E14_comparators/E14_verification.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_extension/results/E14_comparators"
EPS = 0.01
B, SEED = 1000, 20261108


def rule_a(p, la, ld):
    a = np.exp(la + ld)
    return p * a / (p * a + 1 - p) - p


def boot(df, preds, contrasts, rng):
    ex = df.exon.to_numpy()
    groups = [np.flatnonzero(ex == e) for e in pd.unique(ex)]
    y = df.y.to_numpy()
    P = {m: df[m].to_numpy() for m in preds}
    out = {f"{a} - {b}": np.empty(B) for a, b in contrasts}
    for i in range(B):
        pick = rng.integers(0, len(groups), len(groups))
        s = np.concatenate([groups[j] for j in pick])
        r = {m: np.corrcoef(P[m][s], y[s])[0, 1] for m in preds}
        for a, b in contrasts:
            out[f"{a} - {b}"][i] = r[a] - r[b]
    return out


def holm_ci(draws, names):
    p = {n: min(1.0, 2 * min((draws[n] <= 0).mean(), (draws[n] >= 0).mean())) for n in names}
    res = {}
    for k, n in enumerate(sorted(names, key=p.get)):
        lev = 0.05 / (len(names) - k)
        res[n] = {"p_boot": p[n], "ci_holm": [float(np.quantile(draws[n], lev / 2)),
                                              float(np.quantile(draws[n], 1 - lev / 2))]}
    return res


def main():
    src = OUT / "predictions_fresh.tsv"   # not version-controlled (>5 MB); the .gz copy is
    F = pd.read_csv(src if src.exists() else OUT / "predictions_fresh.tsv.gz", sep="\t")
    J = json.loads((OUT / "comparators.json").read_text())["T-fresh"]
    rep = {"label": "Independent verification, added after the single frozen run; does not alter it."}

    p = F.p.to_numpy()
    pc = np.clip(p, EPS, 1 - EPS)
    rec = {"A_SA": rule_a(pc, F.la_SA, F.ld_SA), "A_PG": rule_a(pc, F.la_PG, F.ld_PG),
           "A_DS": rule_a(pc, F.la_DS, F.ld_DS), "N_MM": expit(logit(pc) + F.mmsplice_dlogit_psi) - pc}
    rep["prediction_max_abs_diff"] = {k: float(np.nanmax(np.abs(np.asarray(v) - F[k]))) for k, v in rec.items()}

    preds = ["A_SA", "A_PG", "A_DS", "Rr_SA", "Rr_PG", "Rr_DS", "N_MM", "K", "B1"]
    S = F[np.isfinite(F.N_DS)]
    diffs = {}
    for m in preds:
        diffs[f"full:{m}"] = abs(pearsonr(F[m], F.y).statistic - J["point"][m]["pearson"])
    for m in preds + ["N_DS"]:
        diffs[f"sub:{m}"] = abs(pearsonr(S[m], S.y).statistic - J["N_DS_subsample"]["point"][m]["pearson"])
    rep["pearson_max_abs_diff"] = float(max(diffs.values()))
    rep["n_rows"] = {"full": len(F), "N_DS_subsample": len(S), "families_full": int(F.exon.nunique()),
                     "families_sub": int(S.exon.nunique())}

    rng = np.random.default_rng(SEED)
    p1 = [("A_SA", "N_MM"), ("A_SA", "Rr_PG")]
    p2 = [("A_SA", "Rr_SA"), ("A_PG", "Rr_PG"), ("A_DS", "Rr_DS")]
    d = boot(F, preds, p1 + p2, rng)
    d.update(boot(S, ["A_SA", "N_DS"], [("A_SA", "N_DS")], rng))
    P1 = ["A_SA - N_MM", "A_SA - Rr_PG", "A_SA - N_DS"]
    P2 = [f"{a} - {b}" for a, b in p2]
    h1, h2 = holm_ci(d, P1), holm_ci(d, P2)
    est = {n: J["contrasts"][n]["estimate"] for n in P1 + P2}
    gx = {n: bool(est[n] >= 0.02 and h1[n]["ci_holm"][0] > 0) for n in P1}
    better = {n: bool(h1[n]["ci_holm"][1] < 0) for n in P1}
    rep["rebootstrap"] = {"B": B, "seed": SEED, "P1": h1, "P2": h2,
                          "gate_G-X_pass": all(gx.values()), "per_comparator_pass": gx,
                          "comparator_significantly_better": better,
                          "P2_significant": {n: bool(h2[n]["ci_holm"][0] > 0) for n in P2}}
    frozen = json.loads((OUT / "comparators.json").read_text())["gate_G-X"]
    rep["agrees_with_frozen"] = {
        "gate_pass": frozen["pass"] == rep["rebootstrap"]["gate_G-X_pass"],
        "per_comparator_pass": frozen["per_comparator_pass"] == gx,
        "comparator_significantly_better": frozen["comparator_significantly_better"] == better,
        "P2_significant": {n: (J["contrasts"][n]["ci_holm"][0] > 0) == rep["rebootstrap"]["P2_significant"][n]
                           for n in P2},
    }
    (OUT / "E14_verification.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()

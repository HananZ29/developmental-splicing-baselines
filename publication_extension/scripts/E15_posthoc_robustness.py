"""E15 POST-HOC robustness checks (not pre-registered; added 2026-11-07, after E15_results.json existed).

These checks do not alter the frozen E15 decision, and they write only to a separate file:
publication_extension/results/E15_devstage/posthoc_robustness.json.

Check 1 - independent re-implementation of the frozen bootstrap.
    The frozen analysis (E15_stats.py) uses a vectorised multinomial-weight bootstrap with weighted
    Pearson correlation of fixed full-sample ranks. Here the same estimator is coded differently:
    - truth is re-read by column name with pandas;
    - each replicate explicitly resamples components with rng.choice and concatenates their units;
    - np.corrcoef is applied to the fixed full-sample ranks of the resampled rows.
    Agreement within Monte Carlo error checks the implementation, not the estimator.

Check 2 - conventional re-ranked component bootstrap.
    The same component resampling, but Spearman rho is recomputed with re-ranking inside each
    replicate (scipy.stats.spearmanr on the resampled rows). This checks whether the frozen
    fixed-rank simplification changes the interval.

Both checks use B = 10,000 and a seed (20261107) different from the frozen analysis (20261106), and
cover the three primary contrasts on the primary population (rule P, component_P).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "publication_extension/results/E15_devstage"
SUSIE = ROOT / "data/external/wen_synapse/mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz"
B = 10_000
SEED = 20261107
BASELINES = ["F1", "A1", "G", "U"]
CONTRASTS = [("F1", "A1"), ("F1", "G"), ("F1", "U")]
MEI = 0.0228


def load_units():
    pr = pd.read_csv(RES / "predictions.tsv", sep="\t")
    z = pd.read_csv(SUSIE, sep="\t", usecols=["phenotype_id", "variant_id", "cs_id", "z"])
    z = z.drop_duplicates(["phenotype_id", "cs_id", "variant_id"])
    d = pr.merge(z, left_on=["phenotype_id", "cs_id", "variant_id_P"],
                 right_on=["phenotype_id", "cs_id", "variant_id"], how="left", validate="many_to_one")
    assert len(d) == len(pr)
    d = d.dropna(subset=["z"] + [f"pred_{b}" for b in BASELINES])
    return d.reset_index(drop=True)


def summarise(boot, point):
    out = {}
    for a, b in CONTRASTS:
        k = f"{a}-{b}"
        v = np.asarray(boot[k])
        out[k] = {"delta_rho": float(point[k]),
                  "ci95": [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))],
                  "p_boot_two_sided": float(min(1.0, 2 * min((v <= 0).mean(), (v >= 0).mean()))),
                  "ci95_within_pm_MEI": bool(np.quantile(v, 0.025) > -MEI and np.quantile(v, 0.975) < MEI)}
    return out


def main():
    d = load_units()
    comp = pd.factorize(d.component_P)[0]
    members = [np.flatnonzero(comp == c) for c in range(comp.max() + 1)]
    K = len(members)
    y = d.z.to_numpy()
    P = {b: d[f"pred_{b}"].to_numpy() for b in BASELINES}
    ry = rankdata(y)
    RP = {b: rankdata(P[b]) for b in BASELINES}
    point = {f"{a}-{b}": spearmanr(y, P[a]).statistic - spearmanr(y, P[b]).statistic for a, b in CONTRASTS}
    rho = {b: float(spearmanr(y, P[b]).statistic) for b in BASELINES}

    rng = np.random.default_rng(SEED)
    fixed = {f"{a}-{b}": [] for a, b in CONTRASTS}
    rerank = {f"{a}-{b}": [] for a, b in CONTRASTS}
    for _ in range(B):
        idx = np.concatenate([members[c] for c in rng.choice(K, size=K, replace=True)])
        rf = {b: np.corrcoef(ry[idx], RP[b][idx])[0, 1] for b in BASELINES}
        rr = {b: spearmanr(y[idx], P[b][idx]).statistic for b in BASELINES}
        for a, b in CONTRASTS:
            fixed[f"{a}-{b}"].append(rf[a] - rf[b])
            rerank[f"{a}-{b}"].append(rr[a] - rr[b])

    frozen = json.loads((RES / "E15_results.json").read_text())["primary"]
    out = {
        "label": "POST-HOC (not pre-registered). Does not alter the frozen E15 decision.",
        "n_units": int(len(d)), "n_components": int(K), "B": B, "seed": SEED, "MEI": MEI,
        "rho_point_recomputed": rho,
        "rho_point_frozen": frozen["rho"],
        "check1_independent_reimplementation_fixed_ranks": summarise(fixed, point),
        "check2_conventional_reranked_bootstrap": summarise(rerank, point),
        "frozen_primary": {k: {"delta_rho": v["delta_rho"], "ci95": v["ci95"], "decision": v["decision"]}
                           for k, v in frozen["contrasts"].items()},
    }
    (RES / "posthoc_robustness.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

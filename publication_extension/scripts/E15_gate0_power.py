"""E15 Gate 0: operating characteristics of the pre-registered test, from real predictions and
simulated truth. NO QTL outcome is read.

What this is (and is not):
- It is a power / error-rate calculation for the pre-registered decision rule, given the actual
  eligible pairs, the actual predictions under F1, A1, G and U, and the actual gene-component structure
  used by the bootstrap.
- The "fetal-favourable" scenarios simulate truth that is, by construction, closer to the F1
  prediction. Power under them is the probability of detecting a fetal advantage IF one exists. It is
  NOT evidence that fetal predictions are correct, and Gate 0 passing says nothing about the direction
  of the real result.

Data-generating model (rank-normal scale of each prediction, x_B = Φ^-1(rank)):
    truth y = a·x_F1 + b·x_A1 + c·x_G + noise,   noise ~ N(0,1) unless stated
Scenarios:
  N0  pure null: a = b = c = 0                     (no predictor carries signal)
  N1  equal-signal null: a = b (> 0), c = 0        (F1 and A1 equally good; Δρ = 0 in expectation)
  F   fetal-favourable: (a, b) solved so that ρ(F1, y) = rho_level and Δρ = ρF1 − ρA1 = target
  A   adult-favourable: mirror of F (Δρ = −target)
  G   generic-favourable: truth from x_G; reports the F1−G and F1−A1 decisions
  M   mixture: a = 0.7, b = 0.3, noise calibrated to rho_level
  F-t3, F-het: fetal-favourable with t(3) noise, and with heteroscedastic noise (per-pair SD lognormal,
       σ_log = 0.5), to check robustness of the bootstrap
Effect grid: rho_level ∈ {0.10, 0.20, 0.30}, target Δρ ∈ {0.02, 0.05, 0.10}.

Decision rule simulated (= pre-registered rule):
  family = {F1−A1, F1−G, F1−U}; percentile CI from a paired component bootstrap; Holm at α = 0.05
  (two-sided). "Fetal improvement over B" is declared iff the Holm-adjusted CI for F1−B lies above 0.
Gate 0 passes iff power for F1−A1 ≥ 0.80 in scenario F at rho_level = 0.20, Δρ = 0.05, AND the
Type I error (false declaration of F1−A1 improvement) in N0 and N1 is ≤ 0.05 (simulation SE reported).

Input TSV columns: unit_id, component, pred_F1, pred_A1, pred_G, pred_U
Output: gate0_power.json next to the input.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_stats as S  # noqa: E402

SEED = 20261106
R_SIMS = 500
B_BOOT = 500
ALPHA = 0.05
CONTRASTS = ("A1", "G", "U")


def rank_normal(x):
    return norm.ppf((rankdata(x, method="average") - 0.5) / len(x))


def rho(y, x):
    return np.corrcoef(rankdata(y), rankdata(x))[0, 1]


def simulate(X, comps, coef, noise="normal", sigma=1.0, rng=None, r=R_SIMS, b=B_BOOT):
    n = len(X["F1"])
    sig = coef.get("F1", 0) * X["F1"] + coef.get("A1", 0) * X["A1"] + coef.get("G", 0) * X["G"]
    R = {B: S.ranks(X[B]) for B in X}
    Rd = S.ranks(X["F1"] - X["A1"])          # candidate alternative endpoint: difference correlation
    decl = {k: 0 for k in CONTRASTS}
    worse = {k: 0 for k in CONTRASTS}
    est = {k: [] for k in CONTRASTS}
    rhoF = []
    dif = {"est": [], "pos_005": 0, "pos_0167": 0, "neg_005": 0}
    for _ in range(r):
        if noise == "t3":
            e = rng.standard_t(3, n) / np.sqrt(3.0)
        elif noise == "het":
            e = rng.standard_normal(n) * np.exp(0.5 * rng.standard_normal(n))
        else:
            e = rng.standard_normal(n)
        ry = S.ranks(sig + sigma * e)
        rhoF.append(np.corrcoef(ry, R["F1"])[0, 1])
        W = S.component_weights(comps, b, rng)
        point, boot = S.contrasts(ry, R, others=CONTRASTS, W=W)
        better, worse_set, _, _ = S.holm(boot, point, ALPHA)
        for k in CONTRASTS:
            est[k].append(point[k])
            decl[k] += k in better
            worse[k] += k in worse_set
        pd_ = np.corrcoef(ry, Rd)[0, 1]
        bd = S.wcorr(W, ry, Rd)
        pv = min(1.0, 2 * min((bd <= 0).mean(), (bd >= 0).mean()))
        dif["est"].append(pd_)
        dif["pos_005"] += (pv <= 0.05 and pd_ > 0)
        dif["pos_0167"] += (pv <= 0.05 / 3 and pd_ > 0)
        dif["neg_005"] += (pv <= 0.05 and pd_ < 0)
    out = {"mean_rho_F1_truth": float(np.mean(rhoF)),
           "difference_endpoint": {"mean_rho_diff": float(np.mean(dif["est"])),
                                   "rate_declare_F1_direction_better_alpha0.05": dif["pos_005"] / r,
                                   "rate_declare_F1_direction_better_alpha0.0167": dif["pos_0167"] / r,
                                   "rate_declare_opposite_alpha0.05": dif["neg_005"] / r}}
    for k in CONTRASTS:
        rate = decl[k] / r
        out[f"F1-{k}"] = {"mean_delta_rho": float(np.mean(est[k])),
                          "rate_declare_F1_better": rate, "mc_se": float(np.sqrt(rate * (1 - rate) / r)),
                          "rate_declare_F1_worse": worse[k] / r}
    return out


def solve_coef(X, rho_level, delta):
    """a, b such that (Pearson on rank-normal scale) ρ(F1,y)=rho_level and ρ(A1,y)=rho_level−delta,
    noise variance 1. Pearson on normal scores ≈ Spearman; verified by simulation output."""
    r = np.corrcoef(X["F1"], X["A1"])[0, 1]

    def f(v):
        a, b = v
        sd = np.sqrt(a * a + b * b + 2 * a * b * r + 1)
        return [(a + b * r) / sd - rho_level, (a * r + b) / sd - (rho_level - delta)]

    sol, info, ok, _ = fsolve(f, [rho_level, 0.0], full_output=True)
    if ok != 1 or np.max(np.abs(f(sol))) > 1e-6:
        return None                         # target unattainable given the prediction correlation
    return sol


def calibrate(X, lvl, d, rng, reps=200, iters=6):
    """Coefficients whose simulated MEAN Spearman ρ(F1,y) and Δρ = ρF1 − ρA1 match (lvl, d).
    Starts from the Pearson/normal-scale solution and corrects the targets iteratively."""
    RF, RA = S.ranks(X["F1"]), S.ranks(X["A1"])
    lt, dt = lvl, d
    for _ in range(iters):
        c = solve_coef(X, lt, dt)
        if c is None:
            return None                     # unattainable
        mf, md = [], []
        for _ in range(reps):
            ry = S.ranks(c[0] * X["F1"] + c[1] * X["A1"] + rng.standard_normal(len(RF)))
            f = np.corrcoef(ry, RF)[0, 1]
            mf.append(f)
            md.append(f - np.corrcoef(ry, RA)[0, 1])
        mf, md = float(np.mean(mf)), float(np.mean(md))
        if abs(mf - lvl) < 0.003 and abs(md - d) < 0.0015:
            return c
        lt += lvl - mf
        dt += d - md
    return None                             # could not reach the target: treated as unattainable


def main(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    if "component" not in df and "component_P" in df:
        df = df.rename(columns={"component_P": "component"})
    df = df.dropna(subset=["pred_F1", "pred_A1", "pred_G", "pred_U"])
    X = {B: rank_normal(df[f"pred_{B}"].to_numpy()) for B in ("F1", "A1", "G", "U")}
    comps = pd.factorize(df.component)[0]
    rng = np.random.default_rng(SEED)
    rep = {"what_this_is": "Operating characteristics of the pre-registered test under simulated truth. "
                           "Fetal-favourable scenarios assume a fetal advantage by construction; passing "
                           "Gate 0 is not evidence that fetal predictions are correct.",
           "n_units": int(len(df)), "n_components": int(len(np.unique(comps))),
           "spearman_between_predictions": {f"{p}-{q}": float(rho(df[f"pred_{p}"], df[f"pred_{q}"]))
                                            for p, q in [("F1", "A1"), ("F1", "G"), ("F1", "U"), ("A1", "G")]},
           "settings": {"R_SIMS": R_SIMS, "B_BOOT": B_BOOT, "SEED": SEED, "ALPHA": ALPHA, "family": list(CONTRASTS)}}
    a = np.abs(df.pred_F1.to_numpy()); b = np.abs(df.pred_A1.to_numpy())
    with np.errstate(divide="ignore", invalid="ignore"):
        rep["share_units_F1_A1_differ_>2fold"] = float(np.mean(np.maximum(a, b) / np.minimum(a, b) > 2))
    r_fa = float(np.corrcoef(X["F1"], X["A1"])[0, 1])
    rep["pearson_rank_normal_F1_A1"] = r_fa
    rep["max_achievable_delta_rho"] = 1 - r_fa
    sc = {}

    def run(name, coef, **kw):
        sc[name] = simulate(X, comps, coef, rng=rng, **kw) if coef is not None else {"unattainable": True}

    def fa(c):
        return None if c is None else {"F1": c[0], "A1": c[1]}

    run("N0_pure_null", {})
    run("N1_equal_signal_null(rho=0.2)", fa(calibrate(X, 0.20, 0.0, rng)))
    for lvl in (0.10, 0.20, 0.30):
        for d in (0.02, 0.05, 0.10):
            if d >= 0.05 or lvl == 0.20:
                run(f"F_fetal(rho={lvl},d={d})", fa(calibrate(X, lvl, d, rng)))
    # standardized scenarios: Δρ as a fraction s of the maximum attainable gap (1 − r)
    for s in (0.25, 0.50, 0.75):
        d = round(s * (1 - r_fa), 4)
        run(f"F_fetal(rho=0.2,s={s},d={d})", fa(calibrate(X, 0.20, d, rng)))
    c = fa(calibrate(X, 0.20, 0.05, rng))
    run("F_t3_noise(rho=0.2,d=0.05)", c, noise="t3")
    run("F_hetero_noise(rho=0.2,d=0.05)", c, noise="het")
    d50 = round(0.5 * (1 - r_fa), 4)
    c50 = fa(calibrate(X, 0.20, d50, rng))
    run(f"F_t3_noise(rho=0.2,s=0.5,d={d50})", c50, noise="t3")
    run(f"F_hetero_noise(rho=0.2,s=0.5,d={d50})", c50, noise="het")
    sw = {"F1": X["A1"], "A1": X["F1"]}
    ca = calibrate(sw, 0.20, 0.05, rng)
    run("A_adult(rho=0.2,d=-0.05)", None if ca is None else {"A1": ca[0], "F1": ca[1]})
    ca50 = calibrate(sw, 0.20, d50, rng)
    run(f"A_adult(rho=0.2,s=-0.5,d=-{d50})", None if ca50 is None else {"A1": ca50[0], "F1": ca50[1]})
    s_g = 0.2 / np.sqrt(1 - 0.04)          # ρ(G, y) = 0.20 with unit noise
    run("G_generic(rho_G=0.2)", {"G": s_g})
    s_m = 0.2 / np.sqrt(1 - 0.04)
    run("M_mixture_0.7F_0.3A(rho~0.2)", {"F1": 0.7 * s_m, "A1": 0.3 * s_m})
    rep["scenarios"] = sc
    main_sc = sc["F_fetal(rho=0.2,d=0.05)"]
    t1 = max(sc["N0_pure_null"]["F1-A1"]["rate_declare_F1_better"],
             sc["N1_equal_signal_null(rho=0.2)"]["F1-A1"]["rate_declare_F1_better"])
    if main_sc.get("unattainable"):
        rep["gate0"] = {"power_F1-A1_at_rho0.2_d0.05": None, "max_typeI_F1-A1_N0_N1": t1,
                        "verdict": "FAIL",
                        "reason": f"Δρ = 0.05 is unattainable: max achievable Δρ = 1 − r = {1 - r_fa:.3f}"}
    else:
        p_main = main_sc["F1-A1"]["rate_declare_F1_better"]
        rep["gate0"] = {"power_F1-A1_at_rho0.2_d0.05": p_main, "max_typeI_F1-A1_N0_N1": t1,
                        "verdict": "PASS" if (p_main >= 0.80 and t1 <= 0.05) else "FAIL"}
    rep["standardized_effect_scenarios_delta_rho_test"] = {
        k2: (None if v.get("unattainable") else
             {"mean_delta_rho": v["F1-A1"]["mean_delta_rho"], "power": v["F1-A1"]["rate_declare_F1_better"]})
        for k2, v in sc.items() if ",s=" in k2}
    t1d = max(sc["N0_pure_null"]["difference_endpoint"]["rate_declare_F1_direction_better_alpha0.05"],
              sc["N1_equal_signal_null(rho=0.2)"]["difference_endpoint"]["rate_declare_F1_direction_better_alpha0.05"])
    rep["difference_correlation_endpoint_diagnostic"] = {
        "typeI_max_N0_N1_alpha0.05": t1d,
        "valid": bool(t1d <= 0.05),
        "note": "rank correlation of truth with rank(x_F1 - x_A1); reported only to document whether it is a valid test"}
    return rep


if __name__ == "__main__":
    p = Path(sys.argv[1])
    rep = main(p)
    (p.parent / "gate0_power.json").write_text(json.dumps(rep, indent=1, default=float))
    print(json.dumps(rep["gate0"], indent=1))

"""E15 shared statistics: used identically by Gate 0 (simulated truth) and by the frozen analysis.

Estimator: Spearman ρ between truth and each baseline's prediction = Pearson correlation of
full-sample average ranks.
Uncertainty: paired component bootstrap. Each replicate draws n_components components with
replacement; a unit's weight is its component's draw count; ρ in the replicate is the weighted
Pearson correlation of the FULL-SAMPLE ranks. (Ranks are not recomputed within replicates; this keeps
the computation vectorised and is used identically in Gate 0 and in the analysis.)
Bootstrap p-value: p = min(1, 2·min(P*(Δ ≤ 0), P*(Δ ≥ 0))). Holm across the contrast family.
"""
import numpy as np
from scipy.stats import rankdata


def ranks(x):
    return rankdata(x, method="average")


def component_weights(comps, B, rng):
    """comps: integer component code per unit (0..K-1). Returns B x n weight matrix."""
    K = comps.max() + 1
    counts = rng.multinomial(K, np.full(K, 1.0 / K), size=B)      # B x K
    return counts[:, comps].astype(float)


def wcorr(W, a, b):
    s = W.sum(1)
    ma, mb = W @ a / s, W @ b / s
    cov = W @ (a * b) / s - ma * mb
    va = W @ (a * a) / s - ma ** 2
    vb = W @ (b * b) / s - mb ** 2
    return cov / np.sqrt(va * vb)


def contrasts(ry, R, ref="F1", others=("A1", "G", "U"), W=None):
    """ry: truth ranks; R: dict baseline -> prediction ranks. Point Δρ (ref − other), and bootstrap
    arrays if W is given."""
    def c(a, b):
        return np.corrcoef(a, b)[0, 1]
    point = {o: c(ry, R[ref]) - c(ry, R[o]) for o in others}
    if W is None:
        return point, None
    rref = wcorr(W, ry, R[ref])
    boot = {o: rref - wcorr(W, ry, R[o]) for o in others}
    return point, boot


def holm(boot, point, alpha=0.05):
    p = {k: min(1.0, 2 * min((v <= 0).mean(), (v >= 0).mean())) for k, v in boot.items()}
    order = sorted(p, key=p.get)
    m = len(order)
    adj, running = {}, 0.0
    for i, k in enumerate(order):
        running = max(running, min(1.0, (m - i) * p[k]))
        adj[k] = running
    better = {k for k in p if adj[k] <= alpha and point[k] > 0}
    worse = {k for k in p if adj[k] <= alpha and point[k] < 0}
    return better, worse, p, adj


def holm_ci(boot, adj_alpha):
    """Percentile CI at the Holm step's level for each contrast: {k: (lo, hi)}."""
    return {k: (float(np.quantile(v, adj_alpha[k] / 2)), float(np.quantile(v, 1 - adj_alpha[k] / 2)))
            for k, v in boot.items()}

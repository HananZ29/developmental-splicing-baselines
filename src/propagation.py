"""Propagate SpliceAI site-level changes into isoform proportions under a chosen baseline.

Model (a hypothesis, not a biological law):
  Each annotated splice site s has a usage efficiency proxied by its SpliceAI probability. The variant
  rescales that efficiency by r_s = (p_alt(s) + eps) / (p_ref(s) + eps). Assuming a transcript's
  production rate is proportional to its baseline share times the product of its sites'
  efficiencies (sites act independently and multiplicatively, i.e. additively on the logit scale),
      f_t   = prod_{s in sites(t), |s - v| <= W} r_s
      pi'_t = pi_t * f_t / sum_u pi_u * f_u
  Predicted isoform change: d_pi_t = pi'_t - pi_t.
  Predicted intron usage within its competition set C(i) (annotated introns sharing a donor or
  acceptor with i, as in leafcutter clustering):
      u_i = sum_{t contains i} pi_t / sum_{j in C(i)} sum_{t contains j} pi_t ;  d_u_i = u'_i - u_i
  Not modelled: cryptic sites absent from the annotation, unannotated isoforms, NMD, and
  transcription-level (eQTL) effects. The baseline pi is the only thing that differs between
  conditions.
"""
import numpy as np

EPS = 0.01     # regularizes ratios of low-probability sites; sensitivity analysis varies it
ALPHA = 1e-3   # pseudocount added identically to every baseline before renormalization


def smooth_baseline(raw, transcripts, alpha=ALPHA):
    """raw: {enst: proportion or abundance}; returns array over `transcripts` (missing -> 0)."""
    p = np.array([max(raw.get(t, 0.0), 0.0) for t in transcripts], dtype=float)
    p = p / p.sum() if p.sum() > 0 else np.full(len(p), 1.0 / len(p))
    p = p + alpha
    return p / p.sum()


def transcript_factors(transcripts, var_pos, tracks, window, eps=EPS, max_dist=None):
    """transcripts: list of splice_graph.Transcript; tracks: dict acc_ref/acc_alt/don_ref/don_alt arrays
    covering var_pos-window..var_pos+window. Sites farther than max_dist (default: window) get r_s = 1."""
    max_dist = window if max_dist is None else min(max_dist, window)
    logf = np.zeros(len(transcripts))
    for k, t in enumerate(transcripts):
        for pos, kind in t.sites():
            o = pos - var_pos
            if abs(o) > max_dist:
                continue
            key = "acc" if kind == "acceptor" else "don"
            pr, pa = tracks[f"{key}_ref"][window + o], tracks[f"{key}_alt"][window + o]
            logf[k] += np.log((pa + eps) / (pr + eps))
    return np.exp(logf)


def propagate(pi, f):
    new = pi * f
    return new / new.sum()


def intron_membership(transcripts, mode="component"):
    """Intron x transcript membership and each intron's competition set.

    mode="component" (v2, primary): LeafCutter-style cluster = connected component of introns linked
    by a shared start or end coordinate (transitive), restricted to annotated introns.
    mode="direct" (v1): only introns sharing a start or end with the intron itself.
    """
    introns = sorted({i for t in transcripts for i in t.introns})
    index = {i: k for k, i in enumerate(introns)}
    m = np.zeros((len(introns), len(transcripts)))
    for c, t in enumerate(transcripts):
        for i in t.introns:
            m[index[i], c] = 1.0
    by_start, by_end = {}, {}
    for i in introns:
        by_start.setdefault(i[0], []).append(index[i])
        by_end.setdefault(i[1], []).append(index[i])
    if mode == "direct":
        comp = [sorted(set(by_start[i[0]]) | set(by_end[i[1]])) for i in introns]
        return introns, index, m, comp
    parent = list(range(len(introns)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for grp in list(by_start.values()) + list(by_end.values()):
        for j in grp[1:]:
            parent[find(j)] = find(grp[0])
    members = {}
    for k in range(len(introns)):
        members.setdefault(find(k), []).append(k)
    comp = [sorted(members[find(k)]) for k in range(len(introns))]
    return introns, index, m, comp


def intron_usage(pi, m, comp):
    inc = m @ pi
    return np.array([inc[k] / inc[c].sum() if inc[c].sum() > 0 else np.nan for k, c in enumerate(comp)])

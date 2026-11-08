"""E15 outcome-blind baseline construction and pre-freeze diagnostics.

Reads ONLY junction read counts (BrainSpan MRF extraction, GTEx per-tissue counts), the Wen intron
universe (all.introns.tested.tsv: intron IDs and genes, no statistics), sample metadata and GENCODE
v31lift37. No QTL statistic (z, PIP, slope, p-value) is opened.

Units and groups (fixed by DATA/decision docs):
- BrainSpan unit = donor (reads summed over that donor's neocortex samples).
  F1     = 12 prenatal donors <= 26 pcw not in the Wen cohort (decision D1)
  F1_all = F1 + HSB155 (35 pcw)   [sensitivity]
  A1     = 8 donors >= 18 y
- GTEx unit = tissue (SMTSD; reads summed over samples).  G = all tissues; G_nonbrain [sensitivity]
- U (uniform) needs no data.

Usage per group (rule R3): for each unit with >= MIN_UNIT_READS annotated reads in a cluster, the
within-cluster proportions over annotated introns; the group usage is the unweighted mean over those
units.

Eligibility (rules R1/R2 + PREFREEZE point 3), applied identically to every baseline:
  E1 cluster has >= 2 annotated introns (v31lift37) and the pattern-by-intron incidence matrix has full
     column rank (pattern weights identifiable from junction counts)
  E2 in each of F1, A1, G: pooled annotated reads >= MIN_POOLED and >= MIN_UNITS contributing units
  E3 in each of F1, A1, G: share of the cluster's reads on annotated introns >= MIN_ANN_SHARE

Diagnostics written to publication_extension/results/E15_devstage/baseline_diagnostics.json:
  counts per baseline, eligibility waterfall, F1 stability (leave-one-donor-out and donor bootstrap,
  intron usage and pattern weights), and exclusion bias by expression, annotation and F1-vs-A1
  difference.
Usage tables: data/processed/E15_baselines/usage_<group>.tsv.gz, cluster_eligibility.tsv.gz
"""
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
JDIR = ROOT / "data/processed/E15_junctions"
MAN = ROOT / "data/external/brainspan_mrf/neocortex/manifest.tsv"
META = ROOT / "data/external/brainspan_pec/meta/mRNA-seq_Sample_metadata.xlsx"
WEN = ROOT / "data/external/wen_synapse/all.introns.tested.tsv"
GTF = ROOT / "data/external/gencode/gencode.v31lift37.annotation.gtf.gz"
BDIR = ROOT / "data/processed/E15_baselines"
RES = ROOT / "publication_extension/results/E15_devstage"

WEN_OVERLAP = {"HSB107", "HSB148", "HSB150", "HSB154", "HSB194"}
ADULT = {"HSB105", "HSB123", "HSB126", "HSB127", "HSB130", "HSB135", "HSB136", "HSB145"}
STANDARD_REGIONS = {"A1C", "DFC", "IPC", "ITC", "M1C", "MFC", "OFC", "S1C", "STC", "V1C", "VFC"}
MIN_UNIT_READS = 5
MIN_POOLED = 20
MIN_UNITS = 3
MIN_ANN_SHARE = 0.95
N_BOOT = 200
SEED = 20261106


# ---------------------------------------------------------------- inputs
def wen_universe():
    rows = []
    with open(WEN) as f:
        for line in f:
            iid, gene = line.rstrip("\n").split("\t")[:2]
            c, a, b, clu = iid.split(":")[:4]
            c = "chr" + c.replace("chr", "")
            rows.append((c, int(a) + 1, int(b) - 1, c + ":" + clu, gene))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "cluster", "gene"])
    df = df.drop_duplicates(["chrom", "start", "end", "cluster"])
    return df.sort_values(["cluster", "start", "end"]).reset_index(drop=True)


def transcript_introns():
    ex = defaultdict(list)
    chrom = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "exon":
                continue
            t = p[8].split('transcript_id "')[1].split('"')[0]
            ex[t].append((int(p[3]), int(p[4])))
            chrom[t] = p[0]
    out = {}
    for t, e in ex.items():
        if len(e) < 2:
            continue
        e.sort()
        out[t] = frozenset((chrom[t], e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1))
    return out


def donor_ages():
    m = pd.read_excel(META, header=3).dropna(subset=["Braincode"])
    return dict(zip(m["Braincode"], m["Age"].astype(str)))


def load_brainspan(W, regions=None):
    """Return {donor: counts array over W rows}, list of entries used, per-donor sample counts."""
    key = {(c, s, e): i for i, (c, s, e) in enumerate(zip(W.chrom, W.start, W.end))}
    man = pd.read_csv(MAN, sep="\t")
    by = {}
    nsamp = defaultdict(int)
    for entry, region in zip(man.iloc[:, 1], man.iloc[:, 2]):
        if regions is not None and region not in regions:
            continue
        p = JDIR / f"{entry}.junctions.tsv.gz"
        st = JDIR / f"{entry}.stats.json"
        if not p.exists() or not st.exists() or not json.loads(st.read_text()).get("crc_ok"):
            continue
        donor = "HSB" + entry.split("_")[1]
        v = by.setdefault(donor, np.zeros(len(W), dtype=np.int64))
        with gzip.open(p, "rt") as f:
            f.readline()
            for line in f:
                c, s, e, r = line.split("\t")
                i = key.get((c, int(s), int(e)))
                if i is not None:
                    v[i] += int(r)
        nsamp[donor] += 1
    return by, dict(nsamp)


def load_gtex(W):
    g = pd.read_csv(BDIR / "gtex_tissue_counts.tsv.gz", sep="\t")
    tis = pd.read_csv(BDIR / "gtex_tissue_samples.tsv", sep="\t")
    idx = W.reset_index().merge(g, on=["chrom", "start", "end"], how="left")
    by = {}
    for t in tis.tissue:
        by[t] = idx[t].fillna(0).to_numpy(np.int64)
    present = idx[tis.tissue.iloc[0]].notna().to_numpy()
    return by, tis, present


# ---------------------------------------------------------------- core computations
class Clusters:
    def __init__(self, W, ann):
        self.cl = W.cluster.to_numpy()
        self.starts = np.flatnonzero(np.r_[True, self.cl[1:] != self.cl[:-1]])
        self.size = np.diff(np.r_[self.starts, len(self.cl)])
        self.ids = self.cl[self.starts]
        self.ann = ann.astype(np.int64)

    def tot(self, X):
        return np.add.reduceat(X, self.starts, axis=-1)

    def expand(self, Y):
        return np.repeat(Y, self.size, axis=-1)


def unit_props(C, K):
    """C: units x introns counts. Proportions over annotated introns; NaN if annotated total < MIN."""
    A = C * K.ann
    at = K.tot(A)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = A / K.expand(at)
    P[K.expand(at) < MIN_UNIT_READS] = np.nan
    P[:, K.ann == 0] = np.nan
    return P, at


def group_usage(C, K):
    P, at = unit_props(C, K)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = np.nanmean(P, axis=0)
    n_units = (at >= MIN_UNIT_READS).sum(0)
    pooled_ann = at.sum(0)
    pooled_all = K.tot(C).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        share = pooled_ann / pooled_all
    return u, n_units, pooled_ann, pooled_all, share


def patterns_for_clusters(W, K, tx):
    """For each cluster: list of patterns (tuples of local annotated-intron indices)."""
    keys = list(zip(W.chrom.tolist(), W.start.tolist(), W.end.tolist()))
    pos = {k: i for i, k in enumerate(keys)}
    intron_tx = defaultdict(set)
    for t, ins in tx.items():
        for k in ins:
            i = pos.get(k)
            if i is not None:
                intron_tx[i].add(t)
    pats = []
    for ci, (st, n) in enumerate(zip(K.starts, K.size)):
        idx = list(range(st, st + n))
        ts = set().union(*(intron_tx.get(i, set()) for i in idx))
        P = set()
        for t in ts:
            p = tuple(j - st for j in idx if keys[j] in tx[t])
            if p:
                P.add(p)
        pats.append(sorted(P))
    return pats


def full_rank(pats, n):
    if len(pats) < 2:
        return False
    A = np.zeros((n, len(pats)))
    for j, p in enumerate(pats):
        A[list(p), j] = 1
    return np.linalg.matrix_rank(A) == len(pats)


def em_weights(u_local, pats, n, iters=500, tol=1e-10):
    """Molecule-level pattern weights w from intron usage u (pseudo-counts), EM for the model
    P(read on intron i) ∝ sum_p w_p 1[i in p]."""
    A = np.zeros((n, len(pats)))
    for j, p in enumerate(pats):
        A[list(p), j] = 1
    size = A.sum(0)
    w = np.full(len(pats), 1.0 / len(pats))
    x = np.nan_to_num(u_local)
    for _ in range(iters):
        q = A @ w
        with np.errstate(invalid="ignore", divide="ignore"):
            R = np.where(q[:, None] > 0, A * w[None, :] / q[:, None], 0.0)
        nw = (x[:, None] * R).sum(0) / size
        nw = nw / nw.sum() if nw.sum() > 0 else w
        if np.abs(nw - w).max() < tol:
            w = nw
            break
        w = nw
    return w


def q(a, ps=(5, 25, 50, 75, 95)):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return {f"p{p}": (float(np.percentile(a, p)) if len(a) else None) for p in ps} | {"n": int(len(a))}


# ---------------------------------------------------------------- main
def main():
    RES.mkdir(parents=True, exist_ok=True)
    BDIR.mkdir(parents=True, exist_ok=True)
    W = wen_universe()
    tx = transcript_introns()
    ann_set = set().union(*tx.values())
    ann = np.array([(c, s, e) in ann_set for c, s, e in zip(W.chrom, W.start, W.end)])
    K = Clusters(W, ann)
    ages = donor_ages()

    bs, nsamp = load_brainspan(W)
    pcw = {d: int(a.replace("PCW", "")) for d, a in ages.items() if a.endswith("PCW")}
    F1 = sorted(d for d in bs if d in pcw and pcw[d] <= 26 and d not in WEN_OVERLAP)
    F1_all = sorted(d for d in bs if d in pcw and d not in WEN_OVERLAP)
    A1 = sorted(d for d in bs if d in ADULT)
    gt, tis, gtex_present = load_gtex(W)
    G = list(tis.tissue)
    Gnb = list(tis.tissue[tis.is_brain == 0])

    groups = {"F1": (bs, F1), "F1_all": (bs, F1_all), "A1": (bs, A1), "G": (gt, G), "G_nonbrain": (gt, Gnb)}
    mats, res = {}, {}
    for g, (src, units) in groups.items():
        C = np.vstack([src[u] for u in units]) if units else np.zeros((0, len(W)), np.int64)
        mats[g] = C
        res[g] = group_usage(C, K)

    pats = patterns_for_clusters(W, K, tx)
    n_ann_cl = K.tot(K.ann)
    rank_ok = np.array([full_rank(p, n) for p, n in zip(pats, K.size)])

    def cov_ok(g):
        _, nu, pa, _, sh = res[g]
        return (pa >= MIN_POOLED) & (nu >= MIN_UNITS), (sh >= MIN_ANN_SHARE)

    e1 = (n_ann_cl >= 2) & rank_ok
    cov = {g: cov_ok(g) for g in groups}
    e2 = cov["F1"][0] & cov["A1"][0] & cov["G"][0]
    e3 = cov["F1"][1] & cov["A1"][1] & cov["G"][1]
    # E4: every annotated intron of the cluster exists in the GTEx junction set (GTEx v8 lists only
    # GENCODE v26 junctions; an intron absent there would be a structural zero in G, not a measurement)
    e4 = K.tot((K.ann * ~gtex_present).astype(np.int64)) == 0
    elig = e1 & e2 & e3 & e4

    # ---- per-baseline report
    per = {}
    for g, (src, units) in groups.items():
        u, nu, pa, pall, sh = res[g]
        per[g] = {
            "n_units": len(units),
            "units": units if g.startswith(("F1", "A1")) else f"{len(units)} GTEx tissues",
            "samples_per_unit": {d: nsamp.get(d) for d in units} if g.startswith(("F1", "A1")) else None,
            "reads_on_wen_introns_total": int(mats[g].sum()),
            "reads_on_wen_introns_per_unit": q(mats[g].sum(1)),
            "wen_introns_with_>=1_read": int((mats[g].sum(0) > 0).sum()),
            "clusters_coverage_ok(E2_single)": int(cov[g][0].sum()),
            "clusters_coverage_and_annshare_ok": int((cov[g][0] & cov[g][1]).sum()),
            "clusters_coverage_ok_and_E1": int((cov[g][0] & cov[g][1] & e1).sum()),
            "pooled_annotated_reads_per_eligible_cluster": q(pa[elig]),
        }
    waterfall = {
        "wen_clusters": int(len(K.ids)),
        "E1_>=2_annotated_introns": int((n_ann_cl >= 2).sum()),
        "E1_and_rank_ok": int(e1.sum()),
        "E1_and_E2_coverage_F1_A1_G": int((e1 & e2).sum()),
        "E1_and_E4_all_annotated_introns_in_GTEx": int((e1 & e4).sum()),
        "E1_E4_and_E2_coverage_F1_A1_G": int((e1 & e4 & e2).sum()),
        "E1_E2_E3_E4_eligible": int(elig.sum()),
        "eligible_introns": int(K.expand(elig).sum()),
        "eligible_annotated_introns": int((K.expand(elig) & ann).sum()),
    }

    # ---- F1 stability (on eligible clusters; also on F1-covered clusters)
    rng = np.random.default_rng(SEED)
    u_full = res["F1"][0]
    elig_i = K.expand(elig) & ann
    lodo = {}
    lodo_cluster_max = []
    for k, d in enumerate(F1):
        keep = [j for j in range(len(F1)) if j != k]
        u_k = group_usage(mats["F1"][keep], K)[0]
        dv = np.abs(u_k - u_full)
        cm = np.maximum.reduceat(np.nan_to_num(dv, nan=0.0), K.starts)[elig]
        lodo_cluster_max.append(cm)
        lodo[d] = {"median_cluster_max_abs_du": float(np.median(cm)) if len(cm) else None,
                   "frac_clusters_max_abs_du_>0.05": float((cm > 0.05).mean()) if len(cm) else None,
                   "frac_clusters_max_abs_du_>0.10": float((cm > 0.10).mean()) if len(cm) else None}
    lodo_worst = np.max(np.vstack(lodo_cluster_max), axis=0) if lodo_cluster_max else np.array([])

    def boot_sd(C):
        B = []
        n = C.shape[0]
        for _ in range(N_BOOT):
            B.append(group_usage(C[rng.integers(0, n, n)], K)[0])
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.nanstd(np.vstack(B), axis=0)

    sd_F1 = boot_sd(mats["F1"])
    sd_A1 = boot_sd(mats["A1"])
    diff = np.abs(res["F1"][0] - res["A1"][0])
    z = diff / np.sqrt(sd_F1 ** 2 + sd_A1 ** 2)
    cl_diff = np.maximum.reduceat(np.nan_to_num(diff, nan=0.0), K.starts)
    cl_dist = np.maximum.reduceat(np.nan_to_num(np.where(diff > 0.05, z, 0), nan=0.0), K.starts) > 2

    # pattern-weight stability (LODO) on eligible clusters
    elig_idx = np.flatnonzero(elig)
    w_lodo_max = []
    lodo_u = [group_usage(mats["F1"][[j for j in range(len(F1)) if j != k]], K)[0] for k in range(len(F1))]
    for ci in elig_idx:
        st, n = K.starts[ci], K.size[ci]
        w0 = em_weights(u_full[st:st + n], pats[ci], n)
        m = 0.0
        for uk in lodo_u:
            m = max(m, float(np.abs(em_weights(uk[st:st + n], pats[ci], n) - w0).max()))
        w_lodo_max.append(m)

    stability = {
        "method": f"leave-one-donor-out over {len(F1)} F1 donors; donor bootstrap B={N_BOOT}, seed {SEED}",
        "lodo_per_donor": lodo,
        "lodo_worst_case_cluster_max_abs_du_eligible": q(lodo_worst),
        "frac_eligible_clusters_worst_lodo_max_abs_du_>0.05": float((lodo_worst > 0.05).mean()) if len(lodo_worst) else None,
        "bootstrap_sd_usage_F1_eligible_annotated_introns": q(sd_F1[elig_i]),
        "bootstrap_sd_usage_A1_eligible_annotated_introns": q(sd_A1[elig_i]),
        "lodo_worst_case_max_abs_dw_pattern_weights_eligible": q(w_lodo_max),
        "F1_vs_A1_cluster_max_abs_usage_difference_eligible": q(cl_diff[elig]),
        "frac_eligible_clusters_F1_A1_distinguishable(|du|>0.05_and_>2_bootstrap_SE)": float(cl_dist[elig].mean()) if elig.any() else None,
    }

    # ---- exclusion bias (item 3)
    f1_ok = cov["F1"][0] & cov["F1"][1] & e1   # clusters usable for the fetal baseline on its own
    pa_F1 = res["F1"][2]
    bias = {"universe": "clusters passing E1 and F1's own coverage/annotation-share criteria",
            "n": int(f1_ok.sum()), "n_eligible": int((f1_ok & elig).sum()),
            "frac_lost_E4_annotated_intron_missing_from_GTEx": float((~e4)[f1_ok].mean()) if f1_ok.any() else None,
            "frac_lost_A1_coverage_or_share": float((~(cov["A1"][0] & cov["A1"][1]))[f1_ok].mean()) if f1_ok.any() else None,
            "frac_lost_G_coverage": float((~cov["G"][0])[f1_ok].mean()) if f1_ok.any() else None}
    if f1_ok.any():
        qs = np.quantile(pa_F1[f1_ok], [0.25, 0.5, 0.75])
        strata = np.digitize(pa_F1, qs)
        bias["by_F1_expression_quartile"] = {
            f"Q{s + 1}": {"n": int((f1_ok & (strata == s)).sum()),
                          "frac_retained": float(elig[f1_ok & (strata == s)].mean()) if (f1_ok & (strata == s)).any() else None,
                          "frac_lost_A1_coverage": float((~cov["A1"][0])[f1_ok & (strata == s)].mean()) if (f1_ok & (strata == s)).any() else None,
                          "frac_lost_G_coverage": float((~cov["G"][0])[f1_ok & (strata == s)].mean()) if (f1_ok & (strata == s)).any() else None}
            for s in range(4)}
        bias["F1_pooled_reads_retained_vs_excluded"] = {"retained": q(pa_F1[f1_ok & elig]), "excluded": q(pa_F1[f1_ok & ~elig])}
        both = f1_ok & cov["A1"][0] & cov["A1"][1]
        bias["fetal_only_expressed(F1_ok_but_A1_coverage_fails)"] = {
            "n": int((f1_ok & ~cov["A1"][0]).sum()), "frac_of_F1_ok": float((f1_ok & ~cov["A1"][0]).mean() / f1_ok.mean())}
        if both.any():
            dq = [0.05, 0.10, 0.20]
            bias["by_F1_vs_A1_difference(F1_and_A1_covered)"] = {
                f"max|du|{lo}-{hi}": {"n": int((both & (cl_diff >= lo) & (cl_diff < hi)).sum()),
                                      "frac_retained(G_requirement)": float(elig[both & (cl_diff >= lo) & (cl_diff < hi)].mean()) if (both & (cl_diff >= lo) & (cl_diff < hi)).any() else None}
                for lo, hi in zip([0] + dq, dq + [1.01])}
    # annotation status of all clusters with F1 coverage regardless of E1/E3
    f1_cov_any = cov["F1"][0]
    bias["F1_covered_clusters_lost_to_annotation_rules(E1_or_E3)"] = {
        "n_F1_covered": int(f1_cov_any.sum()),
        "frac_failing_E1": float((~e1)[f1_cov_any].mean()) if f1_cov_any.any() else None,
        "frac_failing_F1_annotated_share": float((~cov["F1"][1])[f1_cov_any].mean()) if f1_cov_any.any() else None}

    out = {"parameters": {"MIN_UNIT_READS": MIN_UNIT_READS, "MIN_POOLED": MIN_POOLED, "MIN_UNITS": MIN_UNITS,
                          "MIN_ANN_SHARE": MIN_ANN_SHARE, "annotation": GTF.name},
           "brainspan_samples_used": int(sum(nsamp.values())),
           "per_baseline": per, "eligibility_waterfall": waterfall, "F1_stability": stability,
           "exclusion_bias": bias}
    tag = sys.argv[1] if len(sys.argv) > 1 else ""
    (RES / f"baseline_diagnostics{tag}.json").write_text(json.dumps(out, indent=1))

    # usage tables + eligibility
    for g in groups:
        u, nu, pa, pall, sh = res[g]
        df = W[["chrom", "start", "end", "cluster"]].copy()
        df["usage"] = u
        df["cluster_pooled_annotated_reads"] = K.expand(pa)
        df["cluster_contributing_units"] = K.expand(nu)
        df.to_csv(BDIR / f"usage_{g}{tag}.tsv.gz", sep="\t", index=False, compression="gzip")
    ce = pd.DataFrame({"cluster": K.ids, "n_introns": K.size, "n_annotated": n_ann_cl, "n_patterns": [len(p) for p in pats],
                       "E1": e1, "E2": e2, "E3": e3, "E4": e4, "eligible": elig})
    ce.to_csv(BDIR / f"cluster_eligibility{tag}.tsv.gz", sep="\t", index=False, compression="gzip")
    print(json.dumps({k: out[k] for k in ("brainspan_samples_used", "eligibility_waterfall")}, indent=1))


if __name__ == "__main__":
    main()

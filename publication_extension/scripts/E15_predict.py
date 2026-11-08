"""E15 units and predictions under every baseline (outcome-blind; PREREGISTRATION_E15 sections 3-5).

Reads: junction counts (BrainSpan, GTEx), Wen intron IDs, Wen fine-mapping columns 0-7 only
(phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index), Wen cluster->gene map, SpliceAI
scores/tracks, GENCODE v31lift37, and cluster eligibility from E15_baseline_diagnostics.py.
Never reads pip, z, posterior_*, cs_log10bf, or any conditional/GTEx sQTL statistic.

Output: publication_extension/results/E15_devstage/predictions.tsv with one row per unit:
  unit_id, phenotype_id, cs_id, cluster, intron (hg19 1-based), rule-P variant + predictions,
  rule-D variant + predictions, component_P, component_D, signal_group, DS_max_P,
  pred_<baseline>[_D] for F1, A1, G, U and sensitivity baselines.
"""
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_baseline_diagnostics as D  # noqa: E402

RES = D.RES
SUSIE = D.ROOT / "data/external/wen_synapse/mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz"
CL2GENE = D.ROOT / "data/external/wen_synapse/leafcutter_clusters_to_genes.txt"
SAI = RES / "spliceai"
EPS = 0.01
EPS_SENS = (0.001, 0.1)
WIN_TRACK = 500
N_SEEDS = 20
SEED0 = 20261106
LD_PROXY_BP = 100_000
SFRAC_MIN = 0.10


# ------------------------------------------------------------------ helpers
class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def pattern_matrix(pats, n):
    A = np.zeros((n, len(pats)))
    for j, p in enumerate(pats):
        A[list(p), j] = 1
    return A


def fit_weights(u, A, iters=2000, tol=1e-12):
    """MLE pattern weights (EM, PREREG section 3), initialised at the NNLS solution."""
    if u is None or not np.isfinite(u).any():
        return None
    x = np.nan_to_num(u)
    size = A.sum(0)
    w0, _ = nnls(A, x)
    w = w0 * size + 1e-6
    w = w / w.sum()
    for _ in range(iters):
        q = A @ w
        R = np.divide(A * w[None, :], q[:, None], out=np.zeros_like(A), where=q[:, None] > 0)
        nw = (x[:, None] * R).sum(0) / size
        nw = nw / nw.sum()
        if np.abs(nw - w).max() < tol:
            return nw
        w = nw
    return w


def usage_from_weights(w, A):
    return (A @ w) / (A.sum(0) @ w)


def unit_mean_usage(C):
    """C: units x local annotated introns (counts). Mean over units with >= MIN_UNIT_READS."""
    tot = C.sum(1)
    ok = tot >= D.MIN_UNIT_READS
    if ok.sum() == 0:
        return None
    return (C[ok] / tot[ok, None]).mean(0)


# ------------------------------------------------------------------ main
def main():
    W = D.wen_universe()
    tx = D.transcript_introns()
    ann_set = set().union(*tx.values())
    keys = list(zip(W.chrom.tolist(), W.start.tolist(), W.end.tolist()))
    ann = np.array([k in ann_set for k in keys])
    K = D.Clusters(W, ann)
    pats_all = D.patterns_for_clusters(W, K, tx)
    elig_tab = pd.read_csv(D.BDIR / "cluster_eligibility.tsv.gz", sep="\t")
    elig = dict(zip(elig_tab.cluster, elig_tab.eligible))
    cl_index = {c: i for i, c in enumerate(K.ids)}

    # strand per cluster (from candidates file, built from annotation)
    cand = pd.read_csv(RES / "candidates.tsv", sep="\t")
    cand = cand[cand.cluster.map(lambda c: bool(elig.get(c, False)))].copy()
    # SpliceAI scores/tracks: Wen candidates (spliceai/) + GTEx crossover candidates (spliceai_gtex/)
    scs, trs = [], []
    for d in (SAI, RES / "spliceai_gtex"):
        if (d / "scores.tsv").exists():
            scs.append(pd.read_csv(d / "scores.tsv", sep="\t"))
            t = np.load(d / "tracks.npz")
            trs.append({k: t[k] for k in ("acc_ref", "acc_alt", "don_ref", "don_alt")})
    sc = pd.concat(scs, ignore_index=True)
    tr = {k: np.concatenate([t[k] for t in trs]) for k in trs[0]}
    sc["row"] = np.arange(len(sc))
    sc = sc.drop_duplicates(["chrom", "pos", "ref", "alt", "strand"])

    def attach(cdf):
        cdf = cdf.merge(sc[["chrom", "pos", "ref", "alt", "strand", "DS_max", "row"]],
                        on=["chrom", "pos", "ref", "alt", "strand"], how="left")
        n_missing = int(cdf.row.isna().sum())
        cdf = cdf.dropna(subset=["row"]).copy()
        cdf["row"] = cdf.row.astype(int)
        return cdf, n_missing

    def select(cdf, key, prefix):
        cols = ["variant_id", "pos", "ref", "alt", "row", "DS_max", "dist_to_site"] + \
               (["allele_flipped"] if "allele_flipped" in cdf else [])
        cdf = cdf.sort_values(key + ["pos", "alt"])
        selP = cdf.sort_values(["DS_max", "pos", "alt"], ascending=[False, True, True]).groupby(key, sort=False).head(1)
        selD = cdf.sort_values(["dist_to_site", "pos", "alt"]).groupby(key, sort=False).head(1)
        u = selP.rename(columns={c: c + "_P" for c in cols})
        u = u.merge(selD[key + cols].rename(columns={c: c + "_D" for c in cols}), on=key)
        u = u.reset_index(drop=True)
        u["unit_id"] = [f"{prefix}{i:05d}" for i in range(len(u))]
        return u

    cand, miss = attach(cand)
    stats = {"candidate_rows_in_eligible_clusters": int(len(cand) + miss), "candidate_rows_without_spliceai": miss}
    units = select(cand, ["phenotype_id", "cs_id"], "U")
    stats["units"] = int(len(units))
    stats["units_same_variant_P_and_D"] = int((units.variant_id_P == units.variant_id_D).sum())
    gunits = None
    if (RES / "gtex_candidates.tsv").exists():
        gc = pd.read_csv(RES / "gtex_candidates.tsv", sep="\t")
        gc = gc[gc.cluster.map(lambda c: bool(elig.get(c, False)))]
        gc, gmiss = attach(gc)
        gunits = select(gc, ["phenotype_id"], "X")
        stats["gtex_crossover_units"] = int(len(gunits))
        stats["gtex_candidate_rows_without_spliceai"] = gmiss

    # ---- signal groups: CSs in the same cluster sharing >= 1 member (all members, columns 0-7 only)
    members = defaultdict(set)
    ukey = set(zip(units.phenotype_id, units.cs_id))
    with gzip.open(SUSIE, "rt") as f:
        f.readline()
        for line in f:
            p = line.split("\t", 8)
            if (p[0], p[6]) in ukey:
                members[(p[0], p[6])].add(p[1])
    uf = UF()
    by_cluster = defaultdict(list)
    for r in units.itertuples():
        by_cluster[r.cluster].append(r.unit_id)
        uf.find(r.unit_id)
    mem_u = {r.unit_id: members[(r.phenotype_id, r.cs_id)] for r in units.itertuples()}
    for c, us in by_cluster.items():
        for i in range(len(us)):
            for j in range(i + 1, len(us)):
                if mem_u[us[i]] & mem_u[us[j]]:
                    uf.union(us[i], us[j])
    units["signal_group"] = [uf.find(u) for u in units.unit_id]

    # ---- components (per rule)
    c2g = pd.read_csv(CL2GENE, sep="\t")
    genes_of = {c: set(str(g).split(",")) for c, g in zip(c2g.clu, c2g.genes)}

    def components(u):
        for rule in ("P", "D"):
            uf = UF()
            link = defaultdict(list)
            for r in u.itertuples():
                uf.find(r.unit_id)
                link[("cl", r.cluster)].append(r.unit_id)
                for g in genes_of.get(r.cluster, {r.cluster}):
                    link[("g", g)].append(r.unit_id)
                link[("v", getattr(r, "variant_id_" + rule))].append(r.unit_id)
            for us in link.values():
                for x in us[1:]:
                    uf.union(us[0], x)
            pos = u[["unit_id", "chrom", "pos_" + rule]].sort_values(["chrom", "pos_" + rule]).to_numpy()
            for a, b in zip(pos[:-1], pos[1:]):
                if a[1] == b[1] and b[2] - a[2] <= LD_PROXY_BP:
                    uf.union(a[0], b[0])
            u["component_" + rule] = [uf.find(x) for x in u.unit_id]

    components(units)
    for rule in ("P", "D"):
        stats["components_" + rule] = int(units["component_" + rule].nunique())
    stats["signal_groups"] = int(units.signal_group.nunique())
    if gunits is not None and len(gunits):
        components(gunits)
        stats["gtex_crossover_components_P"] = int(gunits.component_P.nunique())

    # ---- baseline count matrices (restricted to annotated introns of unit clusters)
    ucl = sorted(set(units.cluster) | (set(gunits.cluster) if gunits is not None else set()))
    loc = {}
    for c in ucl:
        ci = cl_index[c]
        a, n = K.starts[ci], K.size[ci]
        idx = [i for i in range(a, a + n) if ann[i]]
        # patterns re-indexed onto annotated introns only
        pmap = {i - a: k for k, i in enumerate(idx)}
        pats = [tuple(pmap[j] for j in p) for p in pats_all[ci]]
        loc[c] = dict(idx=np.array(idx), A=pattern_matrix(pats, len(idx)), pats=pats)
    all_idx = np.concatenate([loc[c]["idx"] for c in ucl])

    man = pd.read_csv(D.MAN, sep="\t")
    sfrac = {}
    for e in man.iloc[:, 1]:
        p = D.JDIR / f"{e}.stats.json"
        if p.exists():
            s = json.loads(p.read_text())
            sfrac[e] = s["spliced_reads"] / s["reads_or_pairs"]
    bs, _ = D.load_brainspan(W)
    bs_std, _ = D.load_brainspan(W, regions=D.STANDARD_REGIONS)
    bs_sf = load_brainspan_sfrac(W, sfrac)
    gt, tis, _ = D.load_gtex(W)
    ages = D.donor_ages()
    pcw = {d: int(a.replace("PCW", "")) for d, a in ages.items() if a.endswith("PCW")}
    F1 = sorted(d for d in bs if d in pcw and pcw[d] <= 26 and d not in D.WEN_OVERLAP)
    F1_all = sorted(d for d in bs if d in pcw and d not in D.WEN_OVERLAP)
    A1 = sorted(d for d in bs if d in D.ADULT)
    G = list(tis.tissue)
    Gnb = list(tis.tissue[tis.is_brain == 0])

    def mat(src, us):
        us = [u for u in us if u in src]
        return np.vstack([src[u][all_idx] for u in us]) if us else np.zeros((0, len(all_idx)), np.int64)

    groups = {"F1": mat(bs, F1), "A1": mat(bs, A1), "G": mat(gt, G),
              "F1_all": mat(bs, F1_all), "G_nonbrain": mat(gt, Gnb),
              "F1_std": mat(bs_std, F1), "A1_std": mat(bs_std, A1),
              "F1_sf": mat(bs_sf, F1), "A1_sf": mat(bs_sf, A1)}
    stats["units_per_group"] = {g: int(m.shape[0]) for g, m in groups.items()}
    off, o = {}, 0
    for c in ucl:
        off[c] = (o, o + len(loc[c]["idx"]))
        o += len(loc[c]["idx"])

    # pooled totals for thinning target (primary trio)
    def pooled(Mg, c):
        a, b = off[c]
        return Mg[:, a:b].sum()

    T = {c: min(pooled(groups[g], c) for g in ("F1", "A1", "G")) for c in ucl}

    # ---- weights per baseline, per seed
    weights = {}   # (group, cluster) -> list of weight vectors (per seed) ; 'U' and '*_nothin' single
    for g, Mg in groups.items():
        for c in ucl:
            a, b = off[c]
            Cc = Mg[:, a:b]
            tot = Cc.sum()
            ws = []
            for k in range(N_SEEDS):
                rng = np.random.default_rng(SEED0 + k + 7919 * hash_cluster(c))
                frac = min(1.0, T[c] / tot) if tot > 0 else 0.0
                Ct = rng.binomial(Cc, frac) if frac < 1 else Cc
                ws.append(fit_weights(unit_mean_usage(Ct), loc[c]["A"]))
            weights[(g, c)] = ws
            if g in ("F1", "A1", "G"):
                weights[(g + "_nothin", c)] = [fit_weights(unit_mean_usage(Cc), loc[c]["A"])]
    for c in ucl:
        npat = loc[c]["A"].shape[1]
        weights[("U", c)] = [np.full(npat, 1.0 / npat)]

    # ---- predictions
    def site_positions(c, strand):
        """per local annotated intron: (donor_pos, acceptor_pos)"""
        out = []
        for i in loc[c]["idx"]:
            s, e = keys[i][1], keys[i][2]
            out.append((s - 1, e + 1) if strand == "+" else (e + 1, s - 1))
        return out

    def pattern_factors(c, strand, vpos, row, eps):
        sites = site_positions(c, strand)
        acc_r, acc_a = tr["acc_ref"][row], tr["acc_alt"][row]
        don_r, don_a = tr["don_ref"][row], tr["don_alt"][row]

        def ratio(pos, kind):
            j = pos - vpos + WIN_TRACK
            if j < 0 or j > 2 * WIN_TRACK:
                return 1.0
            r, a = (don_r[j], don_a[j]) if kind == "d" else (acc_r[j], acc_a[j])
            return (a + eps) / (r + eps)

        f = []
        for p in loc[c]["pats"]:
            ss = {(sites[i][0], "d") for i in p} | {(sites[i][1], "a") for i in p}
            f.append(np.prod([ratio(pos, k) for pos, k in ss]))
        return np.array(f)

    def predict(w, A, f, li):
        if w is None or not np.all(np.isfinite(w)):
            return np.nan
        w2 = w * f
        w2 = w2 / w2.sum()
        return usage_from_weights(w2, A)[li] - usage_from_weights(w, A)[li]

    base_names = list(groups) + ["F1_nothin", "A1_nothin", "G_nothin", "U"]

    def mean_or_nan(vals):
        vals = np.array(vals, dtype=float)
        return float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan

    def run(u):
        out_cols = defaultdict(list)
        for r in u.itertuples():
            c = r.cluster
            A = loc[c]["A"]
            li = [keys[i] for i in loc[c]["idx"]].index((r.ph_chrom, r.ph_start, r.ph_end))
            for rule in ("P", "D"):
                vpos, row = getattr(r, "pos_" + rule), getattr(r, "row_" + rule)
                f = pattern_factors(c, r.strand, vpos, row, EPS)
                for g in base_names:
                    out_cols[f"pred_{g}" + ("" if rule == "P" else "_D")].append(
                        mean_or_nan([predict(w, A, f, li) for w in weights[(g, c)]]))
                if rule == "P":
                    for eps in EPS_SENS:
                        f2 = pattern_factors(c, r.strand, vpos, row, eps)
                        for g in ("F1", "A1", "G", "U"):
                            out_cols[f"pred_{g}_eps{eps}"].append(
                                mean_or_nan([predict(w, A, f2, li) for w in weights[(g, c)]]))
            # outcome-blind baseline difference for secondary analysis 8
            wFs = [w for w in weights[("F1", c)] if w is not None]
            wAs = [w for w in weights[("A1", c)] if w is not None]
            if wFs and wAs:
                wF, wA = np.mean(wFs, axis=0), np.mean(wAs, axis=0)
                out_cols["cluster_max_abs_usage_diff_F1_A1"].append(
                    float(np.max(np.abs(usage_from_weights(wF, A) - usage_from_weights(wA, A)))))
            else:
                out_cols["cluster_max_abs_usage_diff_F1_A1"].append(np.nan)
        for k, v in out_cols.items():
            u[k] = v
        return sorted(out_cols)

    stats["prediction_columns"] = run(units)
    units.to_csv(RES / "predictions.tsv", sep="\t", index=False)
    stats["units_with_all_primary_predictions"] = int(units[["pred_F1", "pred_A1", "pred_G", "pred_U"]].notna().all(1).sum())
    if gunits is not None and len(gunits):
        run(gunits)
        gunits.to_csv(RES / "predictions_gtex_crossover.tsv", sep="\t", index=False)
        stats["gtex_units_with_all_primary_predictions"] = int(
            gunits[["pred_F1", "pred_A1", "pred_G", "pred_U"]].notna().all(1).sum())
    (RES / "predictions_stats.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps({k: v for k, v in stats.items() if k != "prediction_columns"}, indent=1))


def hash_cluster(c):
    """Deterministic integer from a cluster ID (Python's hash() is salted per process)."""
    return int.from_bytes(c.encode(), "little") % 1_000_003


def load_brainspan_sfrac(W, sfrac):
    """BrainSpan donors' counts excluding samples with spliced-read fraction < SFRAC_MIN."""
    key = {(c, s, e): i for i, (c, s, e) in enumerate(zip(W.chrom, W.start, W.end))}
    man = pd.read_csv(D.MAN, sep="\t")
    by = {}
    for entry in man.iloc[:, 1]:
        if sfrac.get(entry, 0) < SFRAC_MIN:
            continue
        p = D.JDIR / f"{entry}.junctions.tsv.gz"
        donor = "HSB" + entry.split("_")[1]
        v = by.setdefault(donor, np.zeros(len(W), dtype=np.int64))
        with gzip.open(p, "rt") as f:
            f.readline()
            for line in f:
                c, s, e, r = line.split("\t")
                i = key.get((c, int(s), int(e)))
                if i is not None:
                    v[i] += int(r)
    return by


if __name__ == "__main__":
    main()

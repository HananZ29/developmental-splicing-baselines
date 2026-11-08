"""E12a (pre-registered): GSE310131 iPSC -> NGN2-neuron long-read time course.
Integrity, transcript mapping, overlap, local coverage at the 147 Wen sQTL clusters, junction-level
support of lead introns (hg38 -> hg19 lift), and reproducibility across lineages and wells.

Stage labels: gse_t00_iPSC, gse_t04_intermediate, gse_t30_NGN2neuron (NOT fetal / adult brain).
Replication: one donor; 3 lineages (iPSC stocks) = highest-level replicate; wells within a lineage
are derivative/technical replicates.

Inputs:
  data/external/GSE310131_full_expression.tsv/GSE310131_full_expression.tsv   (GEO, MD5 e722a601...)
  data/external/gse310131_zenodo/extracted/V47/                                (Zenodo 16897791, V47.zip MD5 5c845f36...)
Output: publication_extension/results/E12_gse310131/
"""
import gzip
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from isoform_proportions import load_gtf_transcripts, load_sqanti, long_counts  # noqa: E402
from propagation import intron_membership, smooth_baseline  # noqa: E402
from splice_graph import load_lift37  # noqa: E402

EXT = ROOT / "data" / "external"
Z = EXT / "gse310131_zenodo" / "extracted" / "V47"
GEO = EXT / "GSE310131_full_expression.tsv" / "GSE310131_full_expression.tsv"
OUT = ROOT / "publication_extension" / "results" / "E12_gse310131"
STAGE = {"iPSC": "t00", "NPC": "t04", "CN": "t30"}
ATTR = re.compile(r'transcript_id "([^"]+)"')


def sample_info(cols):
    rows = []
    for c in cols:
        p = c.split("_")
        rows.append(dict(sample=c, stage=STAGE[p[0]], lineage=int(p[1]), well=(p[2] if len(p) > 2 else "1")))
    return pd.DataFrame(rows)


def lift37_all_ensts():
    s = set()
    with gzip.open(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if "\ttranscript\t" in line:
                m = re.search(r'transcript_id "(ENST\d+)', line)
                if m:
                    s.add(m.group(1))
    return s


def gse_junctions(expr, stage_cols):
    """hg38 junctions (chrom, strand, first intronic base, last intronic base) -> reads per stage, all isoforms."""
    ex = {}
    with open(Z / "final_transcripts.gtf") as fh:
        for line in fh:
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon":
                continue
            t = ATTR.search(f[8]).group(1)
            ex.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
    st = expr.set_index("isoform")
    sums = {s: st[c].sum(1) for s, c in stage_cols.items()}
    junc = {}
    for t, (c, strand, e) in ex.items():
        e.sort()
        r = {s: float(sums[s].get(t, 0.0)) for s in stage_cols}
        for a, b in zip(e, e[1:]):
            k = (c, strand, a[1] + 1, b[0] - 1)
            j = junc.setdefault(k, {s: 0.0 for s in stage_cols} | {"n_iso": 0})
            for s in stage_cols:
                j[s] += r[s]
            j["n_iso"] += 1
    return junc, ex


def lift_junctions(junc):
    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))
    out, fail = {}, 0
    for (c, strand, s, e), v in junc.items():
        a, b = lo.convert_coordinate(c, s - 1), lo.convert_coordinate(c, e - 1)
        if not a or not b or a[0][0] != b[0][0] or (b[0][1] - a[0][1]) != (e - s):
            fail += 1
            continue
        key = (a[0][0], a[0][1] + 1, b[0][1] + 1)
        prev = out.get(key)
        if prev is None:
            out[key] = dict(v)
        else:
            for kk in v:
                prev[kk] += v[kk]
    return out, fail


def usage(pi, m, comp, k):
    inc = m @ pi
    den = inc[comp[k]].sum()
    return inc[k] / den if den > 0 else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    S = {}
    # ---------------- E12a integrity
    geo = pd.read_csv(GEO, sep="\t")
    full = pd.read_parquet(Z / "full_expression.parquet")
    fin = pd.read_parquet(Z / "final_expression.parquet")
    cls = pd.read_parquet(Z / "final_classification.parquet")
    cols = [c for c in geo.columns if c != "isoform"]
    si = sample_info(cols)
    g2, f2 = geo.set_index("isoform").sort_index(), full.set_index("isoform").sort_index()
    fsub = f2.reindex(fin.isoform)
    passes = ((f2 > 5).sum(1) >= 3)
    gtf_ids = set()
    with open(Z / "final_transcripts.gtf") as fh:
        for line in fh:
            if "\ttranscript\t" in line:
                gtf_ids.add(ATTR.search(line).group(1))
    S["integrity"] = {
        "geo_rows": len(geo), "geo_samples": cols, "geo_equals_zenodo_full": bool(g2.index.equals(f2.index) and (g2[cols].to_numpy() == f2[cols].to_numpy()).all()),
        "final_rows": len(fin), "final_subset_of_full_with_identical_counts": bool(fsub.notna().all().all() and (fsub[cols].to_numpy() == fin.set_index("isoform")[cols].to_numpy()).all()),
        "full_rows_passing_gt5_in_ge3_samples": int(passes.sum()), "final_rows_passing_that_rule": int(passes.reindex(fin.isoform).sum()),
        "final_with_gtf_structure": int(fin.isoform.isin(gtf_ids).sum()), "final_with_classification": int(fin.isoform.isin(cls.isoform).sum()),
        "reads_per_sample_full": {c: int(f2[c].sum()) for c in cols}, "reads_per_sample_final": {c: int(fin[c].sum()) for c in cols},
        "min_nonzero_count_full": int(f2[cols].to_numpy()[f2[cols].to_numpy() > 0].min()),
        "samples": si.to_dict(orient="records"),
    }
    del geo, g2, f2, full, fsub

    # ---------------- mapping
    cls["enst"] = cls.associated_transcript.astype(str).str.split(".").str[0]
    fin = fin.merge(cls[["isoform", "structural_category", "enst", "associated_gene", "chrom"]], on="isoform", how="left")
    stage_cols = {st: si[si.stage == st]["sample"].tolist() for st in ["t00", "t04", "t30"]}
    cat = fin.groupby("structural_category")[cols].sum()
    S["read_share_by_category_per_stage"] = {st: (cat[c].sum(1) / cat[c].sum(1).sum()).round(4).to_dict() for st, c in stage_cols.items()}
    S["isoforms_by_category"] = fin.structural_category.value_counts().to_dict()
    l37 = lift37_all_ensts()
    fsm = fin[fin.structural_category == "full-splice_match"]
    ism = fin[fin.structural_category == "incomplete-splice_match"]
    u_fsm = set(fsm.enst)
    S["fsm_mapping"] = {"fsm_isoforms": len(fsm), "unique_v47_ensts": len(u_fsm), "ensts_in_v31lift37": len(u_fsm & l37),
                        "ensts_not_in_v31lift37": len(u_fsm - l37),
                        "fsm_reads_on_unmappable_ensts_frac": round(float(fsm[~fsm.enst.isin(l37)][cols].to_numpy().sum() / fsm[cols].to_numpy().sum()), 4)}

    def agg(df, samples):
        return df.groupby("enst")[samples].sum().sum(1).to_dict()
    raw = {}
    for st, c in stage_cols.items():
        raw[f"gse_{st}"] = agg(fsm, c)
        raw[f"gse_{st}_fsm_ism"] = agg(pd.concat([fsm, ism]), c)
        for L in (1, 2, 3):
            sc = si[(si.stage == st) & (si.lineage == L)]["sample"].tolist()
            raw[f"gse_{st}_L{L}"] = agg(fsm, sc)
    for c in cols:
        raw[f"well_{c}"] = agg(fsm, [c])
    bt = pd.DataFrame({k: pd.Series(v) for k, v in raw.items()}).fillna(0).rename_axis("enst").reset_index()
    bt.to_csv(OUT / "gse_baseline_counts.tsv.gz", sep="\t", index=False)

    # Leung (counts, as E7) and Patowary (counts)
    gtf38 = load_gtf_transcripts(EXT / "gencode" / "gencode.v31.annotation.gtf.gz")
    leung_samples = {}
    for stage in ["Fetal", "Adult"]:
        c_, samples = load_sqanti(ROOT / "data" / "raw" / f"{stage}CTX_classification.tsv.gz")
        cnt = long_counts(c_, samples, gtf38, "exclude").reset_index()
        cnt = cnt[cnt.feature != "NOVEL"]
        cnt["enst"] = cnt.feature.str.split(".").str[0]
        raw[stage.lower()] = dict(zip(cnt.enst, cnt[samples].sum(1)))
        for s_ in samples:
            raw[f"leung_{s_}"] = dict(zip(cnt.enst, cnt[s_]))
        leung_samples[stage] = samples
    s2b = pd.read_csv(ROOT / "data" / "processed" / "patowary_S2B_counts.tsv.gz", sep="\t")
    libs = [c for c in s2b.columns if c[:3] in ("209", "334", "336")]
    s2b["enst"] = s2b.transcript_id.astype(str).str.split(".").str[0]
    kn = s2b[s2b.enst.str.startswith("ENST")]
    raw["patowary"] = dict(zip(kn.enst, kn[libs].sum(1)))

    # ---------------- universe / overlap
    pooled = pd.read_csv(ROOT / "results" / "proportions" / "fetal_vs_adult_pooled.tsv.gz", sep="\t")
    cand = set(pooled[(pooled.fetal_gene_reads >= 20) & (pooled.adult_gene_reads >= 20)].gene_id)
    bg = load_lift37(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", cand)
    auto = {f"chr{i}" for i in range(1, 23)}
    uni = {g: ts for g, ts in bg.items() if len(ts) >= 2 and ts[0].chrom in auto}
    ev = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    pos = ev[(ev.arm == "sQTL") & (ev.role == "positive")].reset_index(drop=True)
    uni_tx = {t.enst for ts in uni.values() for t in ts}

    def has(src, s):
        return {t for t in s if raw[src].get(t, 0) > 0}
    gene_reads = {st: {g: sum(raw[f"gse_{st}"].get(t.enst, 0) for t in ts) for g, ts in uni.items()} for st in stage_cols}
    ov = {"universe_genes": len(uni), "universe_transcripts": len(uni_tx)}
    for st in stage_cols:
        ov[f"universe_genes_with_reads_{st}"] = int(sum(v > 0 for v in gene_reads[st].values()))
        ov[f"universe_genes_ge20_reads_{st}"] = int(sum(v >= 20 for v in gene_reads[st].values()))
        ov[f"universe_transcripts_with_reads_{st}"] = len(has(f"gse_{st}", uni_tx))
    ov["universe_transcripts_with_reads_any_stage"] = len(set().union(*[has(f"gse_{st}", uni_tx) for st in stage_cols]))
    for src in ["fetal", "adult", "patowary"]:
        ov[f"universe_transcripts_with_reads_{src}"] = len(has(src, uni_tx))
    gse_any = set().union(*[has(f"gse_{st}", uni_tx) for st in stage_cols])
    ov["universe_tx_reads_in_leung_fetal_and_gse_any"] = len(has("fetal", uni_tx) & gse_any)
    ov["universe_tx_reads_in_leung_adult_and_gse_any"] = len(has("adult", uni_tx) & gse_any)
    ov["universe_tx_reads_in_patowary_and_gse_any"] = len(has("patowary", uni_tx) & gse_any)
    ov["universe_tx_reads_in_gse_any_not_leung_fetal"] = len(gse_any - has("fetal", uni_tx))
    ov["universe_tx_reads_in_gse_any_not_patowary"] = len(gse_any - has("patowary", uni_tx))
    for arm in ["sQTL", "isoQTL"]:
        gs = set(ev[ev.arm == arm].ensg)
        ov[f"{arm}_evaluation_genes"] = len(gs)
        ov[f"{arm}_evaluation_genes_with_t30_reads"] = int(sum(gene_reads["t30"].get(g, 0) > 0 for g in gs))
    ov["sQTL_positive_genes"] = pos.ensg.nunique()
    ov["sQTL_positive_genes_ge20_reads_t30"] = int(sum(gene_reads["t30"].get(g, 0) >= 20 for g in set(pos.ensg)))
    S["overlap"] = ov

    # ---------------- junctions
    junc, ex = gse_junctions(fin[["isoform"] + cols], stage_cols)
    j19, fail = lift_junctions(junc)
    S["junctions"] = {"hg38_unique_junctions": len(junc), "lift_failures": fail, "hg19_junctions": len(j19)}

    # ---------------- per-variant local coverage
    by_gene = load_lift37(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", set(pos.ensg))
    srcs = ["gse_t00", "gse_t04", "gse_t30", "fetal", "adult", "patowary"]
    rows, cluster_tx = [], {}
    for r in pos.itertuples():
        ts = sorted(by_gene[r.ensg], key=lambda t: t.enst)
        names = [t.enst for t in ts]
        introns, index, m, comp = intron_membership(ts, "component")
        s, e = map(int, r.target.split("-"))
        k = index[(s, e)]
        in_cluster = m[comp[k]].sum(axis=0) > 0
        tgt = m[k] > 0
        cluster_tx[r.variant_id + "|" + r.ensg] = [n for n, f in zip(names, in_cluster) if f]
        o = dict(variant_id=r.variant_id, ensg=r.ensg, target=r.target, z_obs=r.z_obs, n_tx_cluster=int(in_cluster.sum()))
        for src in srcs:
            vals = np.array([raw[src].get(n, 0.0) for n in names])
            o[f"{src}_cluster_reads"] = float(vals[in_cluster].sum())
            o[f"{src}_target_reads"] = float(vals[tgt].sum())
            o[f"{src}_cluster_tx_zero"] = int((vals[in_cluster] == 0).sum())
            o[f"{src}_u_baseline"] = usage(smooth_baseline(raw[src], names), m, comp, k)
        chrom = ts[0].chrom
        ck = [(chrom, *introns[j]) for j in comp[k]]
        tj = j19.get((chrom, s, e), {})
        o["lead_junction_in_gse"] = bool(tj)
        for st in stage_cols:
            tr_ = tj.get(st, 0.0)
            tot = sum(j19.get(c_, {}).get(st, 0.0) for c_ in ck)
            o[f"gse_{st}_junction_target_reads"] = tr_
            o[f"gse_{st}_junction_cluster_reads"] = tot
            o[f"gse_{st}_u_junction"] = tr_ / tot if tot > 0 else np.nan
        rows.append(o)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "per_variant.tsv", sep="\t", index=False)
    cov = {"n_variants": len(d), "cluster_transcripts_total": len({t for v in cluster_tx.values() for t in v})}
    alltx = sorted({t for v in cluster_tx.values() for t in v})
    for src in srcs:
        cov[f"cluster_tx_zero_reads_{src}"] = int(sum(raw[src].get(t, 0) == 0 for t in alltx))
        cov[f"lead_intron_unsupported_transcript_level_{src}"] = int((d[f"{src}_target_reads"] == 0).sum())
        cov[f"median_cluster_reads_{src}"] = float(d[f"{src}_cluster_reads"].median())
    for st in stage_cols:
        cov[f"lead_intron_zero_junction_reads_{st}"] = int((d[f"gse_{st}_junction_target_reads"] == 0).sum())
        cov[f"lead_intron_unsupported_both_transcript_and_junction_{st}"] = int(((d[f"gse_{st}_target_reads"] == 0) & (d[f"gse_{st}_junction_target_reads"] == 0)).sum())
    cov["lead_junction_present_in_gse_any_isoform"] = int(d.lead_junction_in_gse.sum())
    cov["lead_intron_unsupported_both_leung"] = int(((d.fetal_target_reads == 0) & (d.adult_target_reads == 0)).sum())
    cov["lead_intron_unsupported_t00_and_t30"] = int(((d.gse_t00_target_reads == 0) & (d.gse_t30_target_reads == 0)).sum())
    cov["spearman_u_baseline_vs_u_junction_t30"] = round(float(spearmanr(d.gse_t30_u_baseline, d.gse_t30_u_junction, nan_policy="omit").statistic), 3)
    S["local_coverage"] = cov

    # stage contrasts at lead introns
    def q(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        return {"n": int(len(x)), **{str(p): round(float(np.quantile(x, p)), 4) for p in [.5, .75, .9, .99]}}
    sc = {}
    for a, b in [("gse_t30", "gse_t00"), ("gse_t30", "gse_t04"), ("gse_t04", "gse_t00"), ("fetal", "adult"), ("gse_t30", "fetal"), ("gse_t30", "adult")]:
        dif = (d[f"{a}_u_baseline"] - d[f"{b}_u_baseline"]).abs()
        sup = (d[f"{a}_target_reads"] > 0) & (d[f"{b}_target_reads"] > 0)
        sc[f"{a} vs {b}"] = {"all": q(dif), "lead_supported_in_both": q(dif[sup])}
    sc["junction_usage t30 vs t00 (both clusters >=20 reads)"] = q((d.gse_t30_u_junction - d.gse_t00_u_junction).abs()[(d.gse_t30_junction_cluster_reads >= 20) & (d.gse_t00_junction_cluster_reads >= 20)])
    S["lead_intron_stage_contrast_abs_u"] = sc

    # ---------------- reproducibility (cluster transcripts; abundance strata as E7)
    df = pd.DataFrame({"enst": alltx})
    gene_of = {t.enst: g for g, ts in by_gene.items() for t in ts}
    df["gene"] = df.enst.map(gene_of)
    keys = [f"gse_t30_L{L}" for L in (1, 2, 3)] + [f"gse_t04_L{L}" for L in (1, 2, 3)] + [f"gse_t00_L{L}" for L in (1, 2, 3)] + \
           [f"well_{c}" for c in cols] + ["gse_t00", "gse_t04", "gse_t30", "fetal", "adult", "patowary"] + [f"leung_{s_}" for s_ in leung_samples["Fetal"][:2]]
    for k_ in keys:
        df[k_] = df.enst.map(raw[k_]).fillna(0.0)
        tot = df.groupby("gene")[k_].transform("sum")
        df[f"p_{k_}"] = np.where(tot > 0, df[k_] / tot.replace(0, np.nan), np.nan)
    df["mean_t30_lineage_reads"] = df[[f"gse_t30_L{L}" for L in (1, 2, 3)]].mean(1)
    df["abund_bin"] = pd.cut(df.mean_t30_lineage_reads, [-1, 0, 9, 99, np.inf], labels=["0", "1-9", "10-99", ">=100"])

    def agree(a, b, sub):
        x, y = sub[f"p_{a}"], sub[f"p_{b}"]
        ok = x.notna() & y.notna()
        if ok.sum() < 5:
            return np.nan, np.nan, int(ok.sum())
        return round(float(np.corrcoef(x[ok], y[ok])[0, 1]), 3), round(float(spearmanr(x[ok], y[ok]).statistic), 3), int(ok.sum())
    wells = {st: {L: si[(si.stage == st) & (si.lineage == L)]["sample"].tolist() for L in (1, 2, 3)} for st in ["t04", "t30"]}
    pairs = [("lineage", st, f"gse_{st}_L{a}", f"gse_{st}_L{b}") for st in ["t00", "t04", "t30"] for a, b in itertools.combinations((1, 2, 3), 2)] + \
            [("well_within_lineage", st, f"well_{w[0]}", f"well_{w[1]}") for st in ["t04", "t30"] for w in wells[st].values()] + \
            [("stage", "", "gse_t30", "gse_t00"), ("stage", "", "gse_t30", "gse_t04"), ("cross_study", "", "gse_t30", "fetal"), ("cross_study", "", "gse_t30", "adult"),
             ("cross_study", "", "gse_t30", "patowary"), ("leung_replicate", "", f"leung_{leung_samples['Fetal'][0]}", f"leung_{leung_samples['Fetal'][1]}")]
    out_rows = []
    for bn, sub in [("all", df)] + [(str(k_), g_) for k_, g_ in df.groupby("abund_bin", observed=True)]:
        for kind, st, a, b in pairs:
            pr, sr, n = agree(a, b, sub)
            out_rows.append(dict(abundance_bin_mean_t30_lineage_reads=bn, kind=kind, stage=st, pair=f"{a}~{b}", pearson=pr, spearman=sr, n=n))
    ta = pd.DataFrame(out_rows)
    ta.to_csv(OUT / "transcript_agreement.tsv", sep="\t", index=False)
    summ_ag = ta.groupby(["abundance_bin_mean_t30_lineage_reads", "kind", "stage"]).agg(pearson_min=("pearson", "min"), pearson_max=("pearson", "max"),
                                                                                          spearman_min=("spearman", "min"), spearman_max=("spearman", "max"), n=("n", "max")).reset_index()
    S["cluster_transcript_agreement_summary"] = summ_ag.to_dict(orient="records")
    S["cluster_transcripts_by_t30_abundance_bin"] = df.abund_bin.value_counts().sort_index().astype(int).to_dict()

    # genome-wide lineage structure check on all final isoforms (log CPM)
    X = fin[cols].to_numpy(float)
    L = np.log2(X / X.sum(0) * 1e6 + 1)
    C = np.corrcoef(L.T)
    idx = {c: i for i, c in enumerate(cols)}
    within, between = {}, {}
    for st in ["t04", "t30"]:
        ss = si[si.stage == st]
        w_, b_ = [], []
        for a, b in itertools.combinations(ss.itertuples(), 2):
            (w_ if a.lineage == b.lineage else b_).append(C[idx[a.sample], idx[b.sample]])
        within[st], between[st] = round(float(np.mean(w_)), 4), round(float(np.mean(b_)), 4)
    S["lineage_check_logCPM_corr"] = {"mean_within_lineage_well_pairs": within, "mean_between_lineage_pairs": between,
                                      "t00_lineage_pairs": [round(float(C[idx[a], idx[b]]), 4) for a, b in itertools.combinations(stage_cols["t00"], 2)],
                                      "t30_vs_t00_mean": round(float(np.mean([C[idx[a], idx[b]] for a in stage_cols["t30"] for b in stage_cols["t00"]])), 4)}

    # genome-wide (universe, smoothed) agreement for contrast with E6
    names_u = {g: sorted(t.enst for t in ts) for g, ts in uni.items()}

    def vec(c):
        return np.concatenate([smooth_baseline(raw[c], names_u[g]) for g in names_u])
    V = {k_: vec(k_) for k_ in ["gse_t30_L1", "gse_t30_L2", "gse_t30_L3", "gse_t00_L1", "gse_t00_L2", "gse_t00_L3", "gse_t00", "gse_t04", "gse_t30", "fetal", "adult", "patowary"]}
    S["universe_smoothed_r"] = {f"{a}~{b}": round(float(np.corrcoef(V[a], V[b])[0, 1]), 4) for a, b in
                                [("gse_t30_L1", "gse_t30_L2"), ("gse_t30_L1", "gse_t30_L3"), ("gse_t30_L2", "gse_t30_L3"),
                                 ("gse_t00_L1", "gse_t00_L2"), ("gse_t00_L1", "gse_t00_L3"), ("gse_t00_L2", "gse_t00_L3"),
                                 ("gse_t30", "gse_t00"), ("gse_t30", "gse_t04"), ("gse_t04", "gse_t00"), ("fetal", "adult"),
                                 ("gse_t30", "fetal"), ("gse_t30", "adult"), ("gse_t30", "patowary"), ("gse_t00", "fetal")]}
    S["universe_nonzero_weight_fraction"] = {k_: round(float(np.mean([raw[k_].get(t, 0) > 0 for g in names_u for t in names_u[g]])), 4)
                                             for k_ in ["gse_t00", "gse_t04", "gse_t30", "fetal", "adult", "patowary"]}
    (OUT / "summary_E12a.json").write_text(json.dumps(S, indent=2, default=str))
    print(json.dumps(S, indent=2, default=str))


if __name__ == "__main__":
    main()

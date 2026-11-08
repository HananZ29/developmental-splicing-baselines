"""E7: variant-local audit of baseline coverage and agreement for the 147 Wen sQTL positives.

For each positive, the quantities that actually drive its prediction are the baseline weights of the
transcripts in the lead intron's LeafCutter-style cluster K(i). We measure, per source:
  - reads on cluster transcripts; cluster transcripts with zero reads ("unsupported")
  - the baseline intron usage u_i (smoothed exactly as in the experiment)
and, from Patowary's full TALON annotation (known + novel + ISM transcripts, Bulk.gtf, hg19), the
junction-level read support of the lead intron and its competitors, giving an annotation-independent
observed usage u_junction per donor.

Transcript-level agreement is then stratified by abundance, so that agreement of abundant
transcripts is not mistaken for agreement of low-expression transcripts.

Output: publication_extension/results/E7_local_coverage/{per_variant.tsv, transcript_agreement.tsv, summary.json}
"""
import gzip
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from isoform_proportions import load_gtf_transcripts, load_sqanti, long_counts  # noqa: E402
from propagation import intron_membership, smooth_baseline  # noqa: E402
from splice_graph import load_lift37  # noqa: E402

OUT = ROOT / "publication_extension" / "results" / "E7_local_coverage"
DONORS = ["209", "334", "336"]
ATTR = re.compile(r'(gene_id|transcript_id) "([^"]+)"')


def usage(pi, m, comp, k):
    inc = m @ pi
    den = inc[comp[k]].sum()
    return inc[k] / den if den > 0 else np.nan


def bulk_junction_reads(genes, s2b, libs):
    """Per-junction read support (summed over all TALON transcripts containing it) per donor."""
    tx = {}
    with open(ROOT / "data" / "external" / "patowary" / "DataS1" / "Bulk.gtf") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            if a["gene_id"].split(".")[0] not in genes:
                continue
            tx.setdefault(a["transcript_id"], [f[0], []])[1].append((int(f[3]), int(f[4])))
    reads = s2b.set_index("transcript_id")
    donor_cols = {d: [c for c in libs if c.startswith(d)] for d in DONORS}
    junc = {}
    for t, (c, ex) in tx.items():
        if t not in reads.index:
            continue
        r = {d: float(reads.loc[t, donor_cols[d]].sum()) for d in DONORS}
        ex.sort()
        for a, b in zip(ex, ex[1:]):
            key = (c, a[1] + 1, b[0] - 1)
            j = junc.setdefault(key, {d: 0.0 for d in DONORS} | {"n_tx": 0, "n_known": 0})
            for d in DONORS:
                j[d] += r[d]
            j["n_tx"] += 1
            j["n_known"] += t.startswith("ENST")
    return junc


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ev = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    pos = ev[(ev.arm == "sQTL") & (ev.role == "positive")].reset_index(drop=True)
    genes = set(pos.ensg)
    by_gene = load_lift37(ROOT / "data" / "external" / "gencode" / "gencode.v31lift37.annotation.gtf.gz", genes)

    # raw counts per ENST for every source
    s2b = pd.read_csv(ROOT / "data" / "processed" / "patowary_S2B_counts.tsv.gz", sep="\t")
    libs = [c for c in s2b.columns if c[:3] in DONORS]
    s2b["enst"] = s2b.transcript_id.astype(str).str.split(".").str[0]
    known = s2b[s2b.enst.str.startswith("ENST")]
    raw = {f"d{d}": dict(zip(known.enst, known[[c for c in libs if c.startswith(d)]].sum(1))) for d in DONORS}
    raw["patowary"] = dict(zip(known.enst, known[libs].sum(1)))
    gtf = load_gtf_transcripts(ROOT / "data" / "external" / "gencode" / "gencode.v31.annotation.gtf.gz")
    leung_samples = {}
    for stage in ["Fetal", "Adult"]:
        c, samples = load_sqanti(ROOT / "data" / "raw" / f"{stage}CTX_classification.tsv.gz")
        cnt = long_counts(c, samples, gtf, "exclude").reset_index()
        cnt = cnt[cnt.feature != "NOVEL"]
        cnt["enst"] = cnt.feature.str.split(".").str[0]
        raw[stage.lower()] = dict(zip(cnt.enst, cnt[samples].sum(1)))
        for s in samples:
            raw[f"leung_{s}"] = dict(zip(cnt.enst, cnt[s]))
        leung_samples[stage] = samples

    junc = bulk_junction_reads(genes, s2b, libs)

    rows, cluster_tx = [], {}
    for r in pos.itertuples():
        ts = sorted(by_gene[r.ensg], key=lambda t: t.enst)
        names = [t.enst for t in ts]
        introns, index, m, comp = intron_membership(ts, "component")
        s, e = map(int, r.target.split("-"))
        k = index[(s, e)]
        in_cluster = m[comp[k]].sum(axis=0) > 0
        tgt = m[k] > 0
        cl_names = [n for n, f in zip(names, in_cluster) if f]
        cluster_tx[(r.variant_id, r.ensg, r.target)] = cl_names
        out = dict(variant_id=r.variant_id, ensg=r.ensg, target=r.target, z_obs=r.z_obs,
                   n_tx_gene=len(names), n_tx_cluster=int(in_cluster.sum()), n_tx_target=int(tgt.sum()),
                   n_introns_cluster=len(comp[k]))
        for src in ["fetal", "adult", "patowary"] + [f"d{d}" for d in DONORS]:
            vals = np.array([raw[src].get(n, 0.0) for n in names])
            out[f"{src}_cluster_reads"] = float(vals[in_cluster].sum())
            out[f"{src}_target_reads"] = float(vals[tgt].sum())
            out[f"{src}_cluster_tx_zero"] = int((vals[in_cluster] == 0).sum())
            out[f"{src}_u_baseline"] = usage(smooth_baseline(raw[src], names), m, comp, k)
        # junction-level observed usage (all Patowary transcripts, incl. novel/ISM)
        chrom = ts[0].chrom
        comp_keys = [(chrom, *introns[j]) for j in comp[k]]
        for d in DONORS:
            tr = junc.get((chrom, s, e), {}).get(d, 0.0)
            tot = sum(junc.get(c_, {}).get(d, 0.0) for c_ in comp_keys)
            out[f"d{d}_junction_target_reads"] = tr
            out[f"d{d}_junction_cluster_reads"] = tot
            out[f"d{d}_u_junction"] = tr / tot if tot > 0 else np.nan
        tj = junc.get((chrom, s, e), {})
        out["target_junction_in_patowary_annotation"] = bool(tj)
        out["target_junction_n_talon_tx"] = tj.get("n_tx", 0)
        out["target_junction_n_known_tx"] = tj.get("n_known", 0)
        rows.append(out)
    d = pd.DataFrame(rows)
    d["donor_u_range"] = d[[f"d{x}_u_baseline" for x in DONORS]].max(1) - d[[f"d{x}_u_baseline" for x in DONORS]].min(1)
    d["fetal_adult_u_absdiff"] = (d.fetal_u_baseline - d.adult_u_baseline).abs()
    d["patowary_leung_fetal_u_absdiff"] = (d.patowary_u_baseline - d.fetal_u_baseline).abs()
    d["patowary_u_baseline_minus_u_junction_abs"] = (d.patowary_u_baseline - d[[f"d{x}_u_junction" for x in DONORS]].mean(1)).abs()
    d["pat_bin"] = pd.cut(d.patowary_cluster_reads, [-1, 0, 19, 199, np.inf], labels=["0", "1-19", "20-199", ">=200"])
    d["leung_fetal_bin"] = pd.cut(d.fetal_cluster_reads, [-1, 0, 19, 199, np.inf], labels=["0", "1-19", "20-199", ">=200"])
    d.to_csv(OUT / "per_variant.tsv", sep="\t", index=False)

    def strat(col):
        g = d.groupby(col, observed=True)
        return pd.DataFrame({
            "n_variants": g.size(),
            "median_donor_u_range": g.donor_u_range.median(),
            "median_fetal_adult_u_absdiff": g.fetal_adult_u_absdiff.median(),
            "median_patowary_vs_leung_fetal_u_absdiff": g.patowary_leung_fetal_u_absdiff.median(),
            "median_baseline_vs_junction_usage_absdiff": g.patowary_u_baseline_minus_u_junction_abs.median(),
            "frac_target_unsupported_leung_fetal": g.apply(lambda x: float((x.fetal_target_reads == 0).mean()), include_groups=False),
            "frac_target_unsupported_patowary": g.apply(lambda x: float((x.patowary_target_reads == 0).mean()), include_groups=False),
        }).round(4)
    by_pat = strat("pat_bin")
    by_leung = strat("leung_fetal_bin")

    # transcript-level agreement in the 147 clusters, stratified by abundance
    tx = sorted({t for v in cluster_tx.values() for t in v})
    df = pd.DataFrame({"enst": tx})
    for src in [f"d{x}" for x in DONORS] + [f"leung_{s}" for s in leung_samples["Fetal"]] + ["fetal", "adult", "patowary"]:
        df[src] = df.enst.map(raw[src]).fillna(0.0)
    gene_of = {t.enst: g for g, ts in by_gene.items() for t in ts}
    df["gene"] = df.enst.map(gene_of)

    def within_gene(col):
        tot = df.groupby("gene")[col].transform("sum")
        return np.where(tot > 0, df[col] / tot.replace(0, np.nan), np.nan)
    for c in [c for c in df.columns if c not in ("enst", "gene")]:
        df[f"p_{c}"] = within_gene(c)
    df["mean_donor_reads"] = df[[f"d{x}" for x in DONORS]].mean(1)
    df["abund_bin"] = pd.cut(df.mean_donor_reads, [-1, 0, 9, 99, np.inf], labels=["0", "1-9", "10-99", ">=100"])

    def agree(a, b, sub):
        x, y = sub[f"p_{a}"], sub[f"p_{b}"]
        ok = x.notna() & y.notna()
        if ok.sum() < 5:
            return (np.nan, np.nan, int(ok.sum()))
        return (round(float(np.corrcoef(x[ok], y[ok])[0, 1]), 3), round(float(spearmanr(x[ok], y[ok]).statistic), 3), int(ok.sum()))
    rows2 = []
    pairs = [("d209", "d334"), ("d209", "d336"), ("d334", "d336"), ("patowary", "fetal"), ("fetal", "adult"),
             tuple(f"leung_{s}" for s in leung_samples["Fetal"][:2])]
    for b_, sub in [("all", df)] + [(str(k), g) for k, g in df.groupby("abund_bin", observed=True)]:
        for a, c in pairs:
            pr, sr, n = agree(a, c, sub)
            rows2.append(dict(abundance_bin_mean_donor_reads=b_, pair=f"{a}~{c}", pearson=pr, spearman=sr, n_transcripts=n))
    ta = pd.DataFrame(rows2)
    ta.to_csv(OUT / "transcript_agreement.tsv", sep="\t", index=False)

    summ = {
        "n_variants": len(d),
        "cluster_transcripts_total": len(tx),
        "cluster_transcripts_zero_reads": {s: int((df[s] == 0).sum()) for s in ["fetal", "adult", "patowary", "d209", "d334", "d336"]},
        "target_intron_unsupported_variants": {s: int((d[f"{s}_target_reads"] == 0).sum()) for s in ["fetal", "adult", "patowary", "d209", "d334", "d336"]},
        "target_junction_present_in_patowary_annotation": int(d.target_junction_in_patowary_annotation.sum()),
        "target_junction_reads_zero_all_donors": int((d[[f"d{x}_junction_target_reads" for x in DONORS]].sum(1) == 0).sum()),
        "median_cluster_reads": {s: float(d[f"{s}_cluster_reads"].median()) for s in ["fetal", "adult", "patowary"]},
        "variants_by_patowary_cluster_reads": by_pat.to_dict(orient="index"),
        "variants_by_leung_fetal_cluster_reads": by_leung.to_dict(orient="index"),
        "median_donor_u_range_all": float(d.donor_u_range.median()),
        "median_fetal_adult_u_absdiff_all": float(d.fetal_adult_u_absdiff.median()),
        "spearman_donor_u_across_variants": {f"{a}~{b}": round(float(spearmanr(d[f"d{a}_u_baseline"], d[f"d{b}_u_baseline"]).statistic), 3)
                                              for a, b in itertools.combinations(DONORS, 2)},
        "spearman_u_baseline_vs_u_junction_patowary": round(float(spearmanr(d.patowary_u_baseline, d[[f"d{x}_u_junction" for x in DONORS]].mean(1), nan_policy="omit").statistic), 3),
    }
    (OUT / "summary.json").write_text(json.dumps(summ, indent=2, default=str))
    print(json.dumps(summ, indent=2, default=str))
    print(ta.to_string())


if __name__ == "__main__":
    main()

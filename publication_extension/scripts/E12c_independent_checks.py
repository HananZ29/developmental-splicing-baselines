"""E12c: independent checks of E12a/E12b, using separate code paths (no src/ helpers).

(a) Liftover + FSM mapping: for a random sample of FSM isoforms (seed 20261104), lift every GTF junction
    hg38 -> hg19 and test it is an intron of the associated GENCODE v31lift37 transcript.
(b) Stage baseline counts: recompute gse_t30 / gse_t00 per-ENST FSM counts directly from the parquet
    files and compare with gse_baseline_counts.tsv.gz.
(c) Lead-intron transcript-level support at t30: recompute from the raw GTF text for the 147 Wen sQTL
    positives and compare with E12a (44 unsupported).
(d) Lifted lead junctions: hg19 dinucleotides at lifted GSE junctions that equal a lead intron (GT..AG
    expected on the gene strand).
Output: publication_extension/results/E12_gse310131/independent_checks.json
"""
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
Z = EXT / "gse310131_zenodo" / "extracted" / "V47"
OUT = ROOT / "publication_extension" / "results" / "E12_gse310131"
TID = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'gene_id "([^".]+)')


def lift37_introns(want_tx=None, want_genes=None):
    ex = {}
    gene = {}
    with gzip.open(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            t = TID.search(f[8]).group(1).split(".")[0]
            g = GID.search(f[8]).group(1)
            if (want_tx is not None and t not in want_tx) and (want_genes is None or g not in want_genes):
                continue
            ex.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
            gene[t] = g
    intr = {}
    for t, (c, s, e) in ex.items():
        e.sort()
        intr[t] = (c, s, {(a[1] + 1, b[0] - 1) for a, b in zip(e, e[1:])})
    return intr, gene


def main():
    res = {}
    cls = pd.read_parquet(Z / "final_classification.parquet", columns=["isoform", "structural_category", "associated_transcript"])
    fin = pd.read_parquet(Z / "final_expression.parquet")
    cols = [c for c in fin.columns if c != "isoform"]
    fsm = cls[cls.structural_category == "full-splice_match"].copy()
    fsm["enst"] = fsm.associated_transcript.str.replace(r"\.\d+$", "", regex=True)

    # (b)
    m = fin.merge(fsm[["isoform", "enst"]], on="isoform")
    t30 = m[[c for c in cols if c.startswith("CN_")]].sum(1).groupby(m.enst).sum()
    t00 = m[[c for c in cols if c.startswith("iPSC_")]].sum(1).groupby(m.enst).sum()
    bt = pd.read_csv(OUT / "gse_baseline_counts.tsv.gz", sep="\t").set_index("enst")
    res["b_counts_max_abs_diff"] = {"t30": float((t30 - bt.gse_t30.reindex(t30.index)).abs().max()),
                                    "t00": float((t00 - bt.gse_t00.reindex(t00.index)).abs().max()), "n_enst": int(len(t30))}

    # (a)
    rng = np.random.default_rng(20261104)
    multi = fsm.isoform.to_numpy()
    samp = set(rng.choice(multi, 2000, replace=False))
    iso_ex = {}
    with open(Z / "final_transcripts.gtf") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            t = TID.search(f[8]).group(1)
            if t in samp:
                iso_ex.setdefault(t, [f[0], []])[1].append((int(f[3]), int(f[4])))
    iso2enst = dict(zip(fsm.isoform, fsm.enst))
    intr, _ = lift37_introns(want_tx={iso2enst[i] for i in samp})
    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))
    n_j = n_ok = n_fail = 0
    iso_all_ok = iso_mappable = 0
    for iso, (c, e) in iso_ex.items():
        enst = iso2enst[iso]
        if enst not in intr or len(e) < 2:
            continue
        iso_mappable += 1
        e.sort()
        ok_all = True
        for a, b in zip(e, e[1:]):
            s, t = a[1] + 1, b[0] - 1
            x, y = lo.convert_coordinate(c, s - 1), lo.convert_coordinate(c, t - 1)
            n_j += 1
            if not x or not y:
                n_fail += 1
                ok_all = False
                continue
            hit = (x[0][1] + 1, y[0][1] + 1) in intr[enst][2]
            n_ok += hit
            ok_all &= hit
        iso_all_ok += ok_all
    res["a_fsm_liftover"] = {"sampled_fsm": 2000, "multi_exon_and_in_lift37": iso_mappable, "junctions": n_j,
                             "lift_fail": n_fail, "junctions_matching_lift37_intron": n_ok,
                             "frac_junctions_matching": round(n_ok / max(n_j, 1), 4), "isoforms_all_junctions_match": iso_all_ok}

    # (c)
    ev = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    pos = ev[(ev.arm == "sQTL") & (ev.role == "positive")]
    gi, gene = lift37_introns(want_tx=set(), want_genes=set(pos.ensg))
    by_gene = {}
    for t, g in gene.items():
        by_gene.setdefault(g, []).append(t)
    unsup = 0
    for r in pos.itertuples():
        s, e = map(int, r.target.split("-"))
        txs = [t for t in by_gene[r.ensg] if (s, e) in gi[t][2]]
        if sum(float(t30.get(t, 0.0)) for t in txs) == 0:
            unsup += 1
    res["c_lead_intron_unsupported_t30_transcript_level"] = unsup

    # (d)
    pv = pd.read_csv(OUT / "per_variant.tsv", sep="\t")
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    strand = {r.ensg: gi[by_gene[r.ensg][0]][1] for r in pos.itertuples()}
    dn = []
    for r in pv[pv.lead_junction_in_gse].itertuples():
        s, e = map(int, r.target.split("-"))
        c = gi[by_gene[r.ensg][0]][0]
        d5, d3 = fa[c][s - 1:s + 1].seq.upper(), fa[c][e - 2:e].seq.upper()
        dn.append((d5 + ".." + d3) if strand[r.ensg] == "+" else (d5 + ".." + d3))
    vc = pd.Series(dn).value_counts().to_dict()
    res["d_lead_junction_dinucleotides_hg19"] = vc
    res["d_canonical_frac"] = round(float(sum(v for k, v in vc.items() if k in ("GT..AG", "CT..AC", "GC..AG", "CT..GC")) / max(len(dn), 1)), 4)
    (OUT / "independent_checks.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()

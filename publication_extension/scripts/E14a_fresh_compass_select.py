"""E14a (pre-registered in generalized_model_design/PREREGISTRATION_E14.md): select the never-touched
COMPASS test set ("T-fresh"). Uses counts and coordinates only; no outcome is computed.

Rules (identical QC to E13c):
- Single SNVs with a WT row; >= 20 reads in both replicates for variant and WT, in HEK and K562.
- hg38 -> hg19 liftover with exon length preserved; REF checked against hg19.
- Exclude every event and gene in the E13c subset (train and test).
- Exclude events whose exon is within +-5 bp of any Vex-seq or MFASS exon (identical exons).
- Exclude events whose exon +-100 nt has 8-mer containment >= 0.15 in any E13c exon +-100 nt (a
  paralog proxy for >= ~80% identity; the §4 T2 rule).
- Keep events with >= 5 QC singles; at most CAP variants per event (seed 20261106).

Output: publication_extension/results/E14_comparators/fresh_variants.tsv, fresh_to_score.tsv,
fresh_selection.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
OUT = ROOT / "publication_extension" / "results" / "E14_comparators"
SEED, MIN_READS, MIN_SINGLES, K, CONTAIN = 20261106, 20, 5, 8, 0.15
CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
COMP = str.maketrans("ACGT", "TGCA")


def reads(d, c, r):
    return d[f"{c}_rep{r}_included"] + d[f"{c}_rep{r}_excluded"]


def window(fa, chrom, s, e, strand):
    x = fa[chrom][s - 1 - 100:e + 100].seq.upper()
    return x.translate(COMP)[::-1] if strand == "-" else x


def kmers(x):
    return {x[i:i + K] for i in range(len(x) - K + 1)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(EXT / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt = d[d.snp == "none"].drop_duplicates("event_id").set_index("event_id")
    s = d[(d.seq_type == "single") & d.event_id.isin(wt.index)].copy()
    ok = np.ones(len(s), bool)
    for c in ["HEK", "K562"]:
        for r in (1, 2):
            ok &= (reads(s, c, r) >= MIN_READS).to_numpy()
            ok &= (s.event_id.map(reads(wt, c, r)) >= MIN_READS).to_numpy()
    s = s[ok].copy()
    s["gene"] = s.gene_exon.str.split(" ").str[0]
    e13 = pd.read_csv(ROOT / "publication_extension" / "results" / "E13_compass" / "subset_variants.tsv", sep="\t")
    n0 = s.event_id.nunique()
    s = s[~s.event_id.isin(e13.event_id) & ~s.gene.isin(e13.gene)].copy()
    n1 = s.event_id.nunique()
    p = s.variant_hg38.str.extract(r"^(chr[\w]+):(\d+):([ACGT])>([ACGT])$")
    s["chrom38"], s["pos38"], s["ref"], s["alt"] = p[0], pd.to_numeric(p[1]), p[2], p[3]
    s = s[s.pos38.notna()].copy()
    s["strand"] = s.event_id.str[-1]
    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))

    def lift(c, x):
        r = lo.convert_coordinate(c, int(x) - 1)
        return (r[0][0], r[0][1] + 1) if r else (None, None)
    L = [lift(c, x) for c, x in zip(s.chrom38, s.pos38)]
    s["chrom"], s["pos"] = [a for a, _ in L], [b for _, b in L]
    s["exon_start_hg19"] = [lift(c, x)[1] for c, x in zip(s.chrom38, s.exon_start_hg38)]
    s["exon_end_hg19"] = [lift(c, x)[1] for c, x in zip(s.chrom38, s.exon_end_hg38)]
    s = s[s.chrom.notna() & s.exon_start_hg19.notna() & s.exon_end_hg19.notna()].copy()
    for col in ["pos", "exon_start_hg19", "exon_end_hg19"]:
        s[col] = s[col].astype(int)
    s = s[(s.exon_end_hg19 - s.exon_start_hg19) == (s.exon_end_hg38 - s.exon_start_hg38)]
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    s["hg19_base"] = [fa[c][q - 1:q].seq.upper() for c, q in zip(s.chrom, s.pos)]
    ref_match = float((s.hg19_base == s.ref).mean())
    s = s[s.hg19_base == s.ref].copy()
    n2 = s.event_id.nunique()

    # identical exons (+-5 bp) to Vex-seq or MFASS: same exon definitions and overlap() as E13a
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from E13a_label_audit import overlap  # noqa: E402
    m = pd.read_csv(EXT / "mfass" / "snv_data_clean.txt", sep="\t")
    m = m[m.category != "natural"] if "category" in m else m
    me = m[m.v2_dpsi.notna()].drop_duplicates("ensembl_id")
    mf_ex = pd.DataFrame({"chrom": me.chr,
                          "start": np.where(me.strand == "+", me.start + me.intron1_len, me.end - me.intron1_len - me.exon_len),
                          "end": np.where(me.strand == "+", me.start + me.intron1_len + me.exon_len, me.end - me.intron1_len)})
    vx = pd.read_csv(ROOT / "publication_extension" / "results" / "E10_vexseq" / "vexseq_measured.tsv", sep="\t")
    ve = vx.drop_duplicates(["chrom", "Exon up", "Exon down"])
    vx_ex = pd.DataFrame({"chrom": ve.chrom, "start": ve[["Exon up", "Exon down"]].min(1), "end": ve[["Exon up", "Exon down"]].max(1)})
    both = pd.concat([mf_ex, vx_ex]).reset_index(drop=True)
    ev = s.drop_duplicates("event_id")[["event_id", "chrom", "strand", "exon_start_hg19", "exon_end_hg19"]].reset_index(drop=True)
    q = pd.DataFrame({"chrom": ev.chrom, "start": ev.exon_start_hg19, "end": ev.exon_end_hg19})
    ident_ev = set(ev.event_id[overlap(q, both)])

    # 8-mer containment against all E13c exons
    e13ev = e13.drop_duplicates("event_id")
    ref_k = [kmers(window(fa, r.chrom, r.exon_start_hg19, r.exon_end_hg19, r.strand)) for r in e13ev.itertuples()]
    sim_ev, max_c = set(), {}
    for r in ev.itertuples():
        q = kmers(window(fa, r.chrom, r.exon_start_hg19, r.exon_end_hg19, r.strand))
        c = max(len(q & t) / max(len(q), 1) for t in ref_k)
        max_c[r.event_id] = c
        if c >= CONTAIN:
            sim_ev.add(r.event_id)
    s = s[~s.event_id.isin(ident_ev | sim_ev)]
    cnt = s.groupby("event_id").size()
    keep = cnt[cnt >= MIN_SINGLES].index
    s = s[s.event_id.isin(keep)]
    s = s.groupby("event_id", group_keys=False).apply(lambda g: g.sample(min(len(g), CAP), random_state=SEED))
    s["split"] = "fresh"
    s.drop(columns=[c for c in s.columns if c in ("intron1", "exon", "intron2", "full_seq")]).to_csv(
        OUT / "fresh_variants.tsv", sep="\t", index=False)
    s[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().to_csv(OUT / "fresh_to_score.tsv", sep="\t", index=False)
    rep = {"cap_per_family": CAP, "events_qc_any": n0, "events_after_excluding_E13c_events_genes": n1,
           "events_after_lift_ref": n2, "ref_match_hg19": round(ref_match, 4),
           "excluded_identical_to_vexseq_mfass": len(ident_ev), "excluded_8mer_similar_to_E13c": len(sim_ev),
           "max_8mer_containment_quantiles": [round(float(x), 3) for x in np.quantile(list(max_c.values()), [0.5, 0.9, 0.99, 1])],
           "families": int(s.event_id.nunique()), "genes": int(s.gene.nunique()), "variants": len(s),
           "unique_to_score": int(s[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().shape[0])}
    (OUT / "fresh_selection.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()

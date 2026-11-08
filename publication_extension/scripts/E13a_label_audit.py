"""E13a (pre-registered): audit of measured variant-effect labels on disk and their mutual overlap.

Measured variant-effect datasets: MFASS (HEK293T), Vex-seq (K562, HepG2), COMPASS (HEK293, HeLa, K562,
MCF7, HMC3). Exon overlap is computed in hg19 (COMPASS lifted from hg38).
Output: publication_extension/results/E13_compass/label_audit.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
OUT = ROOT / "publication_extension" / "results" / "E13_compass"
CELLS = ["HEK", "HeLa", "K562", "MCF7", "HMC3"]


def overlap(a, b, tol=5):
    """a, b: DataFrames chrom,start,end -> number of a-intervals overlapping any b-interval (+-tol)."""
    hit = np.zeros(len(a), bool)
    for c, g in a.groupby("chrom"):
        bb = b[b.chrom == c]
        if bb.empty:
            continue
        bs, be = bb.start.to_numpy() - tol, bb.end.to_numpy() + tol
        for i, (s, e) in zip(g.index, zip(g.start, g.end)):
            hit[a.index.get_loc(i)] = bool(((bs <= e) & (be >= s)).any())
    return hit


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    R = {}
    # ---- MFASS
    m = pd.read_csv(EXT / "mfass" / "snv_data_clean.txt", sep="\t")
    mm = m[m.category != "natural"] if "category" in m else m
    meas = mm[mm.v2_dpsi.notna()]
    R["MFASS"] = {"context": ["HEK293T"], "replicates": "2 (R1, R2) per library; natural (reference) index per replicate",
                  "rows": len(m), "variant_rows_with_v2_dpsi": len(meas), "unique_variants": int(meas[["chr", "snp_position", "alt_allele"]].drop_duplicates().shape[0]),
                  "exons": int(meas.ensembl_id.nunique()), "chromosomes": int(meas.chr.nunique()),
                  "reference_psi_per_replicate": bool({"nat_v2_index_R1", "nat_v2_index_R2"} <= set(m.columns)),
                  "abs_dpsi_ge_0.1": int((meas.v2_dpsi.abs() >= 0.1).sum()),
                  "median_reference_inclusion": round(float(meas.nat_v2_index.median()), 3)}
    me = meas.drop_duplicates("ensembl_id")
    mf_ex = pd.DataFrame({"chrom": me.chr, "start": np.where(me.strand == "+", me.start + me.intron1_len, me.end - me.intron1_len - me.exon_len),
                          "end": np.where(me.strand == "+", me.start + me.intron1_len + me.exon_len, me.end - me.intron1_len)}).reset_index(drop=True)
    # ---- Vex-seq
    v = pd.read_csv(ROOT / "publication_extension" / "results" / "E10_vexseq" / "vexseq_measured.tsv", sep="\t")
    R["Vex-seq"] = {"context": ["K562", "HepG2"], "replicates": "3 per cell line; reference exon PSI per replicate",
                    "snvs": len(v), "exons": int(v[["chrom", "Exon up", "Exon down"]].drop_duplicates().shape[0]),
                    "abs_dpsi_ge_0.1": {"K562": int((v.K562_dpsi.abs() >= .1).sum()), "HepG2": int((v.HepG2_dpsi.abs() >= .1).sum())},
                    "median_reference_psi": {"K562": round(float(v.drop_duplicates('Exon up').K562_base.median()), 3), "HepG2": round(float(v.drop_duplicates('Exon up').HepG2_base.median()), 3)}}
    ve = v.drop_duplicates(["chrom", "Exon up", "Exon down"])
    vx_ex = pd.DataFrame({"chrom": ve.chrom, "start": ve[["Exon up", "Exon down"]].min(1), "end": ve[["Exon up", "Exon down"]].max(1)}).reset_index(drop=True)
    # ---- COMPASS
    d = pd.read_csv(EXT / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz")
    wt_ev = set(d[d.snp == "none"].event_id)
    single = d[d.seq_type == "single"]
    R["COMPASS"] = {"context": CELLS, "replicates": "2 per cell line (HEK: 2 more WT-only replicates)",
                    "rows": len(d), "wt_rows": int((d.snp == "none").sum()), "single_rows": len(single), "double_rows": int((d.seq_type == "double").sum()),
                    "events_exons": int(d.event_id.nunique()), "events_with_wt": len(wt_ev), "genes": int(d.gene_exon.str.split(" ").str[0].nunique()),
                    "single_with_wt": int(single.event_id.isin(wt_ev).sum()),
                    "sources": d.source.value_counts().to_dict(),
                    "reference_psi_per_replicate": True,
                    "note": "reference (WT) PSI measured in the same transfections/replicates as variants; replicate split required to avoid shared noise"}
    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))
    ce = d.drop_duplicates("event_id")
    rows = []
    for r in ce.itertuples():
        c = r.event_id.split(":")[0]
        a, b = lo.convert_coordinate(c, int(r.exon_start_hg38) - 1), lo.convert_coordinate(c, int(r.exon_end_hg38) - 1)
        if a and b and a[0][0] == b[0][0]:
            rows.append((a[0][0], min(a[0][1], b[0][1]) + 1, max(a[0][1], b[0][1]) + 1, r.gene_exon.split(" ")[0]))
    cp_ex = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])
    R["COMPASS"]["events_lifted_hg19"] = len(cp_ex)
    R["overlap_identical_exon_(+-5bp)"] = {
        "COMPASS_exons_overlapping_MFASS": int(overlap(cp_ex, mf_ex).sum()),
        "COMPASS_exons_overlapping_VexSeq": int(overlap(cp_ex, vx_ex).sum()),
        "VexSeq_exons_overlapping_MFASS": int(overlap(vx_ex, mf_ex).sum()),
    }

    def near(a, b, w):
        return int(overlap(a, b, tol=w).sum())
    R["overlap_within_100kb"] = {"COMPASS_exons_near_MFASS": near(cp_ex, mf_ex, 100_000), "COMPASS_exons_near_VexSeq": near(cp_ex, vx_ex, 100_000),
                                 "VexSeq_exons_near_MFASS": near(vx_ex, mf_ex, 100_000)}
    R["expression_only_resources_(no_variant_effects)"] = {
        "Leung 2021 Iso-Seq": "fetal/adult cortex transcript counts", "Patowary 2024": "fetal neocortex long-read counts, 3 donors",
        "GSE310131": "iPSC->NGN2 neuron long-read counts, 1 donor", "GTEx v8 RSEM": "adult tissue transcript TPM"}
    R["association_resources_(not_measured_effects)"] = {"Wen 2024": "prenatal cortex sQTL/isoQTL leads (standardized z)", "GTEx v8 Brain_Cortex sQTL": "adult cortex sQTL leads"}
    (OUT / "label_audit.json").write_text(json.dumps(R, indent=2, default=int))
    print(json.dumps(R, indent=2, default=int))


if __name__ == "__main__":
    main()

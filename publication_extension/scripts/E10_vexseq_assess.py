"""E10: feasibility and overlap assessment of Vex-seq (Adamson et al. 2018; GEO GSE113163) BEFORE any
benchmark. No predictions are computed here.

Checks: hg19 REF match; replicate agreement; read depth; baseline (reference-exon) PSI in K562 vs
HepG2 and how much the baseline differs between the two contexts; overlap with evaluated QTL loci.
Output: publication_extension/results/E10_vexseq/assessment.json, vexseq_measured.tsv
"""
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / "data" / "external" / "vexseq"
OUT = ROOT / "publication_extension" / "results" / "E10_vexseq"


def load(cell):
    reps = []
    for f in sorted(glob.glob(str(V / "geo" / f"GSM*_{cell}_rep*_PSI.tsv.gz"))):
        if "KD" in f or "barcode" in f:
            continue
        d = pd.read_csv(f, sep="\t").set_index("variant")
        d["reads"] = d.inclusion_reads + d.exclusion_reads
        reps.append(d)
    return reps


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = {}
    cells = {}
    for cell in ["K562", "HepG2"]:
        reps = load(cell)
        psi = pd.concat([r.PSI for r in reps], axis=1) / 100
        rd = pd.concat([r.reads for r in reps], axis=1)
        cells[cell] = pd.DataFrame({f"{cell}_psi": psi.mean(1), f"{cell}_psi_sd": psi.std(1), f"{cell}_reads": rd.sum(1)})
        c = psi.corr()
        out[f"{cell}_replicates"] = len(reps)
        out[f"{cell}_replicate_psi_r_min"] = round(float(c.values[np.triu_indices(len(reps), 1)].min()), 3)
    m = cells["K562"].join(cells["HepG2"], how="outer").reset_index()
    p = m.variant.str.split("_", expand=True)
    m["chrom"], m["pos"], m["ref"], m["alt"] = p[0], pd.to_numeric(p[1], errors="coerce"), p[2], p[3]
    ctrl = m[m.ref == "Control"].copy()
    var = m[(m.ref != "Control") & m.pos.notna()].copy()
    var["pos"] = var.pos.astype(int)
    ctrl["pos"] = ctrl.pos.astype(int)
    # reference (Original) exon PSI by exon key = chrom_pos of the control row
    orig = ctrl[ctrl.alt == "Original"].drop_duplicates(["chrom", "pos"]).set_index(["chrom", "pos"])
    des = pd.read_excel(V / "Vexseq_AdditionalFile1_design.xlsx")
    des = des[des.reference != "Control"].copy()
    des["pos"] = des.position.astype(int)
    des = des.drop_duplicates(["chromosome", "pos", "reference", "alternative"])
    var = var.merge(des[["chromosome", "pos", "reference", "alternative", "Exon up", "Exon down", "strand", "Upstream Int Adj"]],
                    left_on=["chrom", "pos", "ref", "alt"], right_on=["chromosome", "pos", "reference", "alternative"], how="left")
    # link each variant to its exon's Original row: control rows are keyed by the design's upstream intron start
    ctrl_des = pd.read_excel(V / "Vexseq_AdditionalFile1_design.xlsx")
    ctrl_des = ctrl_des[(ctrl_des.reference == "Control") & (ctrl_des.alternative == "Original")][["chromosome", "position", "Exon up", "Exon down"]]
    ctrl_des = ctrl_des.rename(columns={"position": "ctrl_pos"}).drop_duplicates(["chromosome", "Exon up", "Exon down"])
    var = var.merge(ctrl_des, on=["chromosome", "Exon up", "Exon down"], how="left")
    var = var.merge(orig[["K562_psi", "HepG2_psi"]].rename(columns={"K562_psi": "K562_base", "HepG2_psi": "HepG2_base"}),
                    left_on=["chrom", "ctrl_pos"], right_index=True, how="left")
    var["K562_dpsi"], var["HepG2_dpsi"] = var.K562_psi - var.K562_base, var.HepG2_psi - var.HepG2_base
    fa = Fasta(str(ROOT / "data" / "external" / "genome" / "hg19.fa"))
    snv = var[(var.ref.str.len() == 1) & (var.alt.str.len() == 1)].copy()
    snv["hg19"] = [fa[c][q - 1:q].seq.upper() if c in fa else "N" for c, q in zip(snv.chrom, snv.pos)]
    out["variants_total"] = len(var)
    out["snvs"] = len(snv)
    out["snv_ref_matches_hg19_plus_strand"] = round(float((snv.hg19 == snv.ref).mean()), 4)
    good = snv[(snv.hg19 == snv.ref) & snv.K562_dpsi.notna() & snv.HepG2_dpsi.notna() & snv.K562_base.notna()]
    out["snvs_measured_both_cells_with_baseline"] = len(good)
    out["exons"] = int(good[["chrom", "Exon up", "Exon down"]].drop_duplicates().shape[0])
    ex = good.drop_duplicates(["chrom", "Exon up", "Exon down"])
    out["baseline_psi_quantiles_K562"] = ex.K562_base.quantile([.1, .25, .5, .75, .9]).round(3).to_dict()
    out["baseline_psi_quantiles_HepG2"] = ex.HepG2_base.quantile([.1, .25, .5, .75, .9]).round(3).to_dict()
    bd = (ex.K562_base - ex.HepG2_base).abs()
    out["exon_baseline_absdiff_K562_vs_HepG2_quantiles"] = bd.quantile([.25, .5, .75, .9]).round(3).to_dict()
    out["exons_baseline_diff_ge_0.1"] = int((bd >= 0.1).sum())
    out["exons_baseline_intermediate_0.2_0.8_either_cell"] = int((ex.K562_base.between(.2, .8) | ex.HepG2_base.between(.2, .8)).sum())
    out["dpsi_corr_K562_vs_HepG2"] = round(float(good[["K562_dpsi", "HepG2_dpsi"]].corr().iloc[0, 1]), 3)
    out["abs_dpsi_ge_0.1_K562_HepG2"] = [int((good.K562_dpsi.abs() >= .1).sum()), int((good.HepG2_dpsi.abs() >= .1).sum())]
    ev = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    gt = pd.read_csv(ROOT / "publication_extension" / "results" / "E3_gtex" / "eligible_gtex.tsv", sep="\t")
    qpos = set(zip(ev.chrom, ev.pos)) | set(zip(gt.chrom, gt.pos))
    out["overlap_identical_variants_with_evaluated_QTL"] = int(len(qpos & set(zip(good.chrom, good.pos))))
    qreg = pd.concat([ev[ev.role == "positive"][["chrom", "pos"]], gt[["chrom", "pos"]]])
    near = sum(((qreg.chrom == c) & ((qreg.pos - q).abs() <= 5000)).any() for c, q in zip(ex.chrom, ex["Exon up"]))
    out["vexseq_exons_within_5kb_of_an_evaluated_QTL_lead"] = int(near)
    good.to_csv(OUT / "vexseq_measured.tsv", sep="\t", index=False)
    (OUT / "assessment.json").write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()

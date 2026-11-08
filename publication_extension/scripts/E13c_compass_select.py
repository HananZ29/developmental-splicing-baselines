"""E13c (pre-registered): select the COMPASS subset for SpliceAI scoring.

200 exon families (seed 20261105) with a WT row and >= 5 single SNVs passing QC (>= 20 reads, both
replicates, variant and WT) in HEK and K562; loci within 100 kb of any Vex-seq or MFASS exon excluded
(conservative same-gene/locus proxy). hg38 -> hg19 lift (UCSC chain), REF checked against hg19.
Output: publication_extension/results/E13_compass/{subset_variants.tsv, to_score.tsv, selection.json}
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
OUT = ROOT / "publication_extension" / "results" / "E13_compass"
SEED, N_FAM, MIN_READS, CAP = 20261105, 200, 20, 40


def reads(d, c, r, suffix=""):
    return d[f"{c}_rep{r}_included{suffix}"] + d[f"{c}_rep{r}_excluded{suffix}"]


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
    s = s[ok]
    p = s.variant_hg38.str.extract(r"^(chr[\w]+):(\d+):([ACGT])>([ACGT])$")
    s["chrom38"], s["pos38"], s["ref"], s["alt"] = p[0], pd.to_numeric(p[1]), p[2], p[3]
    s = s[s.pos38.notna()].copy()
    s["strand"] = s.event_id.str[-1]
    # exclusion: loci near Vex-seq / MFASS exons (hg19)
    vx = pd.read_csv(ROOT / "publication_extension" / "results" / "E10_vexseq" / "vexseq_measured.tsv", sep="\t")
    mf = pd.read_csv(EXT / "mfass" / "snv_data_clean.txt", sep="\t", usecols=["chr", "start", "end"])
    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))

    def lift(c, x):
        r = lo.convert_coordinate(c, int(x) - 1)
        return (r[0][0], r[0][1] + 1) if r else (None, None)
    L = [lift(c, x) for c, x in zip(s.chrom38, s.pos38)]
    s["chrom"], s["pos"] = [a for a, _ in L], [b for _, b in L]
    E = [lift(c, x) for c, x in zip(s.chrom38, s.exon_start_hg38)]
    F = [lift(c, x) for c, x in zip(s.chrom38, s.exon_end_hg38)]
    s["exon_start_hg19"], s["exon_end_hg19"] = [b for _, b in E], [b for _, b in F]
    s = s[s.chrom.notna() & s.exon_start_hg19.notna() & s.exon_end_hg19.notna()].copy()
    s["pos"], s["exon_start_hg19"], s["exon_end_hg19"] = s.pos.astype(int), s.exon_start_hg19.astype(int), s.exon_end_hg19.astype(int)
    s = s[(s.exon_end_hg19 - s.exon_start_hg19) == (s.exon_end_hg38 - s.exon_start_hg38)]
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    s["hg19_base"] = [fa[c][q - 1:q].seq.upper() for c, q in zip(s.chrom, s.pos)]
    ref_match = float((s.hg19_base == s.ref).mean())
    s = s[s.hg19_base == s.ref]
    near = set()
    other = pd.concat([pd.DataFrame({"chrom": vx.chrom, "pos": vx["Exon up"]}), pd.DataFrame({"chrom": mf.chr, "pos": mf.start})])
    oth = {c: np.sort(g.pos.to_numpy()) for c, g in other.groupby("chrom")}
    for ev, g in s.groupby("event_id"):
        c, x = g.chrom.iloc[0], g.exon_start_hg19.iloc[0]
        a = oth.get(c)
        if a is not None and len(a) and np.min(np.abs(a - x)) <= 100_000:
            near.add(ev)
    cnt = s.groupby("event_id").size()
    elig = sorted(set(cnt[cnt >= 5].index) - near)
    rng = np.random.default_rng(SEED)
    # deviation D-E13-1: only 155 families are eligible (< 200), so all are used; a cap of 40 randomly
    # chosen variants per family (seed fixed) keeps SpliceAI scoring ~2 h. Decided before any outcome.
    fam = list(rng.choice(elig, min(N_FAM, len(elig)), replace=False))
    sub = s[s.event_id.isin(fam)].copy()
    sub = sub.groupby("event_id", group_keys=False).apply(lambda g: g.sample(min(len(g), CAP), random_state=SEED))
    genes = sorted(sub.gene_exon.str.split(" ").str[0].unique())
    tr_genes = set(rng.choice(genes, len(genes) // 2, replace=False))
    sub["gene"] = sub.gene_exon.str.split(" ").str[0]
    sub["split"] = np.where(sub.gene.isin(tr_genes), "train", "test")
    sub.drop(columns=[c for c in sub.columns if c in ("intron1", "exon", "intron2", "full_seq")]).to_csv(OUT / "subset_variants.tsv", sep="\t", index=False)
    sub[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().to_csv(OUT / "to_score.tsv", sep="\t", index=False)
    rep = {"singles_with_wt_and_QC_HEK_K562": int(ok.sum()), "ref_match_hg19_after_lift": round(ref_match, 4),
           "events_eligible_ge5": int((cnt >= 5).sum()), "events_excluded_near_vexseq_mfass": len(near & set(cnt[cnt >= 5].index)),
           "families_sampled": len(fam), "variants_in_subset": len(sub), "unique_to_score": int(sub[["chrom", "pos", "ref", "alt", "strand"]].drop_duplicates().shape[0]),
           "genes": len(genes), "train_families": int(sub[sub.split == "train"].event_id.nunique()), "test_families": int(sub[sub.split == "test"].event_id.nunique())}
    (OUT / "selection.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()

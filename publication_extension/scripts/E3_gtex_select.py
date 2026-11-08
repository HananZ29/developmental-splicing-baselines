"""E3 (pre-registered): select adult GTEx v8 Brain_Cortex sQTL lead variants for the stage-swap test.

Same eligibility rules as the Wen (fetal) sQTL positives in experiments/06_select_variants.py.
GTEx is hg38; variants and introns are lifted to hg19 with the UCSC hg38ToHg19 chain (MD5-verified)
and the hg19 REF base must equal the GTEx REF. The intron coordinate convention is chosen by GT-AG
rate on hg19 before any prediction.

Output: publication_extension/results/E3_gtex/eligible_gtex.tsv, attrition.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from splice_graph import load_lift37, site_table  # noqa: E402

EXT = ROOT / "data" / "external"
OUT = ROOT / "publication_extension" / "results" / "E3_gtex"
WIN = 500
COMP = str.maketrans("ACGT", "TGCA")


def load06():
    spec = importlib.util.spec_from_file_location("sel06", ROOT / "experiments" / "06_select_variants.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def lift(lo, chrom, pos1):
    """1-based hg38 -> (chrom, 1-based hg19, strand) or None."""
    r = lo.convert_coordinate(chrom, pos1 - 1)
    if not r:
        return None
    c, p0, strand, _ = r[0]
    return c, p0 + 1, strand


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sel = load06()
    G0 = sel.gene_universe()
    by_gene = load_lift37(EXT / "gencode" / "gencode.v31lift37.annotation.gtf.gz", G0)
    G = {g for g, ts in by_gene.items() if len(ts) >= 2 and ts[0].chrom.replace("chr", "").isdigit()}
    sites = site_table({g: by_gene[g] for g in G})
    sites_by_gene = {g: np.sort(d["pos"].to_numpy()) for g, d in sites.groupby("ensg")}
    comp = {g: sel.intron_competitors(by_gene[g]) for g in G}
    info = {g: (by_gene[g][0].chrom, by_gene[g][0].strand) for g in G}

    sg = pd.read_csv(EXT / "drive_new" / "Brain_Cortex.v8.sgenes.txt.gz", sep="\t")
    sg["ensg"] = sg["gene_id"].str.split(".").str[0]
    att = {"sgenes_rows": len(sg), "sig_q05": int((sg.qval < 0.05).sum())}
    sig = sg[(sg.qval < 0.05) & sg.ensg.isin(G)].copy()
    att["sig_in_universe"] = len(sig)

    lo = LiftOver(str(EXT / "genome" / "hg38ToHg19.over.chain.gz"))
    fa = Fasta(str(EXT / "genome" / "hg19.fa"))
    p = sig.phenotype_id.str.split(":", expand=True)
    sig["i_chr"], sig["i_s38"], sig["i_e38"] = p[0], p[1].astype(int), p[2].astype(int)

    # lift introns; keep raw lifted start/end, convention decided below
    rows, why = [], {}
    for r in sig.itertuples():
        def bump(k):
            why[k] = why.get(k, 0) + 1
        if len(r.ref) != 1 or len(r.alt) != 1:
            bump("not_snv"); continue
        lv = lift(lo, r.chr, int(r.variant_pos))
        ls, le = lift(lo, r.i_chr, r.i_s38), lift(lo, r.i_chr, r.i_e38)
        if lv is None or ls is None or le is None:
            bump("liftover_failed"); continue
        g_chr, g_strand = info[r.ensg]
        if lv[0] != g_chr or ls[0] != g_chr or le[0] != g_chr or lv[2] != "+" or ls[2] != "+" or le[2] != "+":
            bump("liftover_chrom_or_strand_change"); continue
        if le[1] - ls[1] != r.i_e38 - r.i_s38:
            bump("intron_length_changed"); continue
        base = fa[g_chr][lv[1] - 1:lv[1]].seq.upper()
        if base != r.ref:
            bump("hg19_ref_mismatch"); continue
        rows.append(dict(ensg=r.ensg, gene_name=r.gene_name, variant_id=r.variant_id, rsid=r.rs_id_dbSNP151_GRCh38p7,
                         chrom=g_chr, pos=lv[1], ref=r.ref, alt=r.alt, strand=g_strand, s19=ls[1], e19=le[1],
                         qval=r.qval, z_obs=r.slope / r.slope_se, slope=r.slope, maf=r.maf, phenotype_id=r.phenotype_id))
    lifted = pd.DataFrame(rows)
    att["lift_outcomes"] = why | {"lifted_ok": len(lifted)}

    # pre-registered convention choice: highest GT-AG rate on hg19
    conv = {"start..end": (0, 0), "start+1..end-1": (1, -1), "start..end-1": (0, -1), "start+1..end": (1, 0)}
    rates = {}
    for name, (ds, de) in conv.items():
        ok = 0
        for r in lifted.itertuples():
            seq = fa[r.chrom][r.s19 + ds - 1:r.e19 + de].seq.upper()
            if r.strand == "-":
                seq = seq.translate(COMP)[::-1]
            ok += seq[:2] == "GT" and seq[-2:] == "AG"
        rates[name] = round(ok / len(lifted), 4)
    best = max(rates, key=rates.get)
    ds, de = conv[best]
    att["intron_convention_gtag_rates"] = rates
    att["intron_convention_chosen"] = best
    lifted["intron"] = list(zip(lifted.s19 + ds, lifted.e19 + de))

    out, why2 = [], {}
    for r in lifted.itertuples():
        nc = comp[r.ensg].get(r.intron, 0)
        if nc == 0:
            k = "intron_not_annotated"
        elif nc < 2:
            k = "no_competitor_intron"
        else:
            d = sel.nearest_site(sites_by_gene, r.ensg, int(r.pos))
            if d > WIN:
                k = "far_from_splice_site"
            else:
                k = "eligible"
                out.append(dict(arm="GTEx_sQTL", role="positive", ensg=r.ensg, variant_id=r.variant_id, chrom=r.chrom,
                                pos=int(r.pos), ref=r.ref, alt=r.alt, strand=r.strand, dist_site=int(d), qval=r.qval,
                                z_obs=r.z_obs, target=f"{r.intron[0]}-{r.intron[1]}", n_competitors=nc,
                                rsid=r.rsid, maf=r.maf))
        why2[k] = why2.get(k, 0) + 1
    att["eligibility_outcomes"] = why2
    elig = pd.DataFrame(out)
    att["final_positives"] = len(elig)
    att["final_genes"] = int(elig.ensg.nunique()) if len(elig) else 0

    wen = pd.read_csv(ROOT / "results" / "variants" / "eligible_variants.tsv", sep="\t")
    wen_pos = wen[(wen.arm == "sQTL") & (wen.role == "positive")]
    att["genes_shared_with_wen_positives"] = int(len(set(elig.ensg) & set(wen_pos.ensg)))
    att["same_lead_variant_pos_as_wen"] = int(len(set(zip(elig.chrom, elig.pos)) & set(zip(wen_pos.chrom, wen_pos.pos))))
    att["same_gene_same_target_intron_as_wen"] = int(len(set(zip(elig.ensg, elig.target)) & set(zip(wen_pos.ensg, wen_pos.target))))
    elig.to_csv(OUT / "eligible_gtex.tsv", sep="\t", index=False)
    (OUT / "attrition.json").write_text(json.dumps(att, indent=2, default=int))
    print(json.dumps(att, indent=2, default=int))


if __name__ == "__main__":
    main()

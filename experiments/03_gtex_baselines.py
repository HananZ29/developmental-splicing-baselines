"""Per-tissue GTEx v8 isoform proportions and the adult-cortex / generic baselines.

Streams the RSEM transcript TPM matrix (rows are grouped by gene) in chunks. For every sample, a
transcript's proportion is its TPM over the gene's summed TPM; samples where the gene sums to 0 are
undefined and skipped. Per tissue we keep the mean and median proportion and the number of samples
expressing the gene.

Baselines (transcript-level, GENCODE v26 ids; join on unversioned ENST):
- gtex_cortex:   Brain - Cortex
- gtex_ba9:      Brain - Frontal Cortex (BA9)   (closest region to the prefrontal SQANTI samples)
- generic_all:   unweighted mean over all 54 tissues of the tissue means (tissues expressing the gene)
- generic_nonbrain: same over the 41 non-brain tissues

Usage: python experiments/03_gtex_baselines.py
Outputs: results/gtex/tissue_mean_props.tsv.gz, tissue_median_props.tsv.gz, tissue_n_expr.tsv.gz,
         baselines.tsv.gz, validation.json
"""
import gzip
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GTEX = ROOT / "data" / "external" / "gtex"
TPM = GTEX / "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
ATTR = GTEX / "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
OUT = ROOT / "results" / "gtex"
CHUNK = 4000
warnings.filterwarnings("ignore", category=RuntimeWarning)  # all-NaN tissues for genes unexpressed there


def tissue_columns(header):
    attr = pd.read_csv(ATTR, sep="\t", usecols=["SAMPID", "SMTSD"]).set_index("SAMPID")["SMTSD"]
    samples = header[2:]
    tissue = attr.reindex(samples)
    groups = {t: np.flatnonzero((tissue == t).to_numpy()) for t in sorted(tissue.dropna().unique())}
    return groups


def summarize_block(block, groups):
    """block: rows of one or more complete genes. Returns per-tissue mean/median/n_expr frames."""
    genes = block["gene_id"].to_numpy()
    x = block.iloc[:, 2:].to_numpy(dtype=np.float32)
    tot = pd.DataFrame(x).groupby(genes).transform("sum").to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(tot > 0, x / tot, np.nan)
    means, medians, nexpr = {}, {}, {}
    for t, idx in groups.items():
        sub = p[:, idx]
        means[t] = np.nanmean(sub, axis=1) if sub.size else np.nan
        medians[t] = np.nanmedian(sub, axis=1)
        nexpr[t] = np.isfinite(sub).sum(axis=1)
    idx = pd.MultiIndex.from_arrays([block["transcript_id"].to_numpy(), genes], names=["transcript_id", "gene_id"])
    return (pd.DataFrame(means, index=idx), pd.DataFrame(medians, index=idx), pd.DataFrame(nexpr, index=idx))


def stream():
    with open_gct() as fh:
        fh.readline(), fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
    groups = tissue_columns(header)
    reader = pd.read_csv(TPM, sep="\t", skiprows=2, chunksize=CHUNK, dtype={"transcript_id": str, "gene_id": str})
    carry, parts = None, []
    for i, chunk in enumerate(reader):
        if carry is not None:
            chunk = pd.concat([carry, chunk])
        last = chunk["gene_id"].iloc[-1]
        carry = chunk[chunk["gene_id"] == last]
        done = chunk[chunk["gene_id"] != last]
        if len(done):
            parts.append(summarize_block(done, groups))
        if i % 10 == 0:
            print(f"chunk {i}: {sum(len(m) for m, _, _ in parts)} transcripts", flush=True)
    parts.append(summarize_block(carry, groups))
    return [pd.concat([p[k] for p in parts]) for k in range(3)]


def open_gct():
    return gzip.open(TPM, "rt")


def baselines(mean, nexpr):
    brain = [t for t in mean.columns if t.startswith("Brain")]
    nonbrain = [t for t in mean.columns if not t.startswith("Brain")]
    b = pd.DataFrame(index=mean.index)
    b["gtex_cortex"] = mean["Brain - Cortex"]
    b["gtex_ba9"] = mean["Brain - Frontal Cortex (BA9)"]
    b["generic_all"] = mean.mean(axis=1, skipna=True)
    b["generic_nonbrain"] = mean[nonbrain].mean(axis=1, skipna=True)
    b["n_tissues_expr"] = (nexpr > 0).sum(axis=1)
    b["n_cortex_expr"] = nexpr["Brain - Cortex"]
    b = b.reset_index()
    b["enst"] = b["transcript_id"].str.split(".").str[0]
    b["ensg"] = b["gene_id"].str.split(".").str[0]
    return b, len(brain), len(nonbrain)


def validate(b):
    """Compare our Brain - Cortex means with the supplied adult_brain_cortex.tsv.gz."""
    ref = pd.read_csv(ROOT / "data" / "raw" / "adult_brain_cortex.tsv.gz", sep="\t")
    j = ref.merge(b[["enst", "gtex_cortex"]], left_on="transcript", right_on="enst", how="left")
    ok = j["gtex_cortex"].notna()
    return {
        "reference_rows": len(ref), "matched": int(ok.sum()),
        "pearson_r": round(float(np.corrcoef(j.loc[ok, "mean_transcript_proportions"], j.loc[ok, "gtex_cortex"])[0, 1]), 5),
        "max_abs_diff": round(float((j.loc[ok, "mean_transcript_proportions"] - j.loc[ok, "gtex_cortex"]).abs().max()), 5),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mean, median, nexpr = stream()
    mean.to_csv(OUT / "tissue_mean_props.tsv.gz", sep="\t", float_format="%.5g")
    median.to_csv(OUT / "tissue_median_props.tsv.gz", sep="\t", float_format="%.5g")
    nexpr.to_csv(OUT / "tissue_n_expr.tsv.gz", sep="\t")
    b, n_brain, n_nonbrain = baselines(mean, nexpr)
    b.to_csv(OUT / "baselines.tsv.gz", sep="\t", index=False, float_format="%.5g")
    v = {"transcripts": len(b), "genes": int(b["ensg"].nunique()), "tissues": mean.shape[1],
         "brain_tissues": n_brain, "nonbrain_tissues": n_nonbrain, "cortex_validation": validate(b)}
    (OUT / "validation.json").write_text(json.dumps(v, indent=2))
    print(json.dumps(v, indent=2))


if __name__ == "__main__":
    main()

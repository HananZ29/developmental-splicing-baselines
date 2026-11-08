# External comparator tools: setup and validation (2026-11-05)

Pangolin, MMSplice and DeltaSplice are set up and validated (§9.3 item 1; §9.6 step 1).

**No comparator outcome on any test set has been computed.** The validation checks below use only:
- three COMPASS variants, for peak positions;
- 200 variants from the COMPASS **training** families, for the MMSplice sign check.

## Locations

- **Code:** `tools/src/<tool>`.
- **Isolated Python 3.10.11 environments:** `tools/envs/<tool>`.
- **Exact package versions:** `tools/envs/requirements_<tool>.lock.txt` (pip freeze).
- **Wrappers:** `publication_extension/scripts/score_{pangolin,deltasplice,mmsplice}.py`.
- **Genome:** all wrappers read the same hg19 FASTA as SpliceAI (`data/external/genome/hg19.fa`), with
  REF checked.

| Tool | Source | Version / commit | Licence | Key packages |
|---|---|---|---|---|
| Pangolin | github.com/tkzeng/Pangolin | 5cf94b8d; weights `final.{1-3}.{0,2,4,6}.3.v2` | GPL-3.0 | torch 2.2.2+cpu, numpy 1.26.4 |
| MMSplice | PyPI `mmsplice` (gagneurlab/MMSplice_MTSplice) | 2.4.0 (bundled models) | MIT | tensorflow 2.15.1, keras 2.15.0, kipoiseq 0.7.1, kipoi 0.8.6 |
| DeltaSplice | github.com/chaolinzhanglab/DeltaSplice | 71c2fd16; `DeltaSplice_models/model.ckpt-{0..4}` (SHA-256 prefixes cbed0a61, ad5f4e4f, fe5b2485, 1f3ccbe9, 1e63b7df) | MIT | torch 2.2.2+cpu, numpy 1.26.4 |

**Not downloaded from DeltaSplice:** `data/`, most of `experiments/`, and the alternative
`DeltaSplice_human` weights. They are not needed for scoring.

## Deviations from the authors' pipelines (none changes the model)

**Pangolin**
- Called through its Python model class instead of the `pangolin` CLI. The CLI needs pyvcf, which no
  longer installs, and a gffutils database.
- Unchanged from the CLI: sequence extraction, 4 tissues × 3 models, output channels [1, 4, 7, 10],
  and strand handling.
- No annotation masking.

**MMSplice**
- Installed without cyvcf2, which does not build on Windows. Only the sequence API is used
  (`predict_on_seq`, `predict_deltaLogitPsi`).
- Overhang is 100/100, the MMSplice default.

**DeltaSplice**
- Runs on CPU. `.cuda()` is made a no-op in the wrapper, without editing the repository.
- hg19 is read with pyfaidx instead of pyfasta.
- Reference-free tracks reproduce `pred_deltassu.py`'s no-exon mode. `exon_delta()` reproduces its
  exon mode, which takes measured reference usage (SSU) as input.

## Validation (setup only, not outcomes)

| Check | Result |
|---|---|
| Pangolin REF P(splice) peaks vs SpliceAI acceptor/donor peaks (3 COMPASS variants, − strand) | Identical genomic positions (103,534,669 acceptor; 103,534,489 / 103,535,497 donors) |
| DeltaSplice reference-free acceptor/donor peaks vs SpliceAI | Identical positions; peak SSU 0.85–0.97 |
| DeltaSplice exon mode: does the predicted reference usage follow the supplied measured usage? | Supplied 0.2 / 0.5 / 0.9 → predicted 0.23 / 0.52 / 0.91 (acceptor) and 0.18 / 0.49 / 0.90 (donor). Positions are placed correctly |
| COMPASS hg19 exon coordinates | 1-based inclusive: 151/157 exons are AG…GT, the rest have GC or non-canonical donors |
| MMSplice Δlogitψ vs measured Δlogit, 196 in-range variants from COMPASS **training** families | r = 0.61 (HEK), 0.62 (K562): correct sign and plausible magnitude |

## CPU throughput (this machine) and implications

| Tool | Speed | E13c subset (3,321 variants) | All COMPASS singles (35,723) |
|---|---|---|---|
| SpliceAI (for reference) | ~0.2 variants/s | 4.8 h (measured) | ~2 days |
| Pangolin | ~5 s/variant | ~4.5 h | ~2 days |
| MMSplice | ~1.5 s/variant (unbatched) | ~1.4 h | ~15 h; much less if batched |
| DeltaSplice | ~30 s/variant (35 kb input × 5 models × 2 sequences) | ~28 h | weeks. **GPU needed** |

## Leakage and comparability caveats for the benchmark

1. **MMSplice's Δlogitψ linear model was fitted on Vex-seq** (CAGI). MMSplice results on Vex-seq are
   therefore not independent and must be flagged or excluded (R3).
2. **DeltaSplice was evaluated, not trained, on Vex-seq and MFASS** (`experiments/eval_mut/run.sh`).
   - Training used a chromosome split; test chromosomes were 1, 3, 5, 7 and 9.
   - Including both sets in the paper's evaluation means the authors' design choices may have seen
     them. State this.
3. **The tools use the reference differently:**
   - SpliceAI and Pangolin are reference-free.
   - MMSplice gives Δlogitψ, to which the reference PSI is applied by logit additivity.
   - DeltaSplice's exon mode takes the measured reference usage as input.
   - A fair comparison feeds every reference-using method the **same replicate-split baseline**
     (replicate r for input, r′ for target), as in E13d.
4. **MMSplice scores only variants within 100 nt of the exon.** Variants farther out get Δ = 0 and
   are flagged `in_range = False`. The fraction must be reported per test set.

## Before any comparator is run on test data

Pre-register the comparator protocol, which is §9.6 step 2. It must specify:
- which outputs feed which prediction rule (raw scores; A-style propagation of Pangolin and DeltaSplice
  site ratios; MMSplice Δlogitψ + baseline; DeltaSplice exon mode with the replicate-split baseline);
- the identical eligible observations;
- the handling of out-of-range variants;
- the MMSplice/Vex-seq exclusion.

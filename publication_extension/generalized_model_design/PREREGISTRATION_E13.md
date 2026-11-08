# Pre-registration E13: diagnostics for a generalizable context-aware splice-effect model

**Registered 2026-11-04**, before any E13 outcome metric was computed.

**What had been inspected before registration:**
- COMPASS (GEO GSE307247) column definitions;
- row counts by `seq_type` and `source`;
- per-replicate read-depth quantiles;
- median wild-type (WT) PSI per cell line;
- the MFASS and Vex-seq column structure, already used in E5 and E11.

No predictor, correlation or contrast had been computed on COMPASS.

**Fixed:**
- PSI is clipped to [0.01, 0.99] before the logit (as COMPASS).
- Exon-family (event) cluster bootstrap, B = 1,000, seed 20261105.
- Pearson r with measured ΔPSI is the primary metric everywhere.
- **Effect-size threshold:** |ΔPSI| ≥ 0.1 defines a "relevant" effect (secondary AUPRC; E11 and
  literature convention). |Δlogit| ≥ 1 is reported as COMPASS's own disruption definition.

## E13a. Label audit (descriptive)

For every dataset on disk, report the number of distinct variants, genes, exons, biological
contexts, replicates, and the independent experimental units.
- **Measured variant-effect datasets:** MFASS, Vex-seq, COMPASS.
- **Association-based:** Wen, GTEx sQTL.
- **Expression only (no variant effects):** Leung, Patowary, GSE310131, GTEx RSEM.

**Overlap:** exon overlap between MFASS, Vex-seq and COMPASS, on hg38 exon coordinates (Vex-seq
lifted from hg19), at the identical-exon and same-gene levels.

**Reference independence:** whether reference PSI is measured in the same replicates as the target.

## E13b. Context decomposition in COMPASS (no sequence model)

**Data.** COMPASS single-SNV rows with a WT row in the same event, in cells
C = {HEK, HeLa, K562, MCF7, HMC3}.
- **QC per cell:** variant reads ≥ 20 and WT reads ≥ 20 in **both** rep1 and rep2. Read = included +
  excluded.
- Δ_{v,c,r} = logit PSI_var − logit PSI_wt in replicate r.

**Replicate split.** Input baseline from replicate r_in; target from replicate r_tg; both orders,
averaged.

**Predictors of the target-replicate ΔPSI in held-out cell c**, using only other cells' measured
effects:

| Predictor | Definition |
|---|---|
| T (transfer) | mean pooled ΔPSI_{v,c′} over c′ ≠ c |
| L (logit-additive, own baseline) | expit(logit p_c^{r_in} + Δ̄_v^{(−c)}) − p_c^{r_in}, where Δ̄^{(−c)} is the mean pooled Δ over c′ ≠ c |
| L_other (logit-additive, wrong baseline) | the same with p̄ = mean pooled WT PSI over c′ ≠ c |
| Ceiling | cor(ΔPSI_c^{r_in}, ΔPSI_c^{r_tg}), replicate reproducibility |

**Primary contrasts** (Pearson, pooled over cells and splits; exon bootstrap): L − T and L − L_other.

**Between-context decomposition (primary for candidate C).**
- For each variant and unordered cell pair (c, c′):
  - D_obs = ΔPSI_c − ΔPSI_c′ in the target replicate.
  - D_L = L_c − L_c′, with Δ̄ estimated from cells other than c and c′, and baselines from r_in.
  - Reliability: rel = cor(D^{rep1}, D^{rep2}).
- **Share of reliable between-context variance explained by baseline-only logit additivity:**
  S = cor(D_L, D_obs)² / rel. **Headroom** for learned interactions: H = 1 − S.
- Stratified by |WT PSI_c − WT PSI_c′| ≥ 0.1.

**Context-invariance of Δ.** Between-cell cor(Δ_c, Δ_c′) versus within-cell replicate
cor(Δ_c^{r1}, Δ_c^{r2}).

**Decision gate G-C.** Candidate C (a learned context interaction) is justified on these contexts only
if the bootstrap 95% lower bound of H exceeds 0.25.

## E13c/d. Sequence-feature models A vs B on unseen exons and independent datasets

**Subset.**
- 200 COMPASS exon families, drawn at random (seed 20261105) from events with a WT row and ≥ 5
  single SNVs passing QC in HEK and K562.
- Excluded: genes present in Vex-seq or MFASS.
- Variants are lifted hg38 → hg19 (UCSC chain) with REF checked against hg19, then scored with
  SpliceAI v1.3.1 exactly as before (genomic context; reporter context is not modelled, a stated
  limitation).
- Split by gene into 100 training and 100 test families (seed 20261105).

**Features.**
- r_acc and r_don: SpliceAI REF/ALT at the reporter exon's acceptor and donor, with ε = 0.01. The
  coordinate convention is fixed by REF-only modal peak offset, as in E11.
- Baseline p: WT PSI from r_in of the same cell. MFASS uses nat_v2_index_R1, with target from
  replicate R2.

**Models:**

| Model | Definition |
|---|---|
| A | p′ = pa/(pa + 1 − p), a = r_acc·r_don; prediction p′ − p; no fitted parameters |
| K | logit p′ = logit p + k·log a (1 parameter) |
| B | logit p′ = logit p + b + k_a·log r_acc + k_d·log r_don (3 parameters; least squares on ΔPSI; ridge λ = 1e-3) |
| R | raw SpliceAI (Δp_acc + Δp_don), for correlation only |

K and B are fitted on the COMPASS training families only, pooled over HEK, HeLa, K562, MCF7 and HMC3.

**Test sets:**
1. COMPASS test families (unseen exons and genes; all 5 cells);
2. Vex-seq (K562, HepG2; replicate split as in E11b);
3. MFASS (HEK293T; replicate split).

**Primary contrast:** B − A, Pearson r, on each test set separately (3 tests, Bonferroni
α = 0.0167). **Secondary:** MAE; AUPRC for |ΔPSI| ≥ 0.1; K − A.

**Gate G-B.** Model B is "validated" only if B − A > 0 with a CI excluding 0 on at least two of the
three test sets, **and** B is not worse than R on any of them.

## Not done here

- No large neural model is trained.
- No model is tuned on any test set.
- No developmental-stage claim is made from these data: all contexts are non-neural or glial cancer
  or immortalized lines.

## Deviations log (E13)

**D-E13-1** (2026-11-04, before any E13c/d outcome):
- **What changed:**
  - Only 155 COMPASS families met the eligibility rule after excluding loci within 100 kb of Vex-seq or
    MFASS exons. 314 had ≥ 5 QC singles; 159 were excluded as nearby. All 155 are used, instead of a
    sample of 200.
  - To keep SpliceAI scoring to about 2 h on CPU, at most 40 randomly chosen variants per family are
    scored (seed 20261105): 3,321 variants.
  - The train/test split is by gene: 78 / 77 families.
- **Why acceptable:** decided from counts only, with no predictor or outcome computed.

**D-E13-2** (secondary analysis added, descriptive). `E13b2_headroom_scale.py` expresses the E13b
headroom on an absolute scale (fraction of total ΔPSI variance). It was added after E13b, before E13d.
It is not part of any gate.

**Code-path test.** E13d was run once in a test mode: COMPASS skipped, parameters fixed at A-equivalent
values, bootstrap of 20, outputs deleted. Its only purpose was to catch errors. Fixed parameters make
A = K = B1 by construction, so it reveals nothing about the contrasts.

**D-E13-3** (2026-11-04, recorded after the E13d run; a clarification, not a change). As implemented
and as pre-registered for logits, baselines p are clipped to [0.01, 0.99] for **all** models,
including the parameter-free A. This affects 31% of COMPASS test rows, 2% of Vex-seq rows and 22% of
MFASS rows. A post-hoc, descriptive scan of the best k on each test set is reported in
`DIAGNOSTIC_RESULTS_E13.md`. It is labelled as such and not used for any gate or claim.

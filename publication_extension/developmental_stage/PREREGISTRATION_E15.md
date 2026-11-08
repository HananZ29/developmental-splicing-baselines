# PREREGISTRATION E15 — FROZEN

**Status.** FROZEN 2026-11-06, approved by the user ("Approve and freeze").
- This is draft v3 unchanged, except for this status block.
- Decision D3-B applies: MEI = 0.0228.
- SHA-256 of this file and of all listed code and inputs are in `PREREG_E15.sha256`.
- No QTL association effect had been read at freeze, apart from the single row disclosed in
  `DATA_PROVENANCE_E15.md`.
- It is frozen only after the user's explicit approval. Freezing consists of:
  - recording the SHA-256 of this file and of every analysis script in `PREREG_E15.sha256`;
  - renaming this file to `PREREGISTRATION_E15.md`.
- **No QTL outcome has been read.** The one exception is the single row disclosed in
  `DATA_PROVENANCE_E15.md`.
- Fields marked **[FILL: …]** are filled from outcome-blind runs (baseline diagnostics, SpliceAI,
  Gate 0) before freezing. Filling them changes no rule below.

E14 is a separate, frozen experiment. Nothing here changes it.

---

## 1. Research question and E15's place in it

**Original research question (unchanged):**
> Does incorporating a fetal/developmental brain isoform baseline improve prediction of isoform-level
> consequences of splice-altering variants compared with adult or generic baselines?

**What E15 tests.** E15 is a **junction-level validation component** of that question. It asks:

> Holding the variant set and the splice-site predictor (SpliceAI) fixed, and changing **only the
> baseline**, does a fetal-brain baseline produce predicted changes in local intron usage that agree
> better with prenatal genetic association effects on intron usage than an adult-brain, a generic
> (GTEx all-tissue) or a uniform baseline?

**What E15 does not test:**
- It does not evaluate full-length isoform proportions. No public prenatal ground truth for
  full-length isoform proportions exists (`E15_PREFREEZE_VERIFICATION.md`, point 6).
- Its truth is a set of **association estimates** (fine-mapped sQTL z-scores), not measured causal
  variant effects.

**Permitted conclusion if positive:**
> "A fetal baseline improved junction-level agreement with prenatal sQTL association effects."

**Not permitted under any result:**
- "fetal baselines predict full-length isoform consequences better";
- "predicted effects were validated";
- any causal wording.

---

## 2. Data (all verified; checksums in `DATA_PROVENANCE_E15.md`)

| Role | Source | Build |
|---|---|---|
| F1/A1 baseline reads | BrainSpan neocortex MRF alignments: 362 samples, 42 donors. Junctions extracted by `E15_mrf_junction_counts.py` (blocks gap ≥ 20 nt, unstranded, CRC-checked) | hg19 |
| G baseline reads | GTEx v8 STAR junction GCT: 17,382 samples, 54 tissues. Annotated (GENCODE v26) junctions only. Lifted to hg19 by `E15_gtex_generic_counts.py` | hg38 → hg19 |
| Intron universe | Wen et al. `all.introns.tested.tsv` (syn51183859). Wen `chr:a:b` = intron first base a+1, last base b−1 (verified) | hg19 |
| Primary truth | Wen SuSiE fine-mapping, `mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz` (syn51183885) | hg19 |
| Secondary truth | Wen QTLtools conditional top variants (syn51183855) | hg19 |
| Crossover truth (adult) | GTEx v8 `Brain_Cortex.v8.sqtl_signifpairs.txt.gz` | hg38 → hg19 |
| Annotation | GENCODE v31lift37 (patterns), v26 (GTEx presence) | hg19 |
| Predictor | SpliceAI v1.3.1 weights. Repository runner `src/spliceai_runner.py`, ±500 nt tracks, unchanged from E14 | hg19 |

---

## 3. Baselines

All baselines are computed on the **same intron universe** (Wen tested introns) by the same code
(`E15_baseline_diagnostics.py`, then `E15_predict.py`).

| ID | Definition | Units |
|---|---|---|
| **F1** fetal (primary) | BrainSpan prenatal donors ≤ 26 pcw **not in the Wen cohort** (decision D1, approved 2026-11-06). Excludes HSB107, HSB148, HSB150, HSB154 and HSB194 (also in Wen's cohort). Units: **HSB112, HSB113, HSB153, HSB103, HSB114, HSB149, HSB178, HSB96, HSB97, HSB98, HSB92, HSB159** | 12 donors, 8–22 pcw |
| **A1** adult (primary comparator) | BrainSpan donors aged ≥ 18 years: **HSB105, HSB123, HSB126, HSB127, HSB130, HSB135, HSB136, HSB145**. Same pipeline and regions as F1 | 8 donors |
| **G** generic | GTEx v8, **all 54 tissues including brain** (as planned) | 54 tissues |
| **U** uniform | Equal weight on every annotated local splicing pattern in a cluster | — |

**GTEx processing rules.**
- A junction listed once per overlapping gene is counted once. This skipped 5,695 repeated rows.
- Liftover keeps a junction only if both ends map uniquely to the same chromosome and strand with
  intron length unchanged.
- hg19 introns receiving two or more distinct hg38 junctions (24 introns; 25 collisions) are dropped from G, so rule E4 excludes
  their clusters.

**Unit-level aggregation (R3).**
- For BrainSpan, a donor's reads are summed over all of that donor's neocortex samples.
- For GTEx, a tissue's reads are summed over its samples.
- For each unit with ≥ 5 annotated reads in a cluster, the within-cluster proportions over annotated
  introns are computed.
- The **group usage** is the unweighted mean of unit proportions.

**Depth equalisation (R4).**
- Per cluster, each unit's counts are binomially thinned by T / pooled_group.
  - T is the smallest pooled annotated read total among F1, A1 and G for that cluster.
  - pooled_group is the group's pooled total.
- This uses 20 fixed seeds (20261106 + k, k = 0…19).
- Predictions are computed per seed and averaged over seeds before evaluation.
- U is not thinned.

**Pattern weights.**
- Local splicing patterns: for each eligible cluster, the distinct sets of the cluster's annotated
  introns used by one GENCODE v31lift37 transcript.
- Pattern weights w_B are fitted to group usage by EM. The model is
  P(read on intron i) ∝ Σ_p w_p·1[i ∈ p].
- The EM recovers known weights exactly on test cases.

---

## 4. Eligible population

### 4.1 Clusters

Each rule is applied identically for every baseline. Clusters are Wen LeafCutter clusters.

| Rule | Criterion |
|---|---|
| E1 | ≥ 2 annotated introns (v31lift37), and the intron × pattern incidence matrix has full column rank, so pattern weights are identifiable from junction counts |
| E2 | Applies in each of F1, A1 and G: pooled annotated reads ≥ 20, from ≥ 3 contributing units (≥ 5 reads each) |
| E3 | Applies in each of F1, A1 and G: ≥ 95% of the cluster's reads are on annotated introns |
| E4 | Every annotated intron of the cluster is present in the GTEx junction set. GTEx lists only v26 junctions, so an absent intron would be a structural zero, not a measurement |
| — | Annotated transcripts using the cluster's introns are all on one strand (2 clusters fail) |

**Eligibility waterfall** (`baseline_diagnostics.json`, all 362 samples):

| Stage | Clusters |
|---|---|
| Wen clusters | 50,311 |
| ≥ 2 annotated introns | 32,777 |
| E1 (also full rank) | 24,573 |
| E1 + E4 | 22,354 |
| E1 + E4 + E2 (coverage in F1, A1 and G) | 20,792 |
| **E1–E4: eligible** | **18,496** (88,471 introns, of which 58,799 annotated) |

### 4.2 Units of analysis

A **unit** is one SuSiE credible set (CS) for one phenotype intron i.

**Inclusion.** Intron i must be an annotated intron of an eligible cluster. The CS must contain ≥ 1
**candidate**. A candidate is a CS member that meets all of the following:
- it is a biallelic SNV;
- its REF allele matches hg19;
- it lies within ±50 nt of an annotated splice site of the cluster.

Only Wen columns `phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index` are used to build
units. PIP, z, posterior and Bayes-factor columns are not used.

Candidates before the coverage rules: 1,917 CS and 1,373 unique SNVs (`candidates_stats.json`).
**After all rules** (`predictions_stats.json`, `prefreeze_report.json`):

| Measure | Count |
|---|---|
| Units (variant–intron pairs) | **1,449**, all with F1, A1, G and U predictions |
| Distinct representative variants (rule P) | 858 |
| Phenotype introns | 1,412 |
| Clusters | 800 |
| Genes | 845 |
| Signal groups | 847 |
| Units with DS_max ≥ 0.10 | 589 |
| Units where rules P and D select the same variant | 1,362 (94%) |

Because rules P and D agree for 94% of units, rule D is a weak test of the selection bias. Stratum
analysis S is the more informative check.

### 4.3 Representative variant per unit

**Primary rule P: maximum predicted splice disruption.**
- Choose the candidate with the largest SpliceAI DS_max = max(DS_AG, DS_AL, DS_DG, DS_DL).
- Ties are broken by lower position, then alphabetical ALT.

*Justification:*
- A CS represents one association signal, and its causal member is unknown.
- E15 asks how well a baseline converts a **predicted** splice disruption into intron-usage change.
  The member with the strongest predicted disruption is the one for which the baseline is actually
  exercised. That is how the predictor would be used in practice.
- The rule uses only the predictor's REF/ALT site scores. These do not depend on any baseline, so the
  same variant is used for all four baselines and selection cannot favour one baseline.

*Known bias:*
- Rule P enriches for large predicted effects, which can inflate every ρ.
- It can also change Δρ, because the baseline matters more or less for strong site effects.

This is why the sensitivity rules below are pre-specified.

**Sensitivity rule D: deterministic genomic position.**
- Choose the candidate closest to any annotated splice site of the cluster.
- Ties are broken by lower position, then alphabetical ALT.
- Coordinates only; no predictor input.

**Sensitivity rule S: strength strata.**
- Run the primary analysis separately in two strata: units whose rule-P variant has DS_max ≥ 0.10,
  and units with DS_max < 0.10.

### 4.4 Overlap and non-independence

**(a) Same signal on several introns.** SuSiE is run per intron, so one genetic signal often appears
as CSs for several introns of the same cluster. For example, the skipping intron and both inclusion
introns of a cassette exon.
- *Primary:* all such units are kept. Dependence is handled by (c).
- *Sensitivity:* one unit per **signal group**. A signal group is the set of CSs in the same cluster
  that share ≥ 1 member variant. The kept unit is the one whose phenotype intron has the lowest
  genomic start.

**(b) Several independent CSs for one intron.** Each CS is its own unit. Dependence is handled by (c).

**(c) Resampling components.** Units are merged (union–find) into **components** if any of the
following holds:
- they share a gene (Wen `leafcutter_clusters_to_genes.txt`);
- they share a LeafCutter cluster;
- they share the representative variant;
- their representative variants lie within 100 kb of each other.

The last criterion is a distance proxy for LD, because no LD reference is used. The component is the
bootstrap resampling unit. Count: **675 components** under rule P (676 under rule D).

**Crossover arm:** 428 GTEx units in 411 clusters, forming 383 components.

---

## 5. Prediction (identical for all baselines; `E15_predict.py`)

For unit (i, v) and baseline B:

1. **Site ratios.** At each annotated splice site s of the cluster within ±500 nt of v, compute
   r_s = (p_ALT(s) + ε) / (p_REF(s) + ε), with ε = 0.01 (as `src/propagation.py`). Sites farther away
   get r_s = 1.
   - Donor and acceptor positions follow `src/splice_graph.py`.
   - On the + strand, the donor is at intron_start − 1 and the acceptor at intron_end + 1.
   - On the − strand, the donor is at intron_end + 1 and the acceptor at intron_start − 1.
2. **Pattern factor.** For each pattern p, f_p = Π r_s over the donor and acceptor sites of p's
   introns.
3. **Reweighting.** w'_p ∝ w_{B,p}·f_p.
4. **Intron usage.** u_i = Σ_{p ∋ i} w_p / Σ_p |p| w_p. The prediction is Δu_i = u'_i − u_i, computed
   REF→ALT.

The prediction is averaged over the 20 thinning seeds (U is not thinned).

**Assumptions** (stated, not tested here):
1. Sites act independently and multiplicatively.
2. Only annotated patterns exist.
3. No NMD or transcriptional effect.
4. The baseline is the relevant stage's.

---

## 6. Truth

- **Primary truth:** the SuSiE `z` of the unit's representative variant for the unit's phenotype
  intron. This is the univariate z of ALT dosage against the rank-normalised, ComBat-corrected
  excision ratio.
- Sign: positive means ALT is associated with higher usage of intron i. This is the same orientation
  as Δu_i, so no sign flip is needed.
- If the representative variant has no z row for that phenotype (not expected, because the units are
  built from the same file), the unit is dropped and counted.
- **Secondary truth:** `bwd_slope` (column 18) from the conditional top-variants file. It is used for
  units whose representative variant matches a conditional top variant by chromosome and position.
  Reported only, with no inferential claim.

---

## 7. Primary endpoint

ρ_B = Spearman correlation, over all units, between the predicted Δu_i under baseline B and the truth
z.

**Primary contrast:** Δρ_A = ρ_F1 − ρ_A1.

---

## 8. Statistical comparison and uncertainty (`E15_stats.py`)

**Family of three contrasts** (all paired: same units, same truth):
- Δρ_A = ρ_F1 − ρ_A1 (primary);
- Δρ_G = ρ_F1 − ρ_G;
- Δρ_U = ρ_F1 − ρ_U.

**Point estimate.** Pearson correlation of full-sample average ranks.

**Uncertainty: paired component bootstrap.**
- B = 10,000 replicates, seed 20261106.
- Each replicate draws all components with replacement. A unit's weight is its component's draw
  count.
- ρ is the weighted correlation of the full-sample ranks. Ranks are not recomputed within a
  replicate. The identical estimator is used in Gate 0.

**p-value.** Bootstrap two-sided p = min(1, 2·min(P*(Δ ≤ 0), P*(Δ ≥ 0))).

**Multiple testing.** Holm step-down across the three contrasts at family-wise α = 0.05. The reported
CIs are percentile intervals at each contrast's Holm-adjusted level.

---

## 9. Criteria for declaring improvement

For each contrast F1 − B:

| Outcome | Criterion |
|---|---|
| **Improvement** | Holm-adjusted p ≤ 0.05 **and** Δρ > 0 |
| **Practically meaningful improvement** | Improvement **and** point Δρ ≥ MEI (minimum effect of interest; **decision D3**, section 9a) |
| **Worse** | Holm-adjusted p ≤ 0.05 and Δρ < 0 |
| **No meaningful difference** | The 95% CI lies entirely within (−MEI, +MEI) |
| **Inconclusive** | Anything else |

### 9a. Decision D3 (user): the minimum effect of interest. The test itself is unchanged.

> **DECIDED 2026-11-06: D3-B.** MEI = **0.0228** = 0.10·√(2(1 − r)), with r = 0.97404 the Pearson
> correlation of rank-normal F1 and A1 predictions. The value is fixed now and is not recomputed
> later. Recorded in `results/E15_devstage/gate0_mei.json` and `E15_analysis.py`
> (`MIN_EFFECT = 0.0228`).
>
> MEI = 0.0228 is used for:
> - "Practically meaningful improvement": point Δρ ≥ 0.0228;
> - "No meaningful difference": 95% CI within (−0.0228, +0.0228).
>
> Interpretation: a Δρ of 0.0228 corresponds to truth correlating 0.10 with the fetal-minus-adult
> prediction difference.
>
> The options considered are kept below for the record.

**Why this is open.**
- The MEI of 0.05 was set (AUDIT D4) before the predictions existed.
- The predictions show that fetal and adult baselines produce nearly identical variant rankings:
  Spearman 0.981, or 0.974 Pearson on the rank-normal scale.
- So Δρ is on a compressed scale. A Δρ of 0.05 requires truth to correlate about 0.22 with the
  fetal-minus-adult prediction difference, i.e. 0.05 / √(2(1 − r)) with r = 0.974. That is a large
  effect.
- At this compression, an equivalence margin of ±0.05 is so wide that "No meaningful difference"
  would almost always be declared. It would exclude only large effects.

This was found from the predictions alone, before any truth was read.

| Option | MEI and margin | Gate 0 power at the MEI | Would establish | Would not establish |
|---|---|---|---|---|
| **D3-A** (as drafted) | 0.05 | 1.00 | A "practically meaningful" fetal advantage only if truth aligns strongly (ρ ≈ 0.22) with the fetal–adult prediction difference | It cannot rule out small or moderate advantages. A "No meaningful difference" verdict excludes only effects ≥ 0.05, which is close to uninformative here |
| **D3-B** (recommended): standardized MEI | 0.10·√(2(1 − r)) = **0.023** (r = 0.974, fixed now from predictions). Truth correlates 0.10 with the fetal–adult prediction difference | ≥ 0.94 (power was 0.94 at Δρ = 0.02) | Whether the fetal-specific component of the prediction carries a modest, detectable alignment with prenatal association effects. A "No meaningful difference" verdict then excludes alignments ≥ 0.10 | Any absolute accuracy gain beyond ≈ 0.02 in ρ. Nothing about full-length isoforms, causality, or clusters outside the annotated, identifiable set |
| D3-C: restrict to units with \|pred_F1 − pred_A1\| > 0.01 | 0.05 | **0.57 → Gate 0 fails** (329 units, 169 components) | — | Not recommended. It is underpowered and is kept only as descriptive analysis 8 |

**Fixed under every option:**
- the statistic (Δρ);
- the paired component bootstrap;
- Holm across F1 − A1, F1 − G and F1 − U;
- the Improvement/Worse rule;
- the population and the robustness rule.

D3 changes only the MEI used for "practically meaningful", for the equivalence margin and for the
Gate 0 pass criterion. Under D3-B, Gate 0 passes if power is ≥ 0.80 at Δρ = MEI.

**E15-level conclusion:**
- "The fetal baseline improved junction-level agreement over adult **and** generic baselines" only if
  F1 − A1 **and** F1 − G are both declared Improvement.
- If only F1 − A1 is declared, the statement is limited to the adult comparison.
- F1 − U is reported but is not required for the conclusion. It checks that any baseline information
  helps.
- A primary-rule result that reverses sign under rule D, or under the signal-group sensitivity
  analysis, is reported as **not robust**. The conclusion is then downgraded to "Inconclusive".

---

## 10. Gate 0 (before any truth is read; `E15_gate0_power.py`)

**What Gate 0 is.**
- It gives the operating characteristics of the exact test in sections 8–9, using the real units, the
  real F1/A1/G/U predictions and the real component structure, with **simulated** truth.
- The "fetal-favourable" scenarios build the fetal advantage into the simulated truth. Their power is
  the chance of detecting a fetal advantage **if one exists**.
- **Gate 0 is not evidence that fetal predictions are correct, and it is not evidence about the
  direction of the real result.**

**Model.** Truth y = a·x_F1 + b·x_A1 + c·x_G + noise, on the rank-normal scale of each prediction.
R = 500 simulated datasets per scenario, B = 500 bootstrap replicates each, seed 20261106.

| Scenario | Truth | Reports |
|---|---|---|
| N0 pure null | a = b = c = 0 | Type I error (false "F1 improves" declarations) for each contrast |
| N1 equal-signal null | ρ_F1 = ρ_A1 = 0.20 | Type I error for F1 − A1 |
| F fetal-favourable | ρ_F1 ∈ {0.10, 0.20, 0.30} × Δρ ∈ {0.05, 0.10}, plus Δρ = 0.02 at ρ = 0.20, plus standardized Δρ = s·(1 − r), s ∈ {0.25, 0.5, 0.75} | Power |
| F-t3, F-het | ρ = 0.20, Δρ = 0.05 and Δρ at s = 0.5, with t(3) noise or heteroscedastic noise | Robustness of power and error control |
| A adult-favourable | Δρ = −0.05 and Δρ at s = −0.5, at ρ_A1 = 0.20 | Rate of falsely declaring F1 better; power to detect "Worse" |
| G generic-favourable | ρ_G = 0.20 | F1 − G and F1 − A1 decisions |
| M mixture | 0.7·x_F1 + 0.3·x_A1 | Partial fetal advantage |

**Pass criterion (D3-B).** Both of the following:
- power ≥ 0.80 for F1 − A1 at ρ_F1 = 0.20, Δρ = MEI = 0.0228;
- Type I error for F1 − A1 ≤ 0.05 in N0 and in N1 (Monte-Carlo SE reported).

The original criterion was the same with Δρ = 0.05. Its result is reported too.

**Result at the MEI** (`gate0_mei.json`; R = 500, B = 500, seed 20261107):

| Measure | Result |
|---|---|
| Simulated Δρ | 0.0226 |
| Power | **0.984** (MC SE 0.006) |
| Maximum Type I error | 0.018 |
| Verdict | **PASS** |

`E15_analysis.py` refuses to read truth unless this record shows PASS with MEI = 0.0228 **and** the
Type I bound in `gate0_power.json` holds.

**If Gate 0 fails:**
- E15 stops **before** any truth is read.
- It is reported as "underpowered with public data".
- No improvement claim is possible, and the null is not interpreted as evidence of no difference.

**Pre-registration check.** The Gate 0 code was validated on synthetic predictions (not E15 data):
n = 2,000, prediction correlation F1–A1 = 0.78, R = 500, B = 500.

| Measure | Result |
|---|---|
| Type I error, F1 − A1 | N0: 0.010; N1: 0.032 |
| Power, F1 − A1, ρ = 0.20, Δρ = 0.05 | 0.90 |
| Power, F1 − A1, ρ = 0.20, Δρ = 0.02 | 0.25 |
| Adult-favourable truth: false "F1 better" rate | 0.000 |
| Adult-favourable truth: "Worse" detected | 0.90 |
| Simulated Δρ against a target of 0.05 | 0.047 (slightly conservative) |

**Interpretive caveat seen in the synthetic run.** When truth was **generic-like**, F1 − A1 was still
declared "Improvement" in 30% of datasets. This is because the synthetic F1 predictions resembled G
more than A1 did. A fetal win over adult alone can therefore arise without any fetal-specific
signal. This is why the E15-level conclusion in section 9 requires both F1 − A1 and F1 − G. The same
check is repeated on the real predictions.

**Calibration.** Each scenario's coefficients are tuned iteratively so that the simulated mean
Spearman ρ_F1 and Δρ hit their targets. A scenario whose target cannot be reached is reported as
unattainable. Standardized scenarios set Δρ = s·(1 − r) with s ∈ {0.25, 0.5, 0.75}.

**Results on the real E15 predictions** (1,449 units, 675 components; `gate0_power.json`):

| Scenario (ρ_F1 = 0.20 unless stated) | Simulated Δρ | F1 > A1 declared | F1 < A1 declared | F1 > G declared |
|---|---|---|---|---|
| N0 pure null | −0.0005 | **0.018** (SE 0.006) | 0.014 | 0.01 |
| N1 equal-signal null | −0.0004 | **0.012** (SE 0.005) | 0.020 | 0.09 |
| Fetal, Δρ = 0.0065 (s = 0.25) | 0.0063 | 0.16 | 0.00 | 0.16 |
| Fetal, Δρ = 0.013 (s = 0.5) | 0.0130 | 0.64 | 0.00 | 0.31 |
| Fetal, Δρ = 0.0195 (s = 0.75) | 0.0191 | 0.94 | 0.00 | 0.42 |
| Fetal, Δρ = 0.02 | 0.0199 | 0.94 | 0.00 | 0.44 |
| Fetal, Δρ = 0.05 | 0.0492 | **1.00** | 0.00 | 0.90 |
| Fetal, Δρ = 0.05, ρ_F1 = 0.10 / 0.30 | 0.049 | 1.00 / 1.00 | 0.00 | 0.77 / 0.98 |
| Fetal, Δρ = 0.10 (all ρ levels) | 0.099 | 1.00 | 0.00 | 1.00 |
| Fetal, Δρ = 0.05, t(3) noise / heteroscedastic noise | 0.061 / 0.047 | 1.00 / 1.00 | 0.00 | 0.99 / 0.89 |
| Fetal, Δρ = 0.013, t(3) noise / heteroscedastic noise | 0.017 / 0.012 | 0.88 / 0.54 | 0.00 | 0.53 / 0.24 |
| Adult, Δρ = −0.05 | −0.0495 | 0.00 | **1.00** | 0.00 |
| Adult, Δρ = −0.013 | −0.0130 | 0.00 | 0.57 | 0.02 |
| Generic-favourable truth (ρ_G = 0.20) | +0.0016 | 0.03 | 0.002 | 0.00 (F1 < G: 0.11) |
| Mixture 0.7 F1 + 0.3 A1 | +0.0019 | 0.04 | 0.004 | 0.11 |

**Verdicts:**
- **Original criterion (Δρ 0.05): PASS.** Power is 1.00, and the maximum Type I error is 0.018.
- **Adopted D3-B criterion (MEI 0.0228): PASS.** Power is 0.984, and Type I error is unchanged.
- Power drops below 0.8 for Δρ < ≈ 0.016. Smaller effects, such as truth = fetal prediction + noise at
  ρ = 0.2, which implies Δρ ≈ 0.005, are not reliably detectable.
- The generic-favourable-truth issue seen with synthetic data does not recur here: F1 > A1 was falsely
  declared in only 3% of datasets.
- The candidate "difference-correlation" statistic is valid on these data (Type I error 0.024). It is
  a rescaling of Δρ, so it is not adopted as a separate endpoint.

**Subset alternative (D3-C)** (329 units with |pred_F1 − pred_A1| > 0.01, 169 components):

| Measure | Result |
|---|---|
| Power at Δρ = 0.05 | 0.57 → **FAIL** |
| Type I error | 0.014 |
| Power at Δρ = 0.0375 | 0.26 |
| Power at Δρ = 0.056 | 0.65 |

Three subset scenarios could not be calibrated and are reported as unattainable.

---

## 11. Pre-specified secondary and sensitivity analyses

These are reported in full and do not change the primary decision, except under the robustness rule
in section 9.

1. Selection rule D (deterministic position).
2. Strength strata S (DS_max ≥ 0.10 vs < 0.10).
3. One unit per signal group.
4. Baseline variants:
   - F1_all (adds HSB155, 35 pcw);
   - G_nonbrain;
   - F1/A1 restricted to the 11 standard neocortical regions;
   - F1/A1 excluding samples with spliced-read fraction < 0.10;
   - no thinning.
5. Secondary truth: `bwd_slope` from the conditional file.
6. ε ∈ {0.001, 0.1}.
7. **Stage-swap crossover (secondary, descriptive).**
   - Truth: GTEx adult cortex sQTL slopes. Significant pairs only, because GTEx all-pairs data are
     not public. Variants and introns are lifted to hg19, and alleles are aligned to hg19 REF/ALT.
   - Restricted to eligible clusters, with variants selected by rule P from SNVs within ±50 nt.
   - Report [ρ_F1 − ρ_A1] for prenatal truth minus the same for adult truth. A component bootstrap
     gives the CI.
   - The variant sets differ between arms, so this is not a matched-pair test. No confirmatory claim.
8. Δρ_A by F1 vs A1 baseline difference (per-cluster max |Δu_baseline| tertiles): does the advantage
   concentrate where the baselines differ?

---

## 12. Outcome-blind diagnostics required before freezing (user requests 1–4)

From `baseline_diagnostics.json` and `prefreeze_report.json`, all outcome-blind:

**Baselines:**

| Baseline | Units | Reads on Wen introns | Clusters with coverage |
|---|---|---|---|
| F1 | 12 donors | 175.2 M | 46,379 |
| F1_all | 13 donors | 175.6 M | 46,388 |
| A1 | 8 donors | 99.1 M | 44,714 |
| G | 54 tissues | 81.8 B | 48,417 |
| G_nonbrain | 41 tissues | 71.7 B | 48,395 |

Sample QC: all 362 samples passed CRC, and all have ≥ 10 M reads. 56 have spliced fraction < 0.10;
they are excluded only in the sensitivity analysis.

**F1 stability (18,496 eligible clusters):**

| Measure | Result |
|---|---|
| Leave-one-donor-out, worst-case per-cluster max \|Δu\| | median 0.010, p95 0.038 |
| Share of clusters with worst-case shift > 0.05 | 2.5% |
| Share of clusters shifting > 0.05 when any single donor is removed | 0.2–0.5% |
| Donor-bootstrap SD of usage, F1 | median 0.010 (p95 0.043) |
| Donor-bootstrap SD of usage, A1 | median 0.013 (p95 0.054) |
| Leave-one-donor-out pattern-weight change | median 0.008, p95 0.059 |

**F1 vs A1 baselines.**
- The median per-cluster max |Δu| is 0.034, and p95 is 0.22.
- 27% of eligible clusters differ by > 0.05 and by > 2 bootstrap SE.
- The developmental baselines therefore differ substantially.
- The *predictions* nonetheless agree closely:

| Prediction comparison | Spearman |
|---|---|
| F1–A1 | 0.981 |
| F1–G | 0.958 |
| F1–U | 0.923 |

| F1 vs A1 prediction difference | Value |
|---|---|
| Median \|Δ\| | 0.0019 |
| Share with \|Δ\| > 0.01 | 23% |
| Share with \|Δ\| > 0.05 | 5.4% |
| Opposite sign | 1.0% |

The variant's site effect dominates the predicted Δu, and the baseline mostly rescales it.

**Effect of excluding samples with spliced fraction < 0.10:**

| Comparison | Spearman |
|---|---|
| F1 vs F1_sf predictions | 0.9998 |
| A1 vs A1_sf predictions | 0.994 |
| F1 − A1 difference preserved | 0.93 |

Other sensitivity baselines barely change the predictions:

| Comparison | Spearman |
|---|---|
| F1 vs F1_all | 0.9995 |
| F1 vs F1_std | 0.9996 |
| F1 vs no thinning | 0.9995 |
| G vs G_nonbrain | 0.998 |

**Exclusion bias** (base: 20,978 clusters that pass E1 and F1's own criteria):

| Exclusion | Share |
|---|---|
| Lost to E4 | 6.7% |
| Lost to A1 coverage or annotated share | 5.8% |
| Lost to G coverage | 0.6% |

| Measure | Result |
|---|---|
| Retention by F1 expression quartile, Q1→Q4 | 77%, 88%, 93%, 95% (low-expression clusters lost more often) |
| Fetal-only-expressed clusters lost | 666 (3.2%) |
| Retention by F1 − A1 usage difference, \|Δu\| < 0.05 → > 0.2 | 94% → 91% (ratio 0.96) |

Developmentally variable clusters are **not** excluded disproportionately: the ratio is far above the
0.5 threshold.

**Main limitation, upstream of all of the above.** 49% of F1-covered clusters fail E1. The propagation
model needs identifiable annotated patterns. A further 10% fail F1's annotated-share rule. E15 cannot
speak to these clusters.

**Pre-declared consequence.** If exclusion removes developmentally variable clusters
disproportionately (a retention ratio < 0.5 between the top and bottom F1 − A1 difference strata),
this is stated as a limitation in every E15 conclusion. The primary analysis is not altered after
seeing it.

---

## 13. Deviations

Any change after freezing is logged in `E15_DEVIATIONS.md` with date and reason, and is labelled
exploratory. Analysis code is hash-locked at freeze.

The sequence is: run Gate 0, fill the fields, obtain user approval, freeze, and only then read truth
columns.

**Files to be hash-locked at freeze** (SHA-256 in `PREREG_E15.sha256`):

*This document:*
- `PREREGISTRATION_E15.md`

*Scripts:*
- `E15_mrf_junction_counts.py`
- `E15_gtex_generic_counts.py`
- `E15_baseline_diagnostics.py`
- `E15_candidates.py`
- `E15_gtex_candidates.py`
- `E15_predict.py`
- `E15_stats.py`
- `E15_gate0_power.py`
- `E15_prefreeze_report.py`
- `E15_analysis.py`
- `score_spliceai.py`
- `src/spliceai_runner.py`
- `src/propagation.py`
- `src/splice_graph.py`

*Outcome-blind inputs to the analysis:*
- `predictions.tsv`
- `predictions_gtex_crossover.tsv`
- `gate0_power.json`
- `gate0_mei.json`
- `cluster_eligibility.tsv.gz`

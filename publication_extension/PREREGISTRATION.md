# Pre-registration of the extension experiments

**Written 2026-11-03T23:23Z**, before any new prediction or evaluation in this extension was run.
Nothing below may be changed after results are seen. Any deviation is logged at the end with its
reason.

The v2 experiment (`results/experiment_v2/`) and the manuscript are frozen. All new outputs go to
`publication_extension/`.

**Fixed for every experiment:**
- SpliceAI v1.3.1, five-model mean, hg19, 10 kb context.
- Propagation as in v2: ε = 0.01, α = 1e-3, sites within ±500 nt, connected-component intron
  clusters.
- 1-Mb block bootstrap, B = 2000.

## E1. Independent replication of v2

- **Method:** rerun `experiments/08_propagate_evaluate.py` (v2 settings) into a fresh directory and
  compare every estimate with `results/experiment_v2/metrics_with_ci.tsv`.
- **Pass criterion:** identical point estimates, because the run is deterministic.
- **Also:** run the independent audit (`13`) on v2.

## E2. Root cause: can the baseline change the predicted sign?

- **Two-isoform algebra:**
  `p' − p = p(1−p)(a−b) / (pa + (1−p)b)`.
  The sign is fixed by a − b for any baseline p.
- **Multi-isoform case:** `Δπ_t = π_t (f_t − F̄)/F̄` with `F̄ = Σ π_u f_u`. The sign is invariant across
  all possible baselines iff f_t ≥ max_u f_u or f_t ≤ min_u f_u (over u ≠ t with any weight).
- **Intron usage:** the sign of Δu_i is invariant iff
  min{f_t : t ∋ i} ≥ max{f_t : t ∋ j ≠ i in K(i)}, or the reverse inequality holds.
- **Quantities (descriptive; no hypothesis test):** for each v2 positive, compute
  1. whether the sign is baseline-invariant by structure;
  2. the range of Δ achievable over 1,000 Dirichlet(1) random baselines on the gene's transcripts;
  3. the observed fetal/adult/generic predictions.

## E3. Stage-swap with adult GTEx v8 Brain_Cortex sQTL

- **Data:** `Brain_Cortex.v8.sgenes.txt.gz` (grouped permutation, LeafCutter, hg38).
- **Positives:**
  - genes with qval < 0.05 that lie in the v2 gene universe (1,969 genes);
  - lead variant an SNV that lifts to hg19 (UCSC hg38ToHg19 chain, MD5-verified) with hg19 REF equal
    to the GTEx REF;
  - lead intron lifted to hg19 with length preserved, annotated in GENCODE v31lift37 and with at
    least one annotated competitor;
  - lead variant within 500 nt of an annotated site of the gene.
- **Intron coordinate convention:** chosen by the GT-AG rate on hg19 (the convention with the highest
  rate) before any prediction.
- **Endpoints:** direction accuracy (all; SpliceAI-informative ≥ 0.05) and signed Spearman, for the
  fetal, adult, generic, gtex_cortex and uniform baselines. No controls (GTEx tested-variant lists
  are not available), so no AUPRC.
- **Primary contrast (pre-specified):** interaction
  I = (fetal − adult)_Wen − (fetal − adult)_GTEx for direction accuracy (all) and for signed
  Spearman. Bootstrap independently within each arm (1-Mb blocks) and difference the draws.
- **Prediction if the developmental hypothesis holds:** I > 0. Two contrasts, Bonferroni α = 0.025.
- **Secondary:** within-GTEx fetal − adult.

## E4. Alternative fetal baseline (Patowary et al. 2024)

- **Data:** `patowary_shared_transcript_counts.tsv.gz`. This was derived by J. Dirawi from Patowary
  Data S2 (Dryad doi:10.5061/dryad.0rxwdbs68) by pooling 24 libraries; GENCODE-matched (ENST)
  transcripts only.
- **Conditions:** `patowary_fetal` (fetal_reads), `patowary_gz`, `patowary_cp`, normalized exactly
  like the other baselines.
- **Endpoints:**
  1. Agreement with the Leung fetal baseline (Pearson r on smoothed weights, evaluation genes), and
     coverage.
  2. The same v2 metrics on the v2 sQTL/isoQTL sets and on the E3 GTEx set.
- **Contrasts:** patowary_fetal − adult, patowary_fetal − leung_fetal; reported with CIs. These are
  exploratory: no correction family, and they are labelled as exploratory.
- **Leakage note:** Patowary and Wen et al. are both PsychENCODE (Gandal lab). Donor overlap is
  possible but cannot be verified. No genotype or QTL information enters the baseline.

## E5. Measured-effect validation of the propagation equation (MFASS)

- **Data:** MFASS `snv_data_clean.txt` (Chong et al. 2019; KosuriLab/MFASS commit 9a8e4f2).
- **Measured quantities:**
  - natural exon inclusion `nat_v2_index` (the baseline p);
  - variant `v2_dpsi` (the measured change);
  - minigene reporter in HEK293T: **not brain, not developmental**.
- **Variant set:**
  - `category == "mutant"` with non-missing `v2_dpsi` and `nat_v2_index`, where the hg19 base at
    `snp_position` equals the reference allele (or its complement for − strand, decided by checking
    both);
  - a random sample of **3,000** variants (numpy seed 20261104), drawn before scoring.
- **Exon sites:** acceptor = first exon base and donor = last exon base. The coordinate convention is
  fixed by checking that the REF SpliceAI acceptor/donor probability peaks at these positions in the
  first 100 scored variants. That check uses only REF tracks, never measured effects.
- **Predictors** (all from the same SpliceAI tracks):
  - P0 `raw`: −(Δp_acc + Δp_don) at the exon's own sites (signed SpliceAI site change);
  - P1 `prop_measured`: two-isoform propagation with p = measured `nat_v2_index`, a = r_acc·r_don,
    b = 1 (the reporter's flanking sites are shared by both isoforms);
  - P2 `prop_constant`: the same with p = the mean `nat_v2_index` of the training exons;
  - P3 `prop_shuffled`: the same with p permuted across exons (seed fixed). This is a negative
    control for the baseline;
  - P4 `prop_half`: p = 0.5.
- **Primary metrics:** Spearman and Pearson between predicted and measured ΔPSI; MAE; direction
  accuracy on variants with |measured ΔPSI| ≥ 0.1.
- **Uncertainty:** bootstrap over exons, so variants of the same exon are resampled together.
- **Primary contrasts:**
  - P1 − P2 (does the true baseline beat a constant baseline?);
  - P1 − P3 (does it beat a shuffled baseline?);
  - for Pearson r and MAE: four tests, Bonferroni α = 0.0125.
- **Calibrated variant (P5, the "improved model"):**
  - Δlogit = k·log(a). k is fitted by least squares on ΔPSI for a random 50% of **exons** (seed);
    all metrics are reported on the held-out 50% only.
  - P5 is compared with P1 (k = 1) on the held-out exons.

## E6. Locus-level evidence (Wen nominal blocks)

- Only about 50 introns have nominal association blocks on Drive. This is too few for inference, so
  they are inventoried only and not analyzed.


## E6. Per-donor Patowary baselines (added 2026-11-04, before any E6 computation)

- **Data:** Patowary et al. 2024 Data S2 sheet B (per-library TALON counts; 214,516 transcripts ×
  24 libraries; donors 209, 334, 336 × 4 libraries × GZ/CP). The user downloaded it into
  `data/external/patowary/` (extracted; archive checksums therefore not verifiable).
- **Integrity checks:**
  - sheet dimensions vs the audit (214,516 × 24);
  - library sums vs J. Dirawi's pooled table (per-ENST fetal/gz/cp reads must match exactly).
- **Conditions:** `patowary_d209`, `patowary_d334`, `patowary_d336` (each donor's 8 libraries
  summed), normalized like every other baseline. Propagation as in v2.
- **Analyses (exploratory, descriptive; no correction family):**
  1. Agreement of smoothed baselines on the evaluation genes for: libraries within donor-zone
     (technical), donors (biological), zones, pooled Patowary vs Leung fetal, Leung fetal vs adult.
  2. Predicted-sign disagreements between donors and vs adult.
  3. Metrics (direction accuracy, Spearman) per donor on Wen sQTL, Wen isoQTL and GTEx sQTL, with
     the between-donor range and fetal-donor − adult contrasts (block bootstrap, B = 2000).


## E11. Context-matched baseline test with measured magnitudes (Vex-seq) (registered 2026-11-04, before any E11 prediction)

- **Data:** Vex-seq (Adamson et al. 2018, *Genome Biol* 19:71; GEO GSE113163 processed PSI,
  3 replicates per cell line; design table = Additional file 1).
- **Variants:** 1,960 SNVs in 110 alternative exons, REF = hg19 for 100%. Measured in K562 and HepG2.
- **Baseline p_c:** the reference ("Control_Original") exon PSI in cell line c (mean of replicates).
- **Outcome:** ΔPSI_c = mutant PSI − p_c.
- **Not brain, not developmental. No overlap with the QTL loci** (assessment E10). This tests
  whether a context-matched baseline improves quantitative prediction in a measured system.
- **Model:** the two-isoform propagation exactly as E5,
  p' = p a/(p a + 1 − p), a = r_acc · r_don, ε = 0.01. Exon sites come from the design
  ("Exon up"/"Exon down", strand); the coordinate convention is fixed by REF-only SpliceAI peaks
  (modal offset).
- **Predictors per cell line c:**
  - M (matched baseline p_c);
  - S (swapped baseline, the other cell line's p);
  - C (constant baseline = mean p_c over exons);
  - R (raw SpliceAI site change);
  - N (null: p_c with permuted a).
- **Primary contrasts** (pooled over both cell lines, exon-cluster bootstrap B = 2000):
  1. M − S, Pearson r with measured ΔPSI;
  2. M − C, Pearson;
  3. the **between-context test**: Pearson between predicted (Δ̂_K562 − Δ̂_HepG2) under matched
     baselines and measured (ΔPSI_K562 − ΔPSI_HepG2), compared with the N null.
  Bonferroni α = 0.0167.
- **Secondary:** MAE; the same contrasts restricted to the 32 exons whose baselines differ by ≥ 0.1;
  calibrated k fitted on 50% of exons (seed 20261104), evaluated on the rest.

## E12. iPSC → NGN2-neuron long-read isoform time course (GSE310131) (registered 2026-11-04, before any overlap, baseline or prediction was computed)

**What had been inspected before registration:**
- the GEO count matrix schema and checksum;
- the GEO SOFT sample metadata;
- the Zenodo record's file list and the schemas of `final_classification.parquet`,
  `final_expression.parquet`, `full_expression.parquet`, `final_transcripts.gtf` and
  `id_to_sample.txt`.

**Data.**
- PacBio Revio Kinnex long-read RNA-seq from one engineered iPSC line (PGPC17_11-CLYBL-TetO-NGN2,
  one neurotypical male donor). Preprint: bioRxiv 10.64898/2025.12.16.693263; processed data:
  Zenodo 10.5281/zenodo.16897791.
- Time points and samples:
  - t00 iPSC: 3 samples, one per iPSC stock (lineage).
  - t04 intermediate: 6 samples, 2 wells per lineage.
  - t30 NGN2-induced neuron: 6 samples, 2 wells per lineage.
- Annotation: isoseq collapse, pigeon classify vs GENCODE v47, GRCh38.

**Not equivalent to fetal or adult brain.**
- t00, t04 and t30 are labelled `gse_t00_iPSC`, `gse_t04_intermediate` and `gse_t30_NGN2neuron`.
- NGN2 induction bypasses the cortical progenitor stage. t04 is not a radial-glia or NPC stage
  (despite the "NPC" column names), and t30 neurons are not fetal or adult cortex.

**Units of replication.**
- One donor, so donor-level biological variance is **not estimable**.
- The three lineages (independent iPSC stocks with staggered differentiations) are the highest-level
  replicate.
- Wells within a lineage are derivative or technical replicates.
- Lineage membership is read from the sample name (`X_L_w`), confirmed against the GEO titles
  ("derived from iPSC replicate L"). It is checked empirically: within-lineage correlation should
  exceed between-lineage correlation.

**E12a. Integrity.**
- The GEO matrix equals the Zenodo `full_expression.parquet`.
- `final_expression` ⊂ `full_expression` with identical counts.
- The published filter (> 5 reads in ≥ 3 samples) is reproduced.
- Every final isoform has GTF structure and a classification.
- Units: raw FLNC read counts per isoform.

**E12b. Overlap**, reported as exact counts:
- **Join:** isoform → transcript by FSM `associated_transcript` with the version stripped. GENCODE v47
  ENSTs absent from v31lift37 are counted as unmappable.
- **Universe and baseline overlap:**
  - the v2 gene universe (1,969 genes);
  - the annotated transcripts in evaluation genes;
  - Leung fetal/adult ENSTs with reads;
  - Patowary ENSTs with reads;
  - the 1,670 cluster transcripts of the 147 Wen sQTL positives.
- **Lead introns:** the 147 lead introns at junction level, using all 182,371 isoforms (including
  novel ones). GTF junctions are lifted hg38 → hg19 with the MD5-verified UCSC chain; both ends
  must lift to the same chromosome with intron length preserved.

**E12c. Coverage and reproducibility.**
- Per stage: cluster transcripts with 0 reads, and lead introns with no transcript-level and no
  junction-level support.
- Agreement of within-gene proportions at the 1,670 cluster transcripts, stratified by abundance
  exactly as in E7:
  - lineage vs lineage, with wells summed per lineage;
  - well vs well within a lineage, labelled technical/derivative.
- Stage contrast at the lead introns: |u_t30 − u_t00| and |u_t30 − u_t04|, compared with Leung
  |u_fetal − u_adult|, restricted to lead introns supported in both stages.

**E12d. Baselines.**
- **Primary:** FSM counts summed per ENST, pooled over all samples of a stage.
- **Secondary:** per-lineage t30 baselines, and an FSM + ISM sensitivity variant.
- Same smoothing and propagation as v2 (α = 1e-3, ε = 0.01, ±500 nt, component clusters).

**E12e. Prediction sensitivity**, for the Wen sQTL (147) and isoQTL (270) positives and GTEx sQTL
(82):
- **Pairs:** t30 vs t00, t30 vs t04, t04 vs t00, t30 vs Leung fetal, t30 vs Leung adult, t30 vs
  Patowary pooled, and lineage vs lineage at t30.
- **Sign:** sign flips counted only where both |pred| > τ, with τ ∈ {1e-12, 1e-6, 1e-3}.
- **Magnitude:** the |a|/|b| ratio among both > 1e-12 (quantiles; fraction outside [0.5, 2] and
  [0.1, 10]); the number of positives with |pred| > 1e-3.
- **Rank:** Spearman of predictions and of |pred|; top-20 overlap by |pred|.

**E12f. Benchmark (descriptive only; no confirmatory test).**
- Direction accuracy with ties for |pred| ≤ 1e-12, and Spearman ρ with Wen z, with block-bootstrap
  CIs for t30 − t00 and t30 − Leung fetal.
- E9 showed this benchmark has about 6% power, so no result here can support or refute a stage
  effect. It is reported for completeness only.

**E12g. Validation feasibility.**
- Determine whether the dataset contains any experimentally measured effect of a specific genetic
  variant: perturbation, allelic series, or allele-resolved counts.
- Differences in expression between differentiation stages are **not** variant effects and will not
  be used as ground truth.

**Interpretation rules fixed now.**
- **Improves the local-sparsity problem:** E12 can be said to improve it only if t30 reduces lead
  introns unsupported at transcript level below Leung fetal (66/147) *and* lineage agreement in the
  1–9-read stratum exceeds the Leung replicate value (Pearson 0.06).
- **Not material for the stage hypothesis:** if sign flips among both > 1e-6 between t30 and t00 are
  ≤ 5% of positives, stage differences in this resource are not material for the direction
  benchmark, whatever the magnitude differences.

## E13. Diagnostics for a generalizable model (registered 2026-11-04)

See `generalized_model_design/PREREGISTRATION_E13.md`, written before any E13 outcome was computed.

## Deviations log

**D1.** 2026-11-04, before any E5 evaluation (only REF/ALT SpliceAI tracks existed, and no
metric had been computed).
- **Change:** predictor P0 is +(Δp_acc + Δp_don), not −(…).
- **Reason:** a predicted site loss (Δp < 0) implies exon skipping, i.e. negative measured ΔPSI. The
  registered sign was a transcription error.
- **Effect:** |r| is unchanged; only the sign convention of direction accuracy is affected.

**D2.** Same time, an additional control. It is not a primary contrast and not part of the
Bonferroni family.
- **New predictor:** P6 `null_permuted_factor` uses the measured baseline p with the SpliceAI factor
  a permuted across variants (seed fixed).
- **Reason:** measured ΔPSI = v2_index − nat_v2_index shares the measured baseline. Noise in p and
  the [−p, 1−p] bounds can create a correlation between any p-dependent predictor and ΔPSI with no
  sequence information. P6 estimates that artifact.

**D3.** Exploration of the Drive folder after registration found the per-library Patowary Data S2
(Dryad) but could not download it (Dryad returned 403 to scripted clients; the Drive copy exceeds the
connector limit). E4 therefore uses only the pooled derived table, as registered.

**D4.** Citation correction only: the MFASS first author is Chong R (Chong et al. 2019, *Mol Cell*
73:183–194.e8, doi:10.1016/j.molcel.2018.10.037; verified via Crossref), not "Cheung". No analysis
change.

**D5.** Correction of a factual statement in the inventory (no analysis change). The Patowary
annotation is hg19 (Bulk.gtf introns 98.4% GT-AG on hg19; GENCODE lift37-style IDs), not hg38 as
first stated. All Patowary joins used ENST IDs only, so no result changes.

**D6.** 2026-11-04, after the E11 primary results.
- **Added:** a post-hoc robustness check, `E11b_vexseq_replicate_split.py`. The baseline comes from
  one replicate and the target ΔPSI from the other two, averaged over the 3 rotations.
- **Reason:** in the registered E11 design, the baseline and the target share replicate noise, which
  could favour the matched predictor.
- **Status:** not part of the Bonferroni family and reported alongside, not instead of, E11. All
  E11 conclusions held; the between-context estimate fell from 0.201 to 0.140.

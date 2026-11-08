# Pre-registration E14: external comparators on identical observations, with a never-touched test set

**Status: FROZEN 2026-11-05** (see the end of this document). Drafted the same day.
- **E14a selection has been run**, to fix the size of T-fresh. Its output files contain label columns,
  but no label was read, summarized or compared with any prediction.
- No scoring for E14b or analysis for E14c runs before the freeze.

**What had been seen before drafting:**
- All E13 results (`DIAGNOSTIC_RESULTS_E13.md`).
- Setup checks of the comparator tools (`COMPARATOR_TOOLS_SETUP.md`):
  - three COMPASS variants (peak positions);
  - 200 COMPASS **training-family** variants (MMSplice sign check, r ≈ 0.6);
  - four training-family variants for DeltaSplice timing and baseline dependence.
- Row counts of the candidate fresh pool (no outcomes).

**No comparator prediction has been compared with any test-set label.**

## E14 aim

Answer §9.3 item 1 and gate G-X. On identical observations:
1. Does the pre-specified method (A_SA) beat each external comparator in its native mode?
2. Does propagation through the measured reference (rule A) improve each sequence model over its raw
   score?

The primary answer comes from a COMPASS test set that has **never been used in any sequence-model
fitting or evaluation**.

## E14a. Test sets

| Set | Role | Definition | Status |
|---|---|---|---|
| **T-fresh** | **Primary** | COMPASS families selected by `E14a_fresh_compass_select.py` (rules below) | Never used for any sequence-model fit or evaluation. Its labels entered E13b only in aggregate, in a variance decomposition with no sequence model. That is stated in any report |
| T-E13 | Secondary (previously seen) | E13c COMPASS test families (77) | Seen in E13d (A, K, B1, R) |
| T-Vex | Secondary (previously seen) | Vex-seq, K562 + HepG2, 3 replicate rotations (as E13d) | Seen in E10, E11, E11b, E13d. **MMSplice's Δlogitψ model was fitted on Vex-seq**, so contrasts involving MMSplice on T-Vex are reported but not interpreted |
| T-MFASS | Secondary (previously seen) | MFASS E5 sample of 3,000, HEK293T, replicate split both directions (as E13d) | Seen in E5 and E13d |

**T-fresh rules** (applied before any outcome; counts only):
1. Single SNVs with a WT row; ≥ 20 reads in both replicates for variant and WT, in HEK and K562
   (identical to E13c).
2. hg38 → hg19 liftover with exon length preserved; REF must match hg19.
3. Exclude every event **and gene** in the E13c subset, train or test.
4. Exclude exons within ±5 bp of any Vex-seq or MFASS exon, using E13a's exon definitions and
   `overlap()`.
5. **Sequence-similarity filter (§4 T2).** Exclude an exon if the 8-mer containment of its ±100 nt
   window in any E13c exon's ±100 nt window is ≥ 0.15, a proxy for paralogs with ≥ ~80% identity.
   In the candidate pool, the maximum containment has median 0.041 and 99th percentile 0.079, so the
   threshold sits well above background.
6. Keep families with ≥ 5 QC singles, at most **CAP** variants per family (seed 20261106). CAP is fixed
   in §E14d from compute alone.

**Resulting size** (`fresh_selection.json`, counts only):

| Step | Events |
|---|---|
| QC-passing events with a WT row | 748 |
| After removing E13c events and genes | 585 |
| Removed as identical to Vex-seq or MFASS exons (rule 4) | 51 |
| Removed by the similarity filter (rule 5) | 2 |
| Kept: ≥ 5 singles, CAP = 20 | **129 families, 125 genes, 1,592 variants** |

Most kept families are ones E13c excluded for lying within 100 kb of Vex-seq or MFASS. Nearness alone
is not leakage for a COMPASS-internal test, because no E14 model was fitted on data at these loci;
rules 4–5 remove identical and paralogous exons.

## E14b. Predictors

### Shared rules
- **Baseline.** Every predictor that uses a baseline receives the same replicate-split baseline:
  - COMPASS: WT PSI from replicate r_in, target from r_tg, both orders, 5 cells.
  - Vex-seq: 3 rotations, as E11b.
  - MFASS: R1→R2 and R2→R1.
  - Baselines are clipped to [0.01, 0.99] for all methods (as E13d and D-E13-3).
- **Exon-site convention.** Fixed per tool and dataset by the REF-only modal-offset rule of E11/E13d,
  applied to that tool's own reference tracks.
- **Labels (COMPASS)** exactly as in E13d:
  - Per-cell QC: ≥ 20 reads in both replicates for variant and WT.
  - y = raw target-replicate PSI_var − PSI_WT.
- **Eligibility.** A row is eligible only if both exon sites have finite SpliceAI ratios, i.e. lie
  within ±500 nt of the variant (as E13d). It must also be scored by every tool (a REF or window
  failure in any tool removes the row for all of them; counts reported). The set is identical for all
  predictors. MMSplice gives Δ = 0 for variants more than 100 nt from the
  exon (`in_range = False`). They stay in the analysis, and the in-range subset is reported as
  secondary.

### Predictor table

| ID | Tool | Mode | Prediction of ΔPSI |
|---|---|---|---|
| R_SA | SpliceAI v1.3.1 | raw | Δacc + Δdon at the exon sites (= E13d R) |
| R_PG | Pangolin | raw | ΔP(splice) at the acceptor + donor, mean over 4 tissues |
| R_DS | DeltaSplice | raw, reference-free | ΔSSU at the acceptor + donor |
| **A_SA** | SpliceAI | rule A | p·a/(p·a + 1 − p) − p, a = r_acc·r_don, ε = 0.01 (= E13d A). **Pre-specified "our method"** (G-B failed, so A stays the default; §9.5) |
| A_PG | Pangolin | rule A | same, using Pangolin ratios (mean of 4 tissues) |
| A_DS | DeltaSplice | rule A | same, using reference-free DeltaSplice ratios |
| N_MM | MMSplice 2.4.0 | **native** | expit(logit p + Δlogitψ) − p |
| N_DS | DeltaSplice | **native** exon mode | acceptor SSU = donor SSU = p (the authors' convention in `generate_mutdata.py`); ΔPSI = mean(Δacc, Δdon) at the exon sites. Computation plan in §E14d |
| K, B1 | SpliceAI | E13d fits | parameters frozen from E13d (k = 2.16; b = −0.25, k_acc = 1.65, k_don = 2.53). Secondary |

Pangolin "brain" and other single tissues, and the alternative DeltaSplice weights, are not used. If
reported at all, they are labelled exploratory.

## E14c. Analysis

**Metrics.**
- Pearson r with measured ΔPSI in the target replicate: primary.
- MAE: co-primary.
- AUPRC for |ΔPSI| ≥ 0.1: secondary. Raw scores (R_*) enter Pearson and AUPRC only.

**Uncertainty.** Exon-family cluster bootstrap, B = 1,000, seed 20261106. Paired: the same resamples
for all predictors.

**Primary contrasts, on T-fresh only.**
- **P1 (gate G-X):** A_SA minus each external comparator in its native mode: N_MM, N_DS and R_PG
  (Pangolin has no reference-using mode). That is 3 contrasts, Holm-adjusted at family-wise α = 0.05.
- **P2:** A_X − R_X for X ∈ {SA, PG, DS}. That is 3 contrasts, Holm-adjusted.

**Gate G-X (from §9.5), evaluated on T-fresh.** Pass only if A_SA beats **each** of N_MM, N_DS and R_PG
by ≥ 0.02 in r, with the Holm-adjusted 95% bootstrap CI of each difference excluding 0. Requiring every
comparator avoids picking a "best comparator" on the test data.
- **Fail**, but no comparator significantly better than A_SA: report A_SA as **not distinguishable**
  from existing tools. The paper is then a benchmark/decomposition (§9.2, rank 1).
- **Fail**, with a comparator significantly better: report that as the finding.

**Secondary.**
- The same contrasts on T-E13, T-Vex and T-MFASS, labelled "previously seen".
- Per-cell results.
- MAE and AUPRC.
- The MMSplice in-range subset.
- Variants inside versus outside the exon.

**Not allowed:**
- choosing the "our method" after seeing T-fresh results;
- tuning any parameter on any test set;
- dropping a comparator after seeing results;
- re-running T-fresh with different rules, unless logged as a deviation with the reason and both
  results reported.

## E14d. Compute plan (fixed from timings only)

**Timing inputs.** Measured on 4 COMPASS training-family variants; 4 CPU threads.

| Item | Measurement |
|---|---|
| DeltaSplice reference-free tracks | 23–26 s per variant |
| DeltaSplice native exon mode | 46–59 s per pass (one baseline value) |
| Exon-site offsets | 0 at both sites |
| Native ΔPSI across baselines 0.05 / 0.2 / 0.5 / 0.8 / 0.95 | Varies strongly and smoothly. Example: −0.07 / −0.16 / −0.27 / −0.29 / −0.26 |
| Pangolin | ~5 s per variant |
| MMSplice | ~1.5 s per variant |
| SpliceAI | ~0.19 variants/s |

**CAP = 20 variants per family.** It was fixed from compute alone, and families, not variants, drive
the cluster-bootstrap precision. Result: **T-fresh = 129 families, 125 genes, 1,592 variants**
(`results/E14_comparators/fresh_selection.json`).

**Native DeltaSplice (N_DS) on CPU.** One pass per observation is infeasible: T-fresh has 5 cells × 2
splits × 1,592 variants, about 16k passes, roughly 220 h. The pre-specified plan is a grid with
interpolation, run on a subsample.

1. **Grid.** For each variant, run the 5 baselines ψ ∈ {0.05, 0.2, 0.5, 0.8, 0.95}. Interpolate
   ΔPSI(ψ) with monotone cubic (PCHIP) interpolation on logit ψ, and evaluate it at each observation's
   replicate-split baseline. Clamp to the grid ends outside [0.05, 0.95].
2. **Accuracy check, before any test use.** Use 20 random variants from **E13c training families**
   (seed 20261106). Add 2 exact passes each at ψ = 0.35 and 0.65. **Pass** if the median absolute
   interpolation error is ≤ 0.005 and the 95th percentile is ≤ 0.02. If it fails, N_DS is computed
   only with a GPU (exact, one pass per observation), and P1 waits. Cost: about 2 h.
3. **Subsample.** 40 T-fresh families, drawn at random (seed 20261106), with 5 variants each (random,
   same seed): 200 variants, about 14 h. The N_DS contrast in P1 is computed on **this subsample only**,
   with all predictors restricted to the same rows. The other P1 and P2 contrasts use all of T-fresh.
   With a GPU, N_DS is run exact on all of T-fresh; that result is primary and the subsample result is
   reported alongside.

**Order and budget on CPU.**

| Step | Content | Time | Cumulative |
|---|---|---|---|
| 1 | Freeze this document | — | — |
| 2 | Score T-fresh with SpliceAI (2.3 h), Pangolin (2.2 h), MMSplice (0.7 h) and DeltaSplice tracks (10.6 h) | 16 h | 16 h |
| 3 | N_DS accuracy check (2 h), then the N_DS subsample (14 h) | 16 h | 32 h |
| 4 | Score the secondary sets T-E13, T-Vex and T-MFASS with Pangolin (~9 h) and MMSplice (~3 h). SpliceAI scores already exist | 12 h | 44 h |
| 5 | Run E14c once | — | — |

DeltaSplice on the secondary sets, for both R_DS/A_DS and N_DS, runs **only with a GPU**. Otherwise it
is reported as not run.

No outcome is examined until all of steps 2–4 are complete, so results are evaluated once, by a single
script.

## Deviations log (E14)

**C-E14-1** (2026-11-05, after the freeze, before any scoring output existed; clarifications, not
changes). Fixed in the docstring of `E14c_comparators_analysis.py`:
- **Gate G-X criterion.** For each comparator, the point estimate of A_SA − comparator must be ≥ 0.02,
  **and** the Holm-adjusted bootstrap CI must exclude 0.
- **Holm-adjusted CI.** The contrast with the k-th smallest bootstrap p-value (out of m) uses a
  two-sided percentile interval at level 1 − 0.05/(m − k + 1).
- **G-X with N_DS missing.** If N_DS is unavailable (failed accuracy check, no GPU), G-X is reported as
  **not evaluable**, not as passed.

**C-E14-2** (2026-11-05). MMSplice inputs for the secondary sets use:
- Vex-seq: 1-based inclusive "Exon up/down" (108/110 exons AG…GT);
- MFASS: E5's acceptor/donor bounds (lengths match; 51% AG…GT, the known MFASS ambiguity seen in E5
  and E13d).

They are used unchanged, as in E13d.

**Code-path test.** It will be run once with `E14_TEST_MODE=1`: labels permuted, B = 20, nothing
written.

**Frozen:** 2026-11-05, on the user's approval in chat ("Freeze and start").
- The SHA-256 of this file as frozen (this line included) is in
  `results/E14_comparators/PREREG_E14.sha256`.
- Later changes go only in the deviations log.

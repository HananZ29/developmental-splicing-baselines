# Developmental splicing baselines for variant-effect prediction

Code, pre-registered protocols and result files for a study testing whether a fetal-brain splicing
baseline improves prediction of splice-altering variant effects compared with adult-brain, generic
(GTEx) and uniform baselines, with the sequence predictor (SpliceAI) and the propagation model held
fixed. Evaluation is at the level of junction usage (prenatal sQTL associations) and exon inclusion
(reporter assays), not full-length isoform proportions.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Shared modules: SpliceAI runner, propagation model, splice graph, isoform proportions |
| `experiments/01`–`08` | Earlier long-read baseline pipeline; its outputs (`results/`) are inputs to later scripts |
| `publication_extension/scripts/` | Experiments E2–E15 and the manuscript figures |
| `publication_extension/results/` | Result files (JSON, TSV, logs) written by those scripts |
| `publication_extension/figures/paper/` | Manuscript Figures 2–6 (`figE7_vexseq.png`: Figure 3a, b) |
| `publication_extension/**/PREREGISTRATION*.md`, `*.sha256` | Pre-registered protocols and their SHA-256 freeze locks |
| `tests/` | Unit tests (propagation model; E15 pattern EM and identifiability) |

Scripts load each other by relative path (for example `E15_predict.py` imports
`E15_baseline_diagnostics.py`; `E14c_comparators_analysis.py` loads the E13d, E11 and E5 scripts).
Keep the layout and run every command from the repository root.

## Installation

Main environment: Python 3.12.6 (developed on Windows 11).

```
python -m pip install -r requirements.lock.txt
python -m pytest tests -q
```

Comparator tools use separate environments, one per tool, pinned in `tools/envs/`:

| Tool | Source | Lock file |
|---|---|---|
| Pangolin | github.com/tkzeng/Pangolin, commit 5cf94b8d | `tools/envs/requirements_pangolin.lock.txt` |
| MMSplice | PyPI (`mmsplice`) | `tools/envs/requirements_mmsplice.lock.txt` |
| DeltaSplice | github.com/chaolinzhanglab/DeltaSplice, commit 71c2fd16 | `tools/envs/requirements_deltasplice.lock.txt` |

Clone Pangolin and DeltaSplice into `tools/src/` and create the environments in `tools/envs/<tool>/`.
Checkpoint hashes are in `publication_extension/generalized_model_design/COMPARATOR_TOOLS_SETUP.md`.
SpliceAI v1.3.1 weights (`spliceai1.h5`–`spliceai5.h5`, from the SpliceAI distribution) go in
`data/external/spliceai/`. They are loaded by `src/spliceai_runner.py`.

## Data (not redistributed)

Download the datasets below into `data/`, using the paths named at the top of each script.

| Dataset | Source | Used in |
|---|---|---|
| BrainSpan MRF alignments, neocortex (362 samples, 42 donors) | download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10 (fetched by `E15_fetch_brainspan_mrf.py`) | E15 |
| BrainSpan/PsychENCODE sample metadata | development.psychencode.org | E15 |
| GTEx v8 junction counts, transcript TPM, Brain_Cortex sQTL | GTEx portal, release v8 | E15, earlier pipeline |
| Wen et al. 2024 SuSiE fine-mapping, conditional variants, intron annotation | Synapse project syn50897018 (login required) | E15 |
| Wen et al. 2024 sQTL/isoQTL leads | github.com/gandallab/devBrain_xQTL | Earlier pipeline, E2–E12 |
| Leung et al. 2021 cortex Iso-Seq | Zenodo deposit cited in Leung et al., *Cell Reports* 2021 | Earlier pipeline |
| Patowary et al. 2024 | Dryad doi:10.5061/dryad.0rxwdbs68 | E6, E7 |
| GSE310131 | GEO; Zenodo 10.5281/zenodo.16897791 | E12 |
| MFASS | github.com/KosuriLab/MFASS (commit 9a8e4f2) | E5, E13, E14 |
| Vex-seq | GEO GSE113163 | E10, E11, E13, E14 |
| COMPASS | GEO GSE307247 | E13, E14 |
| GENCODE v10, v19, v31lift37; hg19 genome; UCSC hg38ToHg19 chain | GENCODE; UCSC | All |

Sizes and MD5 checksums of the E15 inputs: `publication_extension/developmental_stage/DATA_PROVENANCE_E15.md`.
`data/processed/E15_baselines/cluster_eligibility.tsv.gz` is included because the E15 freeze lock
lists it.

## Reproducing the results

All commands are run from the repository root. The `.sh` drivers need a Bash shell (for example Git Bash).

| Result | Commands |
|---|---|
| Earlier long-read pipeline | `python experiments/01_build_proportions.py`, then `02`, `03`, `04`, `05`, `06`, `07`, `08` in order |
| E15: fetal vs adult, generic, uniform (primary) | `bash publication_extension/scripts/run_E15.sh` (step list and manual inputs in its header). Post-hoc bootstrap checks: `python -E publication_extension/scripts/E15_posthoc_robustness.py` |
| E5: MFASS | `python -E publication_extension/scripts/E5_mfass_select.py`, then `E5_mfass_evaluate.py` |
| E10, E11: Vex-seq | `E10_vexseq_assess.py`, `E11_vexseq_evaluate.py`, `E11b_vexseq_replicate_split.py` |
| E13: COMPASS | `E13a_label_audit.py`, `E13b_compass_context_decomposition.py`, `E13b2_headroom_scale.py`, `E13c_compass_select.py`, `E13d_models_A_B_generalization.py` |
| E14: predictor benchmark | `E14a_fresh_compass_select.py`; scoring as in `publication_extension/results/E14_comparators/run_fresh_scoring.sh`, `run_steps3_4.sh`, `run_step4_parallel.sh`; then `python -E publication_extension/scripts/E14c_comparators_analysis.py analyse` and `E14_verify.py` |
| E12: differentiation time course | `E12a_gse310131_overlap_coverage.py`, `E12b_gse310131_predictions.py`, `E12c_independent_checks.py` |
| Earlier-model diagnostics | `E2_root_cause.py`, `E3_gtex_select.py`, `E34_stage_swap_and_patowary.py`, `E6_patowary_donors.py`, `E7_local_coverage.py`, `E9_sensitivity_power.py` |
| Figures | `python -E publication_extension/scripts/make_paper_figures.py` (Figures 2–6); `make_figures.py` (Figure 3a, b) |

Scripts without a path prefix are in `publication_extension/scripts/`.

## Frozen protocols

E15 and E14 were frozen with SHA-256 locks before their outcomes were analysed:

```
sha256sum -c publication_extension/developmental_stage/PREREG_E15.sha256
sha256sum -c publication_extension/results/E14_comparators/PREREG_E14.sha256
```

The E14 check reports one expected mismatch, `PREREGISTRATION_E14.md`, explained in
`publication_extension/results/E14_comparators/LOCK_VERIFICATION_2026-11-07.md`. `.gitattributes`
stores every file byte for byte so that the locks stay valid.

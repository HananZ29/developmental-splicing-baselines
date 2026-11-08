#!/usr/bin/env bash
# E15 driver: the commands that produced publication_extension/results/E15_devstage, in execution order.
#
# Reconstructed on 2026-11-07 from the session transcript and the output and log file names. Each
# command matches what was executed on 2026-11-05/06, with two exceptions:
#   - the SpliceAI commands (steps 8, 10) are inferred from score_spliceai.py's usage line and the log
#     file names, because the original command lines are not in the surviving transcript;
#   - step 14 calls E15_gate0_mei.py, a verbatim transcription of inline code.
# Re-running overwrites the frozen E15 outputs. Do NOT run steps 1-15 against the frozen results; use a
# copy of the workspace. Step 16 is the only step that reads QTL statistics.
#
# Main environment: C:\Python312 (Python 3.12.6); see requirements.lock.txt.
# Manual inputs that must already exist before step 4:
#   - data/external/wen_synapse/  (4 Wen et al. 2024 files, downloaded by the user from Synapse;
#                                  see developmental_stage/DATA_PROVENANCE_E15.md for sizes and MD5s)
#   - data/external/gencode/gencode.v31lift37.annotation.gtf.gz, GTEx v8 junction counts, hg19 FASTA
set -euo pipefail
cd "$(dirname "$0")/../.."
S=publication_extension/scripts
R=publication_extension/results/E15_devstage
D=publication_extension/developmental_stage

# --- Data acquisition and verification (outcome-blind)
python -I $S/E15_fetch_brainspan_mrf.py data/external/brainspan_mrf/neocortex     # 1  BrainSpan MRF (resumable)
python -I $S/E15_mrf_junction_counts.py                                             # 2  junction counts per sample
python -E $S/E15_annotation_harmonization_audit.py                                  # 3  annotation audit -> annotation_audit.json
python -I $S/E15_verify_wen_files.py                                                # 4  Wen file size/readability/MD5

# --- Baselines and eligibility (outcome-blind)
python -E $S/E15_gtex_generic_counts.py                                             # 5  GTEx generic baseline counts
python -E $S/E15_baseline_diagnostics.py > $R/log_diagnostics.txt 2>&1              # 6  E1-E4, LODO, bootstrap, bias

# --- Candidate variants and SpliceAI scores (outcome-blind; fine-mapping columns 0-7 only)
python -E $S/E15_candidates.py                                                      # 7  -> candidates.tsv, to_score.tsv
python -E $S/score_spliceai.py $R/to_score.tsv $R/spliceai > $R/log_spliceai.txt 2>&1   # 8  (inferred command)
python -E $S/E15_gtex_candidates.py > $R/log_gtex_candidates.txt 2>&1               # 9  GTEx cortex crossover candidates
python -E $S/score_spliceai.py $R/gtex_to_score.tsv $R/spliceai_gtex > $R/log_spliceai_gtex.txt 2>&1  # 10 (inferred)

# --- Predictions and pre-freeze checks (outcome-blind)
python -E $S/E15_predict.py > $R/log_predict.txt 2>&1                               # 11 -> predictions.tsv
python -E $S/E15_prefreeze_report.py                                                # 12 -> prefreeze_report.json
python -E $S/E15_gate0_power.py $R/predictions.tsv > $R/log_gate0.txt 2>&1          # 13 Gate 0 -> gate0_power.json
python -E $S/E15_gate0_mei.py                                                       # 14 power at the MEI -> gate0_mei.json
mkdir -p $R/gate0_subset_diff001                                                    # 15 Gate 0 on |pred_F1 - pred_A1| > 0.01
python -E -c "
import pandas as pd;p=pd.read_csv('$R/predictions.tsv',sep='\t');s=p[(p.pred_F1-p.pred_A1).abs()>0.01];s.to_csv('$R/gate0_subset_diff001/predictions_subset.tsv',sep='\t',index=False);print(len(s))"
python -E $S/E15_gate0_power.py $R/gate0_subset_diff001/predictions_subset.tsv > $R/gate0_subset_diff001/log.txt 2>&1

# --- Freeze check and the single confirmatory analysis (reads QTL statistics)
sha256sum -c $D/PREREG_E15.sha256                                                   # lock must verify
python -E $S/E15_analysis.py > $R/log_analysis.txt 2>&1                             # 16 -> E15_results.json

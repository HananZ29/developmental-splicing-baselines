#!/usr/bin/env bash
# E14 steps 3-4 (frozen plan). Run after run_fresh_scoring.sh has finished.
set -e
cd /c/Users/Hanan/Desktop/bioInformatics-research
O=publication_extension/results/E14_comparators
S=publication_extension/scripts
DS=tools/envs/deltasplice/Scripts/python
echo "step3 start $(date)"
$DS -E $S/E14b_ds_native_grid.py check > $O/log_ds_check.txt 2>&1
python -E $S/E14c_comparators_analysis.py check > $O/log_ds_check_eval.txt 2>&1
if python -E -c "import json,sys; sys.exit(0 if json.load(open('$O/ds_check.json'))['pass'] else 1)"; then
  echo "N_DS check passed; running fresh subsample $(date)"
  $DS -E $S/E14b_ds_native_grid.py fresh > $O/log_ds_fresh.txt 2>&1
else
  echo "N_DS check FAILED: N_DS deferred to GPU per pre-registration $(date)"
fi
echo "step4 start $(date)"
for n in e13test vexseq mfass; do
  tools/envs/pangolin/Scripts/python -E $S/score_pangolin.py $O/${n}_to_score.tsv $O/pangolin_${n} tools/src/Pangolin > $O/log_pangolin_${n}.txt 2>&1; echo "pangolin $n done $(date)"
  tools/envs/mmsplice/Scripts/python -E $S/score_mmsplice.py $O/${n}_mmsplice_in.tsv $O/mmsplice_${n}.tsv > $O/log_mmsplice_${n}.txt 2>&1; echo "mmsplice $n done $(date)"
done
echo "steps 3-4 done $(date)"

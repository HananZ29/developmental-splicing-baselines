#!/usr/bin/env bash
# E14 step 4 run in parallel with the step-3 N_DS subsample (execution-schedule change only; 2026-11-07).
# Same scripts, inputs and output paths as run_steps3_4.sh step 4. Pangolin is resumable from parts/.
cd /c/Users/Hanan/Desktop/bioInformatics-research
O=publication_extension/results/E14_comparators
S=publication_extension/scripts
echo "parallel step4 start $(date)"
for n in e13test vexseq mfass; do
  ( OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 tools/envs/pangolin/Scripts/python -E $S/score_pangolin.py $O/${n}_to_score.tsv $O/pangolin_${n} tools/src/Pangolin > $O/log_pangolin_${n}.txt 2>&1; echo "pangolin $n done exit=$? $(date)" ) &
  ( OMP_NUM_THREADS=1 tools/envs/mmsplice/Scripts/python -E $S/score_mmsplice.py $O/${n}_mmsplice_in.tsv $O/mmsplice_${n}.tsv > $O/log_mmsplice_${n}.txt 2>&1; echo "mmsplice $n done exit=$? $(date)" ) &
done
wait
echo "parallel step4 done $(date)"

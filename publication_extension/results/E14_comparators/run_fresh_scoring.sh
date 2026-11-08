#!/usr/bin/env bash
# E14 step 2 (frozen plan): score T-fresh with SpliceAI, Pangolin, MMSplice, DeltaSplice tracks, sequentially.
set -e
cd /c/Users/Hanan/Desktop/bioInformatics-research
O=publication_extension/results/E14_comparators
S=publication_extension/scripts
echo "start $(date)"
python -E $S/score_spliceai.py $O/fresh_to_score.tsv $O/spliceai_fresh > $O/log_spliceai_fresh.txt 2>&1; echo "spliceai done $(date)"
tools/envs/pangolin/Scripts/python -E $S/score_pangolin.py $O/fresh_to_score.tsv $O/pangolin_fresh tools/src/Pangolin > $O/log_pangolin_fresh.txt 2>&1; echo "pangolin done $(date)"
tools/envs/mmsplice/Scripts/python -E $S/score_mmsplice.py $O/fresh_mmsplice_in.tsv $O/mmsplice_fresh.tsv > $O/log_mmsplice_fresh.txt 2>&1; echo "mmsplice done $(date)"
tools/envs/deltasplice/Scripts/python -E $S/score_deltasplice.py $O/fresh_to_score.tsv $O/deltasplice_fresh tools/src/DeltaSplice > $O/log_deltasplice_fresh.txt 2>&1; echo "deltasplice done $(date)"
echo "all done $(date)"

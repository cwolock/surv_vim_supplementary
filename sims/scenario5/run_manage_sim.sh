#!/bin/bash
ml fhR
# num_combos is number of parameter combinations (see manage_sim.R)
# 2 is the number of reps per combos
# 3 is number of reps performed per job
# 1 is sim name
num_combos=36

njobs=`expr 500 \* $num_combos`

sbatch --array=1-$njobs -p short -t 00:45:00 -e ./iotrash/s-%A_%a.out -o ./iotrash/s-%A_%a.out /home/cwolock/surv_vim_supplementary/sims/scenario5/call_manage_sim.sh

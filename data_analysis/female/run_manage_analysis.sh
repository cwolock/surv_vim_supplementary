#!/bin/bash
num_combos=2
njobs=`expr 50 / 1 \* $num_combos`
sbatch --array=1-$njobs -e ./iotrash/s-%A_%a.out -o ./iotrash/s-%A_%a.out /home/cwolock/surv_vim_supplementary/data_analysis/female/call_manage_analysis.sh

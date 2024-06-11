#!/bin/bash
ml R/4.4.0-gfbf-2023b

num_combos=2

njobs=`expr 10 \* $num_combos`

sbatch --array=1-$njobs -e ./iotrash/s-%A_%a.out -o ./iotrash/s-%A_%a.out /home/cwolock/surv_vim_supplementary/data_analysis/female/call_manage_analysis.sh

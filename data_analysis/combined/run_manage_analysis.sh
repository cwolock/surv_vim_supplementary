#!/bin/bash

num_combos=2

njobs=`expr 10 / 1 \* $num_combos`

sbatch --array=1-$njobs -e ./iotrash/s-%A_%a.out -o ./iotrash/s-%A_%a.out /home/cwolock/surv_vim_supplementary/data_analysis/combined/call_manage_analysis.sh

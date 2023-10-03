#!/bin/bash

_now=$(date '+%Y_%m_%d__%H_%M_%S')

srun -p xeon -w jacob --exclusive --cpu-freq=performance -n 1 -c 48 -t $2:0:0 ./$1.sh > bench_out/$1_$_now.txt 2>&1
#srun -p xeon -w bouillon --exclusive --cpu-freq=performance -n 1 -c 48 -t $2:0:0 ./$1.sh > bench_out/$1_$_now.txt 2>&1

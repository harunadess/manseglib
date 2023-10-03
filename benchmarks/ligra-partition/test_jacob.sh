#!/bin/bash

_now=$(date '+%Y_%m_%d__%H_%M_%S')

srun -p xeon -w jacob -n 1 -c 8 -t 6:0:0 ./$1.sh > test_out/$1_$_now.txt
#srun -p xeon -w jacob -n 1 -c 8 -t 10 ./$1.sh > test_out/$1_$_now.txt

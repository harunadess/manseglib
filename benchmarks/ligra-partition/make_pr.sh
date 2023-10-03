#!/bin/bash
# setup environment
# module load compilers/gcc-4.9.0

#export LD_LIBRARY_PATH="../../cilk-swan/lib/"  # req. lib
export CILK_NWORKERS=4                                             # no. cpu threads
#export LD_PRELOAD="./bin/interposer_cilk.so"                       # req. for cilk

# note: large graphs should take about ~15 minutes to run
# any longer than that and something has went wrong

PR="PageRankUpdate PageRankUpdate_Floats PageRankManSeg PageRankManSeg_f PageRankUpdate_F2D PageRankManSeg_hybrid"
#PR="PageRankManSeg_dev"

BASE_PATH="graphs"

date
make ${PR}
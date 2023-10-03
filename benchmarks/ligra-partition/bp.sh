#!/bin/bash
# setup environment
module load compilers/gcc-4.9.0

export LD_LIBRARY_PATH="/home/jsun/swanModel/jacob/cilk-swan/lib/" # i think we need these 3 lines to run things
export CILK_NWORKERS=8                                             # i.e. no. cpu threads to use (cilk is like intel omp)
export LD_PRELOAD="./bin/interposer_cilk.so"                       # and also this (shared object for use with cilk?)

# todo: run on larger graphs - either twitter/orkut/others that are large -> use more rounds/splits/more cores when running
# note: large graphs should take about ~15 minutes to run
# any longer than that and something has went wrong

make BP
echo "made bp"
make BPUpdate
echo "made bpupdate"
make BPManSeg
echo "made bpmanseg"

# todo: create benchmarking script
# i.e. jacob --all cores script (take param .sh)
# try with -c ([48/96]*16)

# use -rounds 10 and -c (no. cores * 16)

date
echo "ref bp start"
./BP -c 128 -v edge -rounds 2 /var/shared/projects/asap/graphs/adj/realworld/recent/USAroad_undir > bp_out/ref_USAroad_undir.txt 2>&1
date
echo "ref bp finish"

date
echo "bp update start"
./BPUpdate -c 128 -v edge -rounds 2 /var/shared/projects/asap/graphs/adj/realworld/recent/USAroad_undir > bp_out/USAroad_undir.txt 2>&1
date
echo "bp update end"

date
echo "bp manseg start"
./BPManSeg -c 128 -v edge -rounds 2 /var/shared/projects/asap/graphs/adj/realworld/recent/USAroad_undir > bp_out/USAroad_undir_msa.txt 2>&1
date
echo "bp manseg end"

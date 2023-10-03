#!/bin/bash
# setup environment
module load compilers/gcc-4.9.0

export LD_LIBRARY_PATH="/home/jsun/swanModel/jacob/cilk-swan/lib/" # i think we need these 3 lines to run things
export CILK_NWORKERS=8                                             # i.e. no. cpu threads to use (cilk is like intel omp)
export LD_PRELOAD="./bin/interposer_cilk.so"                       # and also this (shared object for use with cilk?)

# todo: run on larger graphs - either twitter/orkut/others that are large -> use more rounds/splits/more cores when running
# note: large graphs should take about ~15 minutes to run
# any longer than that and something has went wrong

make PageRank
echo "made pagerank"
make PageRankManSeg
echo "made pagerank manseg"

# todo: create benchmarking script
# i.e. jacob --all cores script (take param .sh)
# try with -c ([48/96]*16)


#
# So as it appears -c 768 is okay-ish. try the higher number
# What we have also discovered is how incredibly slow your library is
# The library was 60s slower in twitter_undir where the majority of the iterations
#  were completed as "head only" operations, which should be ~2x faster than doubles

# use -rounds 10 and -c (no. cores * 16)

#                            other graph paths
#           BEFORE USE: create pr_out/${name}.txt file for results
#

# /var/shared/projects/asap/graphs/adj/realworld/recent/RMAT27_undir_rmdup_undir
# /var/shared/projects/asap/graphs/adj/realworld/recent/twitter_undir_xl
# /var/shared/projects/asap/graphs/adj/realworld/uk_dir

date
echo "pagerank start"
./PageRank -c 768 -v edge -rounds 2 /var/shared/projects/asap/graphs/adj/realworld/recent/USAroad_undir > pr_out/USAroad_undir.txt 2>&1
date
echo "pagerank finish"

date
echo "pagerank manseg start"
./PageRankManSeg -c 768 -v edge -rounds 2 /var/shared/projects/asap/graphs/adj/realworld/recent/USAroad_undir > pr_out/USAroad_undir_msa.txt 2>&1
date
echo "pagerank manseg finish"

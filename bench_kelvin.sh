#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=16

#SBATCH -o /users/40150688/meng_project/pr_jacobi_pcg.txt
#SBATCH -e /users/40150688/meng_project/pr_jacobi_pcg.err

## load compiler we want to use
module add compilers/gcc/7.2.0

PROJ_DIR=/users/40150688/meng_project

echo "starting pagerank"
date
echo ""
./$PROJ_DIR/ligra-partition/pagerank_bench.sh > $PROJ_DIR/ligra-partition/pagerank.txt 2>&1
echo ""
echo "finish pagerank"
date
echo ""

echo ""
echo "starting jacobi_stencil"
date
echo ""
./$PROJ_DIR/jacobi_stencil/jacobi_bench.sh > $PROJ_DIR/jacobi_stencil/jacobi.txt 2>&1
echo ""
echo "finish jacobi_stencil"
date
echo ""

echo ""
echo "finish conjugate_gradient"
date
echo ""
./$PROJ_DIR/conjugate_gradient/conjugate_bench.sh > $PROJ_DIR/conjugate_gradient/conjugate_gradient.txt 2>&1
echo ""
echo "finish conjugate_gradient"
date
echo ""

date

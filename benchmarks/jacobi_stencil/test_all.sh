#!/bin/bash

export OMP_NUM_THREADS=8

ITER=5000

JACOBI="jacobi_omp jacobi_mod_omp"

echo ""
for prog in $JACOBI;
do
    echo "start ${prog}"
    ./${prog} $ITER > ${prog}_out.txt 2>&1
    echo "end ${prog}"
    echo ""
done

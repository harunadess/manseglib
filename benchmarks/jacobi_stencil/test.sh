#!/bin/bash

make J_jacobi_up J_jacobi_omp

export OMP_NUM_THREADS=32

echo ""
echo "seq:"
./J_jacobi_up $1

echo ""
echo "omp:"
./J_jacobi_omp $1

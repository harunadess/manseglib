#!/bin/bash

export OMP_NUM_THREADS=8

make jacobi_mod jacobi_mod_omp
echo ""
echo "manseg:"
./jacobi_mod $1

echo "manseg w/omp:"
./jacobi_mod_omp $1

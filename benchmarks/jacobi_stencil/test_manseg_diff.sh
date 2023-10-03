#!/bin/bash

echo "diff jacobi_mod | jacobi_mod_omp:"
diff --speed-large-files --strip-trailing-cr --suppress-common-lines \
    ./jacobi_mod.txt ./jacobi_mod_omp.txt | grep '>' -c
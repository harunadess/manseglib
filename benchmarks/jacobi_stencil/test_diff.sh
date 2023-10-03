#!/bin/bash

echo "diff jacobi_u | jacobi_omp : "
diff --speed-large-files --strip-trailing-cr --suppress-common-lines \
    ./jacobi_u.txt ./jacobi_omp.txt | grep '>' -c
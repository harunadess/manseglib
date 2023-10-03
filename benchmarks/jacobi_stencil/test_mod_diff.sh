#!/bin/bash

echo "diff jacobi_up | jacobi_mod:"
diff --speed-large-files --strip-trailing-cr --suppress-common-lines \
    ./jacobi_u.txt ./jacobi_mod.txt | grep '>' -c
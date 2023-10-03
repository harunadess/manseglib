#!/bin/bash

make J_jacobi_up jacobi_mod
echo ""
echo "up:"
./J_jacobi_up $1

echo "mod:"
./jacobi_mod $1

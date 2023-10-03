#!/bin/bash

## default small graphs ##
# std

make pagerank omp_pagerank

export OMP_NUM_THREADS=8

BASE_PATH="/mnt/d/Users/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/mantissa-segmentation/snap_graphs"

echo "std version"

## SNAP ##
echo "serial version"
# CSC
./pagerank SNAP CSC "$BASE_PATH/roadNet-PA.txt"

echo "parallel version"
# CSC
./omp_pagerank SNAP CSC "$BASE_PATH/roadNet-PA.txt"


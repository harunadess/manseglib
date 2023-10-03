#!/bin/bash

## todo: figure out how to do fixed number of iterations: i.e. only max 20k total

MTX= #"bcsstk01 bcsstk06 bcsstk13 bcsstk16 bcsstk19 bcsstk24 bcsstk28 bcsstm25 bcsstm27"
MTX="bcsstk01 bcsstk06 bcsstk13 " #bcsstk16 bcsstk19 bcsstk24 bcsstk28 bcsstm25 bcsstm27"
TOL=1e-7
CG_ITER=2000
IR_ITER=100
STEP_CHECK=1000

for mat in $MTX
do
    echo "start ${mat}"
    ./sparsesolve ../data/${mat}.mtx $IR_ITER $TOL $CG_ITER $TOL $STEP_CHECK > ../results/manseg_${mat}.out
    echo "end ${mat}"
    echo ""
done
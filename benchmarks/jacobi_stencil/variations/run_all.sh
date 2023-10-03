#!/bin/bash

VAR="nb16_b512 nb32_b128 nb32_b256 nb32_b512 nb64_b64 nb64_b128 nb128_b32 nb128_b64 nb256_b32 nb512_b16"
JS="jacobi_omp jacobi_mod_omp"

ITER=10000

for var in $VAR;
do
	cd ${var}
	make $JS
	cd ..
done

echo "finished making"
echo ""

for var in $VAR;
do
	echo "runs for ${var}"
	date
	for j in $JS;
	do
		echo "start ${j}"
		date
		./${var}/${j} $ITER > ./${var}/${j}_out.txt
		date
		echo "end ${j}"
		echo ""
	done
	date
	echo "finished runs for ${var}"
	echo "----------------------------------------"
done

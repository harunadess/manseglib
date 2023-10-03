#/!bin/bash

cd mpir && make && cd ..

cd mpir_class && make && cd ..

cd mpir_class_manseg && make && cd ..

./run_both.sh

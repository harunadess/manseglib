CCX = g++
CCXFLAGS = -std=c++11 -O3 -mavx2 -fopenmp
LDFLAGS = -lm -lrt -fopenmp

.SECONDARY:

.PHONY: build

build: sparsesolve

sparsesolve: sparsesolve.o mmio.o matrix.o cg.o ir.o
	$(CCX) $^ $(LDFLAGS) -o $@

clean:
	rm -rf *.o sparsesolve


test: sparsesolve
	./sparsesolve ../data/bcsstk01.mtx 1 1e-7 10000 1e-7 100

%.o: %.cpp cg.h vector.h matrix.h
	$(CCX) $(CCXFLAGS) -c $< -o $@

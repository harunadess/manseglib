CC = gcc
CFLAGS = -std=gnu99 -Wall -O3 -mavx2 -fopenmp
LDFLAGS = -lm -lrt -fopenmp

.SECONDARY:

.PHONY: build

build: sparsesolve

sparsesolve: sparsesolve.o mmio.o matrix.o cg.o ir.o
	$(CC) $^ $(LDFLAGS) -o $@

clean:
	/bin/rm -rf tags core *.o sparsesolve

tags: *.c *.h
	ctags *.c *.h

test: sparsesolve
	./sparsesolve ../data/bcsstk01.mtx 1 1e-7 10000 1e-7 100

%.o: %.c cg.h vector.h matrix.h
	$(CC) $(CFLAGS) -c $< -o $@


CCX=g++
CCXFLAGS=-O3 -std=c++17 -mavx2

all: 
	make main assem

main: main.o matrix.o
	${CCX} -o main main.o matrix.o -lnuma

.PHONY: main.o
main.o: main.cpp ../manseglib.hpp matrix.h vector.h
	${CCX} ${CCXFLAGS} -c main.cpp -lnuma

matrix.o: matrix.cpp ../manseglib.hpp
	${CCX} ${CCXFLAGS} -c matrix.cpp


assem: assem.o
	${CCX} -o assem assem.o

.PHONY: assem.o
assem.o: assem.cpp ../manseglib.hpp
	${CCX} ${CCXFLAGS} -c assem.cpp

.PHONY: clean
clean:
	rm main  *.o

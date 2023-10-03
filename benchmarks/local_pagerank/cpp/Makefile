
CCX=g++
CCXFLAGS=-O3 -std=c++11 -mavx2

all: 
	make pagerank pagerank_check msa_pagerank

pagerank: pagerank.o
	${CCX} -o pagerank pagerank.o

.PHONY: pagerank.o
pagerank.o: pagerank.cpp quicksort.h
	${CCX} ${CCXFLAGS} -c pagerank.cpp

pagerank_check: pagerank_check.o
	${CCX} -o pagerank_check pagerank_check.o

.PHONY: pagerank_check.o
pagerank_check.o: pagerank_check.cpp quicksort.h
	${CCX} ${CCXFLAGS} -c pagerank_check.cpp

msa_pagerank: msa_pagerank.o
	${CCX} -o msa_pagerank msa_pagerank.o

msa_pagerank.o: msa_pagerank.cpp ../../../manseglib.hpp quicksort.h
	${CCX} ${CCXFLAGS} -c msa_pagerank.cpp

omp_pagerank: omp_pagerank.o
	${CCX} -fopenmp -o omp_pagerank omp_pagerank.o

.PHONY: omp_pagerank.o
omp_pagerank.o: omp_pagerank.cpp quicksort.h
	${CCX} ${CCXFLAGS} -fopenmp -c omp_pagerank.cpp

omp_pagerank_manseg: omp_pagerank_manseg.o
	${CCX} -fopenmp -o omp_pagerank_manseg omp_pagerank_manseg.o

.PHONY: omp_pagerank_manseg.o 
omp_pagerank_manseg.o: omp_pagerank_manseg.cpp quicksort.h ../../../manseglib.hpp
	${CCX} ${CCXFLAGS} -fopenmp -c omp_pagerank_manseg.cpp

.PHONY: clean
clean:
	rm pagerank pagerank_check msa_pagerank omp_pagerank omp_pagerank_manseg*.o
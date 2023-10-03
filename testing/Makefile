
CXX=g++
CXXFLAGS=-O3 -std=c++17 -mavx2

SEQ=array_element_copy array_operator_functions basic_read_write head_pair_basic_sum precision_switch type_conversion pico_pagerank pico_random_read pico_random_write
PARALLEL=pico_parallel_pagerank pico_parallel_random_read pico_parallel_random_write

all:
	make $(SEQ) $(PARALLEL)

sequential:
	make $(SEQ)

parallel:
	make $(PARALLEL)

pico_parallel_%: pico_parallel_%.o
	$(CXX) -fopenmp $^ -o $@

pico_parallel_.o: pico_parallel_.cpp util.h
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

%: %.o
	$(CXX) $^ -o $@

%.o: %.cpp util.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm $(SEQ) $(PARALLEL)
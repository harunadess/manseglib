#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "util.h"
#include "../manseglib.hpp"

using namespace ManSeg;
using namespace std;

template<class Arr>
void helper(Arr& ms, int size, int repeats, int* loc)
{
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            double val = (static_cast<double>(rand()) / RAND_MAX);
            ms[loc[j]] = val;
        }
    }
}

int main()
{	
	const int size = 1000000;
    const int repetitions = 1000;
    HeadsArray* f = new HeadsArray(size);
    PairsArray* t = new PairsArray(size);
    double* d = new double[size];
    
    srand(2517);
    int* loc = new int[size];
    for(int i = 0; i < size; ++i) loc[i] = rand() % size;

    auto d_start = chrono::_V2::high_resolution_clock::now();
    helper(d, size, repetitions, loc);
    auto d_end = chrono::_V2::high_resolution_clock::now();

    auto f_start = chrono::_V2::high_resolution_clock::now();
    helper(*f, size, repetitions, loc);
    auto f_end = chrono::_V2::high_resolution_clock::now();

    auto p_start = chrono::_V2::high_resolution_clock::now();
    helper(*t, size, repetitions, loc);
    auto p_end = chrono::_V2::high_resolution_clock::now();

    auto d_time = chrono::duration_cast<chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = chrono::duration_cast<chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = chrono::duration_cast<chrono::nanoseconds>(p_end - p_start).count()*1e-9;

    cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:      " << p_time << "s" << endl;

    delete[] d, f, t;

	return 0;
}

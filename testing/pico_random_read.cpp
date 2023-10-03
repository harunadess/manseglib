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

template<class T1, class T2>
void helper_2(double* d, T1& a, T2& b, int size, int repeats)
{
    for(int i = 0; i < repeats; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            a[j] = d[j];
			b[j] = d[j];
        }
    }
}

// benchmark reading
template<class Arr>
double test(Arr& ms, int size, int repetitions, int* loc)
{
    double t = 0.0;
    for(int i = 0; i < repetitions; ++i)
    {
        for(int j = 0; j < size; ++j)
            t += ms[loc[j]];
    }
    return t;
}

int main()
{	
	cout << fixed << setprecision(16);

	const int size = 1000000;
    const int repetitions = 1000;
    HeadsArray* f = new HeadsArray(size);
    PairsArray* t = new PairsArray(size);
    double* d = new double[size];
    
    srand(2517);
    int* loc = new int[size];
    for(int i = 0; i < size; ++i) loc[i] = rand() % size;

    helper(d, size, 1, loc);
    helper_2(d, *f, *t, size, 1);

    auto f_start = chrono::_V2::high_resolution_clock::now();
    double ftotal = test(*f, size, repetitions, loc);
    auto f_end = chrono::_V2::high_resolution_clock::now();

    auto p_start = chrono::_V2::high_resolution_clock::now();
    double ptotal = test(*t, size, repetitions, loc);
    auto p_end = chrono::_V2::high_resolution_clock::now();

    auto d_start = chrono::_V2::high_resolution_clock::now();
    double dtotal = test(d, size, repetitions, loc);
    auto d_end = chrono::_V2::high_resolution_clock::now();

    auto d_time = chrono::duration_cast<chrono::nanoseconds>(d_end - d_start).count()*1e-9;
    auto f_time = chrono::duration_cast<chrono::nanoseconds>(f_end - f_start).count()*1e-9;
    auto p_time = chrono::duration_cast<chrono::nanoseconds>(p_end - p_start).count()*1e-9;

    cout << "std double: " << d_time << "s\nheads only: " << f_time << "s\npairs:      " << p_time << "s\ntotals" << endl;
    cout << "std double: " << dtotal << "\nheads only: " << ftotal << "\npairs:      " << ptotal << endl;

    f->del(); t->del();
    delete[] d, f, t, loc;

	return 0;
}

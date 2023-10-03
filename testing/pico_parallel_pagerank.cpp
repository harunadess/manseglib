#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>
#include <omp.h>

#include "util.h"
#include "../manseglib.hpp"

using namespace ManSeg;
using namespace std;

constexpr double val = 0.85;

template<class T1, class T2, class T3>
void fill_pr(T1& a1, T1& a2, T2& b1, T2& b2, T3& c1, T3& c2, int* v, const int& size)
{
	for(int i = 0; i < size; ++i)
	{
		v[i] = (i + 1);

		b1[i] = (double)(i+1)/0.1;
        b2[i] = 0.0;

        c1[i] = (double)(i+1)/0.1;
        c2[i] = 0.0;
        
        a1[i] = (double)(i+1)/0.1;
        a2[i] = 0.0;    
	}
}

template<class Arr>
double pr(Arr& a, Arr& b, int* v, int* loc, const int& size)
{
	auto dst = chrono::_V2::high_resolution_clock::now();

	#pragma omp parallel for
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) b[i] += val*(a[loc[j]]/v[j]);
    }

    auto ded = chrono::_V2::high_resolution_clock::now();
    auto dt = chrono::duration_cast<chrono::nanoseconds>(ded - dst).count()*1e-9;

	return dt;
}

template<class Arr>
double pr_with_acc(Arr& a, Arr& b, int* v, int* loc, const int& size)
{
    double* acc = new double[size];

    auto mxst = chrono::_V2::high_resolution_clock::now();
	#pragma omp parallel for
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j) acc[i] += val*(a[loc[j]]/v[j]);

        b[i] = acc[i];
    }
    auto mxed = chrono::_V2::high_resolution_clock::now();
    auto mt = chrono::duration_cast<chrono::nanoseconds>(mxed - mxst).count()*1e-9;

    delete[] acc;
	return mt;
}

int main()
{
	const int size = 200000;

	int num_threads_avail;
	num_threads_avail = omp_get_max_threads();

	cout << "running with: " << num_threads_avail << " threads" <<  endl;
	omp_set_num_threads(num_threads_avail);

	HeadsArray* x = new HeadsArray(size);
    HeadsArray* xy = new HeadsArray(size);       // x, xy = heads
    PairsArray* y = new PairsArray(size);
    PairsArray* yy = new PairsArray(size);       // y, yy = pairs

	int* v = new int[size];
    double* d = new double[size];
    double* dy = new double[size];              // d, dy = doubles

	srand(2517);
	// randomise data accesses
    int* loc = new int[size];
    for(int i = 0; i < size; ++i)
        loc[i] = (rand() % size);

	fill_pr(d, dy, *x, *xy, *y, *yy, v, size);

	double dt;
    std::cout << "round 1 - read and write direct \n";
    // doubles
    dt = pr(d, dy, v, loc, size);
    std::cout << "std double = " << dt << "s\n";
    // heads
    dt = pr(d, dy, v, loc, size);
	std::cout << "heads ver = " << dt << "s\n";
    // pairs
    dt = pr(d, dy, v, loc, size);
	std::cout << "pairs ver = " << dt << "s\n";

    std::cout << "round 2 - read from manseg, write to double accumulator \n";
    // reset
    fill_pr(d, dy, *x, *xy, *y, *yy, v, size);

    // doubles
    dt = pr_with_acc(d, dy, v, loc, size);
    std::cout << "std double = " << dt << "s\n";
    // heads
    dt = pr_with_acc(d, dy, v, loc, size);
	std::cout << "heads ver = " << dt << "s\n";
    // pairs
    dt = pr_with_acc(d, dy, v, loc, size);
	std::cout << "pairs ver = " << dt << "s\n";

    delete[] d, dy, v, loc;
    x->del(); xy->del();
    y->del(); yy->del();

	return 0;
}

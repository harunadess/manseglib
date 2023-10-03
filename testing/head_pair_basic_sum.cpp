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
double test(Arr& a, Arr& b, size_t size)
{
	double sum = 0;
	for(size_t i = 0; i < size; ++i)
		sum += (a[i] + b[i]);

	cout << "sum =\t" << sum << "        ";
	printBinary(sum);

	return sum;
}

int main()
{
	cout << fixed << setprecision(16);
	const size_t size = 10000000UL;

	HeadsArray arr(size);
    HeadsArray arr2(size);

    double* darr = new double[size];
    double* darr2 = new double[size];

    srand(2517);
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        darr[i] = val;
        arr.set(i, val);

        val += 2.0;

        arr2.set(i, val);
        darr2[i] = val;
    }

    cout << "std double" << endl;
	cout << "round 1:  ";
	double dsum1 = test(darr, darr2, size);
	cout << "round 2:  ";
	double dsum2 = test(darr, darr2, size);

	cout << "heads only" << endl;
	cout << "round 1:  ";
	double hsum1 = test(arr, arr2, size);
	cout << "round 2:  ";
	double hsum2 = test(arr, arr2, size);

	PairsArray parr = arr.createFullPrecision();
	PairsArray parr2 = arr2.createFullPrecision();
	for(size_t i = 0; i < size; ++i)
    {
        parr.set(i, darr[i]);
        parr2.set(i, darr2[i]);
    }

	cout << "pairs" << endl;
	cout << "round 1:  ";
	double psum1 = test(parr, parr2, size);
	cout << "round 2:  ";
	double psum2 = test(parr, parr2, size);

	int return_code = 0;
	if(dsum1 != dsum2)
	{
		cerr << "dsum1 != dsum2 -- test failed" << endl;
		return_code = 1;
	}

	if(hsum1 != hsum1)
	{
		cerr << "hsum1 != hsum2 -- test failed" << endl;
		return_code = 1;
	}

	if(psum1 != psum2)
	{
		cerr << "psum1 != psum2 -- test failed" << endl;
		return_code = 1;
	}

	if((dsum1 != psum1) || (dsum2 != psum2))
	{
		cerr << "psum != dsum -- test failed" << endl;
		return_code = 1;
	}

	delete[] darr, darr2;
	arr.del(); arr2.del();

	if(return_code == 0)
		cout << "test passed!" << endl;

	return return_code;
}

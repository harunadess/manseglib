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

double v1 = 1.5;				   		// 1.5
double v2 = 1.100000000000000089;  		// 1.1
double v3 = -2.333333000000000101;  	// -2.333333
int v4 = 4;
uint64_t v5 = (1UL << 33);		   		// 8589934592
bool v6 = 1;							// 1
float v7 = 3.4f;						// 3.4f
char v8 = 'B';							// 66

template<class ArrType>
void func(ArrType& p)
{
	cout << "(double)v1  = " << (double)v1 << " --> "; printBinary((double)v1);
	p.set(0, v1);
	cout << "(ArrType)v1 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v2  = " << (double)v2 << " --> "; printBinary((double)v2);
	p.set(0, v2);
	cout << "(ArrType)v2 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v3  = " << (double)v3 << " --> "; printBinary((double)v3);
	p.set(0, v3);
	cout << "(ArrType)v3 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v4  = " << (double)v4 << " --> "; printBinary((double)v4);
	p.set(0, v4);
	cout << "(ArrType)v4 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v5  = " << (double)v5 << " --> "; printBinary((double)v5);
	p.set(0, v5);
	cout << "(ArrType)v5 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v6  = " << (double)v6 << " --> "; printBinary((double)v6);
	p.set(0, v6);
	cout << "(ArrType)v6 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v7  = " << (double)v7 << " --> "; printBinary((double)v7);
	p.set(0, v7);
	cout << "(ArrType)v7 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

	cout << "(double)v8  = " << (double)v8 << " --> "; printBinary((double)v8);
	p.set(0, v8);
	cout << "(ArrType)v8 = " << p.read(0) << " --> "; printBinary(p.read(0));
	cout << endl;

}

int main()
{
	cout << fixed << setprecision(16);

	HeadsArray h(1);
	PairsArray p(1);

	cout << "ArrType = HeadsArray\n";
	func(h);
	cout << "=============================================================================================================================\n" << endl;
	cout << "ArrType = PairsArray\n";
	func(p);

	h.del();
	p.del();

	return 0;
}
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

// should demonstrate use of all array operators for double precision number 

template<class ArrType, typename T>
int eq_func(ArrType& a, T* t)
{
	int return_code = 0;
	
	a[0] = t[0];
	if(a[0] != t[0])
	{
		cerr << "a[0] != t[0]" << endl;
		cerr << "expected " << t[0] << ", actual = " << a[0] << endl;
		return_code = 1;
	}
	a[1] = t[1];
	if(a[1] != t[1])
	{
		cerr << "a[1] != t[1]" << endl;
		cerr << "expected " << t[1] << ", actual = " << a[1] << endl;
		return_code = 1;
	}
	a[2] = t[2];
	if(a[2] != t[2])
	{
		cerr << "a[2] != t[2]" << endl;
		cerr << "expected " << t[2] << ", actual = " << a[2] << endl;
		return_code = 1;
	}

	if(return_code == 1)
		cerr << "failed = op" << endl;

	return return_code;		
}

template<class ArrType, typename T>
int op_func(ArrType& a, T* t)
{
	int return_code = 0;

	// testing the global operators also tests the compound operators
	// as they are implemented using the global operators

	t[0] = t[0] + 2.0;
	a[0] = a[0] + 2.0;
	if(a[0] != t[0])
	{
		cerr << "a[0] != t[0]" << endl;
		cerr << "expected " << t[0] << ", actual = " << a[0] << endl;
		cerr << "failed + op" << endl;
		return_code = 1;
	}

	t[1] = t[1] - 1.0;
	a[1] = a[1] - 1.0;
	if(a[1] != t[1])
	{
		cerr << "a[1] != t[1]" << endl;
		cerr << "expected " << t[1] << ", actual = " << a[1] << endl;
		cerr << "failed - op" << endl;
		return_code = 1;
	}

	t[2] = t[2] * 2.0;
	a[2] = a[2] * 2.0;
	if(a[2] != t[2])
	{
		cerr << "a[2] != t[2]" << endl;
		cerr << "expected " << t[2] << ", actual = " << a[2] << endl;
		cerr << "failed * op" << endl;
		return_code = 1;
	}

	t[0] = t[0] / 3.0;
	a[0] = a[0] / 3.0;
	if(a[0] != t[0])
	{
		cerr << "a[0] != t[0]" << endl;
		cerr << "expected " << t[0] << ", actual = " << a[0] << endl;
		cerr << "failed / op" << endl;
		return_code = 1;
	}

	return return_code;		
}

int main()
{
	double d[] = {1.0, 1.5, 3.0};
	PairsArray p(3);
	HeadsArray h(3);

	cout << "PairsArray\n";
	int return_code = eq_func(p, d);
	if(return_code == 0)
		cout << "test passed!" << endl;
	else
		cerr << "test failed!" << endl;

	return_code = op_func(p, d);
	if(return_code == 0)
		cout << "test passed!" << endl;
	else
		cerr << "test failed!" << endl;

	// reset values
	d[0] = 1.0; d[1] = 1.5; d[2] = 3.0;
	
	cout << "\nHeadsArray\n";
	return_code = eq_func(h, d);
	if(return_code == 0)
		cout << "test passed!" << endl;
	else
		cerr << "test failed!" << endl;

	return_code = op_func(h, d);
	if(return_code == 0)
		cout << "test passed!" << endl;
	else
		cerr << "test failed!" << endl;
	
	p.del();
	h.del();
	return 0;
}
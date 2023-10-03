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

int main()
{
	cout << fixed << setprecision(16);
	const size_t size = 100UL;

	PairsArray arr(size);
    PairsArray arr2(size);
    double d[size];

    srand(2517);
    for(size_t i = 0; i < size; ++i)
    {
        double val = (static_cast<double>(rand()) / RAND_MAX);
        d[i] = val;
        arr.setPair(i, val);
		arr2[i] = val;
    }

	int return_code = 0;
    for(size_t i = 0; i < size; ++i)
    {
		if(d[i] != arr.read(i))
		{
			cout << "d[" << i << "](" << d[i] << ") != arr.read(" << i << ")(" << arr.read(i) << ")!\n";
			return_code = 1;
		}
		if(d[i] != arr2[i])
		{
			cout << "d[" << i << "](" << d[i] << ") != arr2[" << i << "](" << arr2[i] << ")!\n";
			return_code = 1;
		}
		if(arr.read(i) != arr2[i])
		{
			cout << "arr.read(" << i << ")(" << arr.read(i) << ") != arr2[" << i << "](" << arr2[i] << ")!\n";
			return_code = 1;
		}
    }

    arr.del();
    arr2.del();

	

	if(return_code == 0)
		cout << "test passed!" << endl;

	return return_code;
}

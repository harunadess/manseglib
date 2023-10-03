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

	int size = 10;

	int return_code = 0;

	HeadsArray h(size);
	PairsArray p(size);
	PairsArray dest(size);

	srand(2517);
	for(int i = 0; i < size; ++i)
	{
		h[i] = (static_cast<double>(rand()) / RAND_MAX);
		p[i] = (static_cast<double>(rand()) / RAND_MAX);
	}

	for(int i = 0; i < size; ++i)
	{
		double cache = h[i];
		dest[i] = h[i];
		if(h[i] != dest[i])
		{
			cerr << "h[" << i << "] != dest[" << i << "]\n";
			cerr << "expected = " << dest[i] << ", actual = " << h[i] << "\n";
			return_code = 1;
		}
		if(cache != h[i])
		{
			cerr << "h[" << i << "] != cache\n";
			cerr << "expected = " << cache << ", actual = " << h[i] << "\n";
			return_code = 1;
		}

		cache = p[i];
		dest[i] = p[i];
		if(p[i] != dest[i])
		{
			cerr << "p[" << i << "] != dest[" << i << "]\n";
			cerr << "expected = " << dest[i] << ", actual = " << p[i] << "\n";
			return_code = 1;
		}
		if(cache != p[i])
		{
			cerr << "p[" << i << "] != cache\n";
			cerr << "expected = " << cache << ", actual = " << p[i] << "\n";
			return_code = 1;
		}
	}

    h.del();
    p.del();
	dest.del();

	if(return_code == 0)
		cout << "test passed!" << endl;
	else
		cerr << "test failed!" << endl;

	return return_code;
}
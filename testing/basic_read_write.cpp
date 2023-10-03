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
	double d[4] = {1.1, 2.0, 3.14, 5.5};
	HeadsArray h(4);
	PairsArray p(4);

	for(int i = 0; i < 4; ++i)
	{
		h[i] = d[i];
		p[i] = d[i];
	}

	cout << "ieee binary\n";
	for(int i = 0; i < 4; ++i)
	{
		cout << "actual value: " << d[i] << endl;
		cout << "d[" << i << "] = ";
		printBinary(d[i]);

		cout << "p[" << i << "] = ";
		printBinary(p[i]);

		cout << "h[" << i << "] = ";
		printBinary(h[i]);
	}

	h.del();
	p.del();

	return 0;
}
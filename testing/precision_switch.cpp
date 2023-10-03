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

	double d[] = {1.1, 2.0, 3.4, 6.25};

	ManSegArray m(4);

	for(int i = 0; i < 4; ++i)
		m.heads.setPair(i, d[i]);
	
	int return_code = 0;

	for(int i = 0; i < 4; ++i)
	{
		if(d[i] != m.pairs[i])
		{
			cerr << "d[" << i << "] != m.pairs[" << i << "]\n";
			cerr << "expected = " << d[i] << ", actual = " << m.pairs[i] << "\n";
		}
	}

	m.copytoIEEEdouble();
	for(int i = 0; i < 4; ++i)
	{
		if(m.heads[i] != m.full[i])
		{
			cerr << "m.heads[" << i << "] != m.full[" << i << "]\n";
			cerr << "expected = " << m.heads[i] << ", actual = " << m.full[i] << "\n";
		}
	}


	m.delSegments();
	m.del();

	if(return_code == 0)
		cout << "test passed !" << endl;
	else
		cerr << "test failed !" << endl;

	return 0;
}
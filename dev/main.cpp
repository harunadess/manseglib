#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <chrono>
#include <math.h>

#include "../manseglib.hpp"
// #include "matrix.h"
// #include "vector.h"


using namespace ManSeg;

void printBinary(double d)
{
    unsigned long l = *reinterpret_cast<unsigned long*>(&d);
    for(int i = 63; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 63 || i == 52)
            std::cout << " ";
        else if(i == 32)
            std::cout << "|";
    }
    std::cout << std::endl;
}

// print binary of float
void printBinary(float f)
{
    unsigned int l = *reinterpret_cast<unsigned int*>(&f);
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 31 || i == 23)
            std::cout << " ";
    }
    std::cout << std::endl;
}

// print binary of int
void printBinary(int l)
{
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i % 8 == 0 && i > 0)
            std::cout << " ";
    }
}

/* int main5(int argc, char** argv)
{
    const int m_size = 4;
    double x[m_size*m_size] = {
        0.1, 0.1, 2.5, 1.5,
        1.6, 0.2, 0.21, 0.71,
        2.5, 1.34, 1.22, 1.88,
        0.21, 0.11, 0.691, 2.
    };

    matr_ m(m_size);
    
    m.fill(x);

    double y[m_size] = {3.25, 3.1, 3.001, 3.77};
    double z[m_size] = {2.03, 2.44, 2.98, 2.25};
    
    m.mm(y, z);

    std::cout << "m = \n";

    ManSegArray msa(m_size);
    for(int i = 0; i < m_size; ++i)
    {
        for(int j = 0; j < m_size; ++j)
        {
            std::cout << m.v[i].pairs[j] << "  ";

            msa.pairs[j] = m.v[i].pairs[j];
        }
        std::cout << '\n';
    }

    std::cout << "\nz=\n";
    for(int i = 0; i < m_size; ++i)
        std::cout << z[i] << "  ";
    std::cout << '\n';

    std::cout << "\nmsa.pairs=\n";
    for(int i = 0; i < m_size; ++i)
        std::cout << msa.pairs[i] << "  ";
    std::cout << '\n';

    std::cout << "prod of m= " << prod(&m, m_size) << '\n';

    return 0;
} */

template<class TwoSegArray>
double sum(TwoSegArray& t, int n)
{
    double sum = 0.0;
    for(int i = 0; i < n; ++i)
    {
        sum += t[i];
    }
    return sum;
}
double sum(double* t, int n)
{
    double sum = 0.0;
    for(int i = 0; i < n; ++i)
    {
        sum += t[i];
    }
    return sum;
}

#define nbsize 2
ManSegArray m1[nbsize][nbsize];
ManSegArray m2[nbsize][nbsize];
double* d1[nbsize][nbsize];
double* d2[nbsize][nbsize];
int main()
{
	// allocate and fill
	for(int i = 0; i < nbsize; ++i)
	{
		for(int j = 0; j < nbsize; ++j)
		{
			d1[i][j] = new double[nbsize*nbsize];
			d2[i][j] = new double[nbsize*nbsize];
			
			m1[i][j].alloc(nbsize*nbsize);
			m1[i][j].full = new double[nbsize*nbsize];

			m2[i][j].alloc(nbsize*nbsize);
			m2[i][j].full = new double[nbsize*nbsize];

			for(int k = 0; k < nbsize; ++k)
			{
				for(int l = 0; l < nbsize; ++l)
				{
					d1[i][j][k*nbsize + l] = l+1;
					d2[i][j][k*nbsize + l] = l+1 + (nbsize*nbsize);
					m1[i][j].pairs[k*nbsize + l] = l+1;
					m2[i][j].pairs[k*nbsize + l] = l+1 + (nbsize*nbsize);
				}
			}
		}
	}
	
	for(int iteration = 0; iteration < 10; ++iteration)
	{
		// make some change to some values
		for(int i = 0; i < nbsize; ++i)
		{
			for(int j = 0; j < nbsize; ++j)
			{
				for(int k = 0; k < nbsize; ++k)
				{
					for(int l = 0; l < nbsize; ++l)
					{
						d2[i][j][k*nbsize + l] += l + (nbsize*nbsize);
						m2[i][j].pairs[k*nbsize + l] += l + (nbsize*nbsize);
					}
				}
			}
		}

		for(int i = 0; i < nbsize; ++i)
		{
			for(int j = 0; j < nbsize; ++j)
			{
				for(int k = 0; k < nbsize; ++k)
				{
					for(int l = 0; l < nbsize; ++l)
					{
						std::cout << k*nbsize + l << ":d1 " << d1[i][j][k*nbsize+l] << "| "; 
						std::cout << k*nbsize + l << ":m1 " << m1[i][j].pairs[k*nbsize+l] << "\t\t"; 
						std::cout << k*nbsize + l << ":d2 " << d2[i][j][k*nbsize+l] << "| "; 
						std::cout << k*nbsize + l << ":m2 " << m2[i][j].pairs[k*nbsize+l] << "\n"; 
					}
				}
				std::cout << std::endl;
			}
		}

		// swap does not work
		// causes double deallocation for some reason
		// and idk how 

		// this works here, but does not work
		// for jacobi, so..?
		std::cout << "*swap*" << std::endl;
		std::swap(*m1, *m2);
		std::swap(d1, d2);

		// print and things should be swapped
		for(int i = 0; i < nbsize; ++i)
		{
			for(int j = 0; j < nbsize; ++j)
			{
				for(int k = 0; k < nbsize; ++k)
				{
					for(int l = 0; l < nbsize; ++l)
					{
						std::cout << k*nbsize + l << ":d1 " << d1[i][j][k*nbsize+l] << "| "; 
						std::cout << k*nbsize + l << ":m1 " << m1[i][j].pairs[k*nbsize+l] << "\t\t"; 
						std::cout << k*nbsize + l << ":d2 " << d2[i][j][k*nbsize+l] << "| "; 
						std::cout << k*nbsize + l << ":m2 " << m2[i][j].pairs[k*nbsize+l] << "\n"; 
					}
				}
				std::cout << "\n" << std::endl;
			}
		}
	}

	// clean up
	for(int i = 0; i < nbsize; ++i)
	{
		for(int j = 0; j < nbsize; ++j)
		{
			m1[i][j].delSegments();
			m2[i][j].delSegments();
		}
	}

	return 0;
}

int main002()
{
    const int size = 10;
    ManSegArray a(size);

    std::cout << "heads\n";
    // heads only
    for(int i = 0; i < size; ++i)
    {
        a.heads[i] = (i + 1);// * 0.2;
        std::cout << a.heads[i] << " ";
    }
    // interim
    std::cout << "\ndoing interim assignments\n";
    for(int i = 0; i < size; ++i)
    {
        a.pairs[i] = a.heads[i] + 0.1; //(a.heads[i] * 1.25);
    }
    std::cout << "sum (heads)=" << sum<HeadsArray>(a.heads, size) << "\n";
    std::cout << "sum (pairs)=" << sum<PairsArray>(a.pairs, size) << "\n";
    // copy values
    a.copytoIEEEdouble();
    std::cout << "sum (full)=" << sum(a.full, size) << "\n";
    // use full
    std::cout << "regular doubles\n";
    for(int i = 0; i < size; ++i)
    {
        std::cout << a.full[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

int main_55()
{
    HeadsArray h1(5);
    HeadsArray h2(5);

    h1[0] = 0.11;
    h1[1] = 1.4f;
    h1[2] = 22;
    h1[3] = 45L;
    h1[4] = true;
    std::cout << "doing assignments now\n";
    h2[0] = h1[0];
    h2[1] = h1[1];
    h2[2] = h1[2];
    h2[3] = h1[3];
    h2[4] = h1[4];
    std::cout << "pairs time\n";
    PairsArray p1(5);
    PairsArray p2(5);

    p1[0] = 0.11;
    p1[1] = 1.4f;
    p1[2] = 22;
    p1[3] = 45L;
    p1[4] = true;
    std::cout << "doing assignments now\n";
    p2[0] = p1[0];
    p2[1] = p1[1];
    p2[2] = p1[2];
    p2[3] = p1[3];
    p2[4] = p1[4];

    
    printf("heads\n");
    for(int i = 0; i < 5; ++i)
    {
        printf(" %.2f = ", (double)h1[i]);
        printBinary((double)h1[i]);
        printf(" %.2f = ", (double)h2[i]);
        printBinary((double)h2[i]);
        printf("\n");
    }
    printf("\n");

    printf("pairs\n");
    for(int i = 0; i < 5; ++i)
    {
        printf(" %.2f = ", (double)p1[i]);
        printBinary((double)p1[i]);
        printf(" %.2f = ", (double)p2[i]);
        printBinary((double)p2[i]);
        printf("\n");
    }

    return 0;
}

int main01(int argc, char** argv)
{
    // std::cout << std::setprecision(16);
    // std::cout.setf(std::ios::fixed, std::ios::floatfield);
    const int size = 6;

    TwoSegArray<false> f(size);
    TwoSegArray<true> t(size);

    // f.setPair(0, 0.1);
    // f.setPair(1, 1.4f);
    // f.setPair(2, 22);
    // f.setPair(3, 45L);
    // f.setPair(4, true);

    f[0] = 0.11;
    f[1] = 1.4f;
    f[2] = 22;
    f[3] = 45L;
    f[4] = true;
    f[5] = 1.25e-308;
 
    // t[0] = f[0]; 
    // t[1] = f[1]; 
    // t[2] = f[2];
    // t[3] = f[3];
    // t[4] = f[4];
    // t[5] = f[5];
    t[0] = 0.11;
    t[1] = 1.4f;
    t[2] = 22;
    t[3] = 45L;
    t[4] = true;
    t[5] = 1.25e-308;

    for(int i = 0; i < size; ++i)
    {
        double d = f[i];
        std::cout << i << ": " << d << "\n";
        // printBinary(d);
        std::cout << "\n";
    }
    std::cout << "\n";
    for(int i = 0; i < size; ++i)
    {
        double d = t[i];
        std::cout << i << ": " << d << "\n";
        // printBinary(d);
        std::cout << "\n";
    }

    f.del();
    t.del();

    return 0;
}

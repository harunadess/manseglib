/*
* Copyright (c) 2008, BSC (Barcelon Supercomputing Center)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY BSC ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include "../../../../manseglib.hpp"

using namespace ManSeg;

#define NB 128
#define B 64
#define FALSE (0)
#define TRUE (1)

typedef double fp_type;

typedef fp_type* vin;
typedef fp_type* vout;

typedef ManSegArray bin;
typedef ManSegArray binout;

ManSegArray A[NB][NB];
ManSegArray A_new[NB][NB];

enum Precision { HEADS, PAIRS, INTERIM }; // INTERIM = read heads, write pairs

Precision MatrixPrecision = Precision::HEADS;

void alloc_and_genmat()
{
    int init_val, i, j, ii, jj;
    fp_type *p, *p_new;

    init_val = 1325;

    for (ii = 0; ii < NB; ii++)
    {
        for (jj = 0; jj < NB; jj++)
        {
            A[ii][jj].alloc(B * B);
			A[ii][jj].full = new double[B * B];
            A_new[ii][jj].alloc(B * B);
			A_new[ii][jj].full = new double[B * B];

            for (i = 0; i < B; i++)
            {
                for (j = 0; j < B; j++)
                {
                    init_val = (3125 * init_val) % 65536;
                    A[ii][jj].heads[i * B + j] = (fp_type)((init_val - 32768.0) / 16384.0);
                }
            }
        }
    }
}

long usecs(void)
{
    struct timeval t;

    gettimeofday(&t, NULL);
    return t.tv_sec * 1000000 + t.tv_usec;
}

void clear(vout v)
{
    int i, j, k;
    for (i = 0; i < B; i++)
        v[i] = (fp_type)0.0;
}

template<enum Precision = Precision::PAIRS>
void getlastrow(bin& A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A.full[(B - 1) * B + j];
}
template<>
void getlastrow<Precision::HEADS>(bin& A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A.heads[(B - 1) * B + j];
}

template<enum Precision = Precision::PAIRS>
void getlastcol(bin& A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A.full[i * B + B - 1];
}
template<>
void getlastcol<Precision::HEADS>(bin& A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A.heads[i * B + B - 1];
}

template<enum Precision = Precision::PAIRS>
void getfirstrow(bin& A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A.full[0 * B + j];
}
template<>
void getfirstrow<Precision::HEADS>(bin& A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A.heads[0 * B + j];
}

template<enum Precision = Precision::PAIRS>
void getfirstcol(bin& A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A.full[i * B + 0];
}
template<>
void getfirstcol<Precision::HEADS>(bin& A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A.heads[i * B + 0];
}

template<enum Precision = Precision::PAIRS>
void jacobi(vin lefthalo, vin tophalo, vin righthalo, vin bottomhalo, 
                bin& A, binout& A_new, int ii, int jj, int iter)
{
    int i, j;
    fp_type tmp;
    fp_type left, top, right, bottom;

    for (i = 0; (i < B); i++)
    {
        for (j = 0; j < B; j++)
        {
            tmp = A.full[i * B + j];
            left = (j == 0 ? lefthalo[j] : A.full[i * B + j - 1]);
            top = (i == 0 ? tophalo[i] : A.full[(i - 1) * B + j]);
            right = (j == B - 1 ? righthalo[i] : A.full[i * B + j + 1]);
            bottom = (i == B - 1 ? bottomhalo[i] : A.full[(i + 1) * B + j]);

            A_new.full[i * B + j] = 0.2 * (A.full[i * B + j] + left + top + right + bottom);
        }
    }
}
template<>
void jacobi<Precision::HEADS>(vin lefthalo, vin tophalo, vin righthalo, vin bottomhalo, 
                bin& A, binout& A_new, int ii, int jj, int iter)
{
    int i, j;
    fp_type tmp;
    fp_type left, top, right, bottom;

    for (i = 0; (i < B); i++)
    {
        for (j = 0; j < B; j++)
        {
            tmp = A.heads[i * B + j];
            left = (j == 0 ? lefthalo[j] : A.heads[i * B + j - 1]);
            top = (i == 0 ? tophalo[i] : A.heads[(i - 1) * B + j]);
            right = (j == B - 1 ? righthalo[i] : A.heads[i * B + j + 1]);
            bottom = (i == B - 1 ? bottomhalo[i] : A.heads[(i + 1) * B + j]);

            A_new.heads[i * B + j] = 0.2 * (A.heads[i * B + j] + left + top + right + bottom);
        }
    }
}

template<class Arr>
inline double blockdelta(Arr& A_new, Arr& A)
{
	double delta = -__DBL_MAX__;
	for (int i = 0; (i < B); ++i)
	{
		for (int j = 0; j < B; ++j)
		{
			double diff = fabs(A_new[i * B + j] - A[i * B + j]);
			if(diff > delta) delta = diff;
		}
	}
	return delta;
}

double maxdelta(int iters)
{
	double dmax = -__DBL_MAX__;

	#pragma omp parallel for schedule(static) reduction(max: dmax)
	for (int ii = 0; ii < NB; ++ii)
	{
		for (int jj = 0; jj < NB; ++jj)
		{
			// if(BlockPrecision[ii][jj] == Precision::HEADS)
			if(MatrixPrecision == Precision::HEADS)
			{
				double blockmax = blockdelta(A_new[ii][jj].heads, A[ii][jj].heads);
				if(dmax < blockmax) dmax = blockmax;

				/* // do copy from heads to full if block delta is small enough
				if(blockmax < AdaptivePrecisionBound)
				{
					printf("block[%d][%d] switch at iter %d\n", ii, jj, iters);
					BlockPrecision[ii][jj] = Precision::PAIRS;
					for (int i = 0; (i < B); ++i)
						for (int j = 0; j < B; ++j)
							A_new[ii][jj].full[i * B + j] = A_new[ii][jj].heads[i * B + j];
				} */
			}
			// else if(BlockPrecision[ii][jj] == Precision::PAIRS)
			else if(MatrixPrecision == Precision::PAIRS)
			{
				double blockmax = blockdelta(A_new[ii][jj].full, A[ii][jj].full);
				if(dmax < blockmax) dmax = blockmax;
			}
		}
	}

	return dmax;
}

void compute(int niters)
{
    int iters;
    int ii, jj;
    fp_type lefthalo[B], tophalo[B], righthalo[B], bottomhalo[B];

	double delta = 2.0;
	double epsilon = 1e-7;

	iters = 0;
    // for (iters = 0; iters < niters; iters++)
	// while((delta > epsilon) && (iters < niters))
	while(iters < niters)
    {
		++iters;

        #pragma omp parallel \
            private(iters, ii, jj, lefthalo, tophalo, righthalo, bottomhalo) \
            shared(A, A_new) 
        {
            #pragma omp for schedule(static)
            for (ii = 0; ii < NB; ii++)
            {
                for (jj = 0; jj < NB; jj++)
                {
                //    if(BlockPrecision[ii][jj] == Precision::HEADS)
                   if(MatrixPrecision == Precision::HEADS)
				   {
					    if (ii > 0)
                        	getlastrow<Precision::HEADS>(A[ii - 1][jj], tophalo);
						else
							clear(tophalo);

						if (jj > 0)
							getlastcol<Precision::HEADS>(A[ii][jj - 1], lefthalo);
						else
							clear(lefthalo);

						if (ii < NB - 1)
							getfirstrow<Precision::HEADS>(A[ii + 1][jj], bottomhalo);
						else
							clear(bottomhalo);

						if (jj < NB - 1)
							getfirstcol<Precision::HEADS>(A[ii][jj + 1], righthalo);
						else
							clear(lefthalo);

                        jacobi<Precision::HEADS>(lefthalo, tophalo, righthalo, bottomhalo, A[ii][jj], A_new[ii][jj], ii, jj, iters);
				   }
                    else
					{
						if (ii > 0)
                        	getlastrow<Precision::PAIRS>(A[ii - 1][jj], tophalo);
						else
							clear(tophalo);

						if (jj > 0)
							getlastcol<Precision::PAIRS>(A[ii][jj - 1], lefthalo);
						else
							clear(lefthalo);

						if (ii < NB - 1)
							getfirstrow<Precision::PAIRS>(A[ii + 1][jj], bottomhalo);
						else
							clear(bottomhalo);

						if (jj < NB - 1)
							getfirstcol<Precision::PAIRS>(A[ii][jj + 1], righthalo);
						else
							clear(lefthalo);

                        jacobi<Precision::PAIRS>(lefthalo, tophalo, righthalo, bottomhalo, A[ii][jj], A_new[ii][jj], ii, jj, iters);
					}
                } // jj
            } // ii
        } // end parallel

		delta = maxdelta(iters);
		printf("iteration %d: delta = %e\n", iters, delta);

		// yes, this is an inefficient copy
		// however, this is required to avoid segmentation fault
		// because of some issue with constructor/deconstructor and
		// memory allocation/deallocation i have not had time to
		// diagnose.
        // memcpy(A, A_new);
		/* #pragma omp parallel for schedule(static) shared(A, A_new)
		for(int i = 0; i < NB; ++i)
		{
			for(int j = 0; j < NB; ++j)
			{
				if(BlockPrecision[i][j] == Precision::HEADS)
				{
					for(int k = 0; k < B; ++k)
						for(int l = 0; l < B; ++l)
							A[i][j].heads[k * B + l] = A_new[i][j].heads[k * B + l];
				}
				else
				{
					for(int k = 0; k < B; ++k)
						for(int l = 0; l < B; ++l)
							A[i][j].full[k * B + l] = A_new[i][j].full[k * B + l];
				}
			}
		} */
		
		if(MatrixPrecision == Precision::HEADS)
		{
			// precision switch
			if(delta <= AdaptivePrecisionBound)
			{
				MatrixPrecision = Precision::PAIRS;
				printf("precision switch at iter %d\n", iters);
				#pragma omp parallel for schedule(static) shared(A, A_new)
				for(int i = 0; i < NB; ++i)
				{
					for(int j = 0; j < NB; ++j)
					{
						for(int k = 0; k < B; ++k)
							for(int l = 0; l < B; ++l)
								A[i][j].full[k * B + l] = A_new[i][j].heads[k * B + l];
					}
				}
			}
			else
			{
				#pragma omp parallel for schedule(static) shared(A, A_new)
				for(int i = 0; i < NB; ++i)
				{
					for(int j = 0; j < NB; ++j)
					{
						for(int k = 0; k < B; ++k)
							for(int l = 0; l < B; ++l)
								A[i][j].heads[k * B + l] = A_new[i][j].heads[k * B + l];
					}
				}
			}
		}
		else if(MatrixPrecision == Precision::PAIRS)
		{
			#pragma omp parallel for schedule(static) shared(A, A_new)
			for(int i = 0; i < NB; ++i)
			{
				for(int j = 0; j < NB; ++j)
				{
					for(int k = 0; k < B; ++k)
						for(int l = 0; l < B; ++l)
							A[i][j].full[k * B + l] = A_new[i][j].full[k * B + l];
				}
			}
		}
		
    } // iter

	if(iters >= niters)
		printf("hit max iters\n");
	if(delta <= epsilon)
		printf("converged to %e\n", delta); 
}

int main(int argc, char *argv[])
{
    int niters;
    struct timespec start, end;

    if (argc > 1)
    {
        niters = atoi(argv[1]);
    }
    else
        niters = 1;

    alloc_and_genmat();

    clock_gettime(CLOCK_MONOTONIC, &start);
    compute(niters);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("Running time  = %g %s\n", time_taken, "s");

    /* FILE *outFile;
    outFile = fopen("./jacobi_mod_omp_values.txt", "w");
    if (outFile == NULL)
    {
        fprintf(stderr, "Error writing to file\n");
    }
    else
    {
        int ii, jj, i, j;
        for (ii = 0; ii < NB; ++ii)
            for (jj = 0; jj < NB; ++jj)
                // if(BlockPrecision[ii][jj] == Precision::PAIRS)
                if(MatrixPrecision == Precision::PAIRS)
					for (i = 0; i < B; ++i)
						for (j = 0; j < B; ++j)
							fprintf(outFile, "%.15f\n", (double)(A[ii][jj].full[i * B + j]));
				else
					for (i = 0; i < B; ++i)
						for (j = 0; j < B; ++j)
							fprintf(outFile, "%.15f\n", (double)(A[ii][jj].heads[i * B + j]));

        fclose(outFile);
    } */

    return 0;
}
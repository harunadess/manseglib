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
#include <math.h>
#include <time.h>

#define NB 64
#define B 128
#define FALSE (0)
#define TRUE (1)

typedef double fp_type;

typedef fp_type *vin;
typedef fp_type *vout;

typedef fp_type *bin;
typedef fp_type *binout;

fp_type *A[NB][NB];
fp_type *A_new[NB][NB];
fp_type *tmp[NB][NB];

void alloc_and_genmat()
{
    int init_val, i, j, ii, jj;
    fp_type *p, *p_new;

    init_val = 1325;

    for (ii = 0; ii < NB; ii++)
    {
        for (jj = 0; jj < NB; jj++)
        {
            A[ii][jj] = (fp_type *)malloc(B * B * sizeof(fp_type));
            A_new[ii][jj] = (fp_type *)malloc(B * B * sizeof(fp_type));
            tmp[ii][jj] = (fp_type *)malloc(B * B * sizeof(fp_type));

            if (A[ii][jj] == NULL || A_new[ii][jj] == NULL || tmp[ii][jj] == NULL)
            {
                printf("Out of memory\n");
                exit(1);
            }
            p = A[ii][jj];
            p_new = A_new[ii][jj];
            for (i = 0; i < B; i++)
            {
                for (j = 0; j < B; j++)
                {
                    init_val = (3125 * init_val) % 65536;
                    (*p) = (fp_type)((init_val - 32768.0) / 16384.0);
                    (*p_new) = (*p);
                    p++;
                    p_new++;
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

void getlastrow(bin A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A[(B - 1) * B + j];
}

void getlastcol(bin A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A[i * B + B - 1];
}

void getfirstrow(bin A, vout v)
{
    int j;
    for (j = 0; j < B; j++)
        v[j] = A[0 * B + j];
}

void getfirstcol(bin A, vout v)
{
    int i;
    for (i = 0; i < B; i++)
        v[i] = A[i * B + 0];
}

void jacobi(vin lefthalo, vin tophalo, vin righthalo, vin bottomhalo, bin A, binout A_new)
{
    int i, j;
    fp_type tmp;
    fp_type left, top, right, bottom;

    for (i = 0; (i < B); i++)
    {
        for (j = 0; j < B; j++)
        {
            tmp = A[i * B + j];
            left = (j == 0 ? lefthalo[j] : A[i * B + j - 1]);
            top = (i == 0 ? tophalo[i] : A[(i - 1) * B + j]);
            right = (j == B - 1 ? righthalo[i] : A[i * B + j + 1]);
            bottom = (i == B - 1 ? bottomhalo[i] : A[(i + 1) * B + j]);

            A_new[i * B + j] = 0.2 * (A[i * B + j] + left + top + right + bottom);
        }
    }
}

double maxdelta()
{
	double dmax = -__DBL_MAX__;

	int ii, jj, i, j;
	#pragma omp parallel for schedule(static) reduction(max: dmax)
	for (ii = 0; ii < NB; ii++)
	{
		for (jj = 0; jj < NB; jj++)
		{
			for (i = 0; (i < B); i++)
			{
				for (j = 0; j < B; j++)
				{
					double diff = fabs(A_new[ii][jj][i * B + j] - A[ii][jj][i * B + j]);
					if(diff > dmax) dmax = diff;
				}
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
	while(iters < niters)
    {
		++iters;
        #pragma omp parallel \
            private(ii, jj, lefthalo, tophalo, righthalo, bottomhalo) \
            shared(A, A_new)
        {
            #pragma omp for schedule(static)
            for (ii = 0; ii < NB; ii++)
            {
                for (jj = 0; jj < NB; jj++)
                {
                    if (ii > 0)
                        getlastrow(A[ii - 1][jj], tophalo);
                    else
                        clear(tophalo);

                    if (jj > 0)
                        getlastcol(A[ii][jj - 1], lefthalo);
                    else
                        clear(lefthalo);

                    if (ii < NB - 1)
                        getfirstrow(A[ii + 1][jj], bottomhalo);
                    else
                        clear(bottomhalo);

                    if (jj < NB - 1)
                        getfirstcol(A[ii][jj + 1], righthalo);
                    else
                        clear(lefthalo);

                    jacobi(lefthalo, tophalo, righthalo, bottomhalo, A[ii][jj], A_new[ii][jj]);
                } // jj
            } // ii
        } // end parallel

		delta = maxdelta();
		printf("iteration %d: delta = %e\n", iters, delta);

		// yes, this is an inefficient copy
		// however, the library version requires you to do a copy in this way
		// on all of the component parts to avoid segmentation fault
		#pragma omp parallel for schedule(static) shared(A, A_new)
		for(int i = 0; i < NB; ++i)
		{
			for(int j = 0; j < NB; ++j)
			{
				for(int k = 0; k < B; ++k)
					for(int l = 0; l < B; ++l)
						A[i][j][k * B + l] = A_new[i][j][k * B + l];
			}
		}
			
    } // iter
}

int main(int argc, char *argv[])
{
    int niters;
    //  pp_time_t tm;
    //  memset( &tm, 0, sizeof(tm) );
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
    outFile = fopen("./jacobi_omp_values.txt", "w");
    if (outFile == NULL)
    {
        fprintf(stderr, "Error writing to file\n");
    }
    else
    {
        int ii, jj, i, j;
        for (ii = 0; ii < NB; ++ii)
            for (jj = 0; jj < NB; ++jj)
                for (i = 0; i < B; ++i)
                    for (j = 0; j < B; ++j)
                        fprintf(outFile, "%.15f\n", A[ii][jj][i * B + j]);

        fclose(outFile);
    } */

    return 0;
}
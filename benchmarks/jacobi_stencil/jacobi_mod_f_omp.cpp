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
#include "../../manseglib_vector_dev.hpp"

using namespace ManSeg;

#define NB 128
#define B 32
#define FALSE (0)
#define TRUE (1)

typedef double fp_type;

typedef fp_type* vin;
typedef fp_type* vout;

typedef ManSegArray bin;
typedef ManSegArray binout;

ManSegArray A[NB][NB];
ManSegArray A_new[NB][NB];
fp_type blockDelta[NB][NB];

enum Precision { HEADS, PAIRS, INTERIM }; // INTERIM = read heads, write pairs

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
            A_new[ii][jj].alloc(B * B);

            blockDelta[ii][jj] = 0.0;
            for (i = 0; i < B; i++)
            {
                for (j = 0; j < B; j++)
                {
                    init_val = (3125 * init_val) % 65536;
                    A[ii][jj].pairs[i * B + j] = (fp_type)((init_val - 32768.0) / 16384.0);
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

template<enum Precision = Precision::PAIRS>
void clear(vout v)
{
    int i, j, k;
    for (i = 0; i < B; i++)
        v[i] = (fp_type)0.0;
}
template<>
void clear<Precision::HEADS>(vout v)
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
        v[j] = A.pairs[(B - 1) * B + j];
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
        v[i] = A.pairs[i * B + B - 1];
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
        v[j] = A.pairs[0 * B + j];
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
        v[i] = A.pairs[i * B + 0];
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
                bin& A, binout& A_new, fp_type& blockDelta)
{
    int i, j;
    fp_type tmp;
    fp_type left, top, right, bottom;
    fp_type fullResult, deltaErr = 0.0;

    for (i = 0; (i < B); i++)
    {
        for (j = 0; j < B; j++)
        {
            tmp = A.pairs[i * B + j];
            left = (j == 0 ? lefthalo[j] : A.pairs[i * B + j - 1]);
            top = (i == 0 ? tophalo[i] : A.pairs[(i - 1) * B + j]);
            right = (j == B - 1 ? righthalo[i] : A.pairs[i * B + j + 1]);
            bottom = (i == B - 1 ? bottomhalo[i] : A.pairs[(i + 1) * B + j]);

            // A_new.pairs[i * B + j] = 0.2 * (A.pairs[i * B + j] + left + top + right + bottom);
            
            fullResult = 0.2 * (A.pairs[i * B + j] + left + top + right + bottom);
            A_new.pairs[i * B + j] = fullResult;

            // record difference between full result and stored value
            // blockDelta += std::fabs(fullResult - A_new.pairs[i * B + j])
            fp_type deltaTmp = blockDelta;
            fp_type y = std::fabs(fullResult - A_new.pairs[i * B + j]) + deltaErr;
            
            blockDelta = deltaTmp + y;
            deltaErr = deltaTmp - blockDelta;
            deltaErr += y;
        }
    }
}
template<>
void jacobi<Precision::HEADS>(vin lefthalo, vin tophalo, vin righthalo, vin bottomhalo, 
                bin& A, binout& A_new, fp_type& blockDelta)
{
    int i, j;
    fp_type tmp;
    fp_type left, top, right, bottom;
    fp_type fullResult, deltaErr = 0.0;

    for (i = 0; (i < B); i++)
    {
        for (j = 0; j < B; j++)
        {
            tmp = A.heads[i * B + j];
            left = (j == 0 ? lefthalo[j] : A.heads[i * B + j - 1]);
            top = (i == 0 ? tophalo[i] : A.heads[(i - 1) * B + j]);
            right = (j == B - 1 ? righthalo[i] : A.heads[i * B + j + 1]);
            bottom = (i == B - 1 ? bottomhalo[i] : A.heads[(i + 1) * B + j]);

            // A_new.pairs[i * B + j] = 0.2 * (A.pairs[i * B + j] + left + top + right + bottom);
            
            fullResult = 0.2 * (A.heads[i * B + j] + left + top + right + bottom);
            A_new.heads[i * B + j] = fullResult;

            // record difference between full result and stored value
            // blockDelta += std::fabs(fullResult - A_new.heads[i * B + j])
            fp_type deltaTmp = blockDelta;
            fp_type y = std::fabs(fullResult - A_new.heads[i * B + j]) + deltaErr;
            
            blockDelta = deltaTmp + y;
            deltaErr = deltaTmp - blockDelta;
            deltaErr += y;
        }
    }
}

void compute(int niters)
{
    int iters;
    int ii, jj;
    fp_type lefthalo[B], tophalo[B], righthalo[B], bottomhalo[B];

    for (iters = 0; iters < niters; iters++)
    {
        #pragma omp parallel for schedule(static) \
            private(ii, jj, lefthalo, tophalo, righthalo, bottomhalo) shared(A, A_new, blockDelta)
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

                if(blockDelta[ii][jj] > MaxSingleSegmentPrecision) // suggested delta
                    jacobi<Precision::PAIRS>(lefthalo, tophalo, righthalo, bottomhalo, A[ii][jj], A_new[ii][jj], blockDelta[ii][jj]);
                else
                    jacobi<Precision::HEADS>(lefthalo, tophalo, righthalo, bottomhalo, A[ii][jj], A_new[ii][jj], blockDelta[ii][jj]);
            } // jj
        } // ii

        // std::copy(A, A_new);
        for(int i = 0; i < NB; ++i)
            for(int j = 0; j < NB; ++j)     // for these inner two loops we could split out into a separate funciton
                for(int k = 0; k < B; ++k)  // if blockDelta[i][j] > MaxSingle, then do pairs (or full) 
                    for(int l = 0; l < B; ++l)
                        A[i][j].pairs[k * B + l] = A_new[i][j].pairs[k * B + l];
    } // iter
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

    FILE *outFile;
    outFile = fopen("./jacobi_mod_f_omp_values.txt", "w");
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
                        fprintf(outFile, "%.15f\n", (double)(A[ii][jj].pairs[i * B + j]));

        fclose(outFile);
    }
    
    // write out block deltas
    outFile = fopen("./jacobi_mod_f_omp_deltas.txt", "w");
    if (outFile == NULL)
    {
        fprintf(stderr, "Error writing to file\n");
    }
    else
    {
        int ii, jj, i, j;
        for (ii = 0; ii < NB; ++ii)
            for (jj = 0; jj < NB; ++jj)
                fprintf(outFile, "%.15f\n", blockDelta[ii][jj]);

        fclose(outFile);
    }

    return 0;
}

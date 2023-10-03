
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#include "mmio.h"
#include <tgmath.h> // interferes with mmio.h

#include "cg.h"
#include "matrix.h"

static int compar(const void *pa, const void *pb)
{
    matrix_coo *a = (matrix_coo*)pa;
    matrix_coo *b = (matrix_coo*)pb;
    if (a->i < b->i) return -1;
    if (a->i > b->i) return 1;
    if (a->j < b->j) return -1;
    if (a->j > b->j) return 1;
    return 0;
}

matrix_coo* coo_load(const char *fname, int *n, int *nz)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error opening file: %s\n", fname);
        exit(1);
    }
    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner\n");
        exit(1);
    }
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || mm_is_complex(matcode) || !mm_is_symmetric(matcode)) {
        fprintf(stderr, "This application does not support the Market Market type: %s\n",
                mm_typecode_to_str(matcode));
        exit(1);
    }
    int M, N, NZ;
    if (mm_read_mtx_crd_size(f, &M, &N, &NZ) != 0) {
        fprintf(stderr, "Could not parse matrix size\n");
        exit(1);
    }
    if (M != N) {
        fprintf(stderr, "Matrix is not square\n");
        exit(1);
    }
    matrix_coo *coo = ALLOC(matrix_coo, 2 * NZ);

    int k = 0;
    for (int l = 0; l < NZ; l++) {
        double real;
        if (mm_read_mtx_crd_entry(f, &coo[k].i, &coo[k].j, &real, NULL, matcode) != 0) {
            fprintf(stderr, "Error reading matrix element %i\n", l);
            exit(1);
        }
        coo[k].i--;
        coo[k].j--;
        coo[k].a = real;
        if (coo[k].i == coo[k].j) k++;
        else {
            coo[k + 1].i = coo[k].j;
            coo[k + 1].j = coo[k].i;
            coo[k + 1].a = coo[k].a;
            k += 2;
        }
    }
    fclose(f);
    qsort(coo, k, sizeof(matrix_coo), compar);
    *n = N; *nz = k;

    return coo;
}

double coo_norm_inf(int n, int nz, matrix_coo *coo)
{
    double *amax = CALLOC(double, n);
    for (int i = 0; i < nz; i++) {
        int row = coo[i].i;
        int val = coo[i].a;
        amax[row] += abs(val);
    }
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        if (amax[i] > norm) norm = amax[i];
    }
    FREE(amax);
    return norm;
}

double coo_max_nz(int n, int nz, matrix_coo *coo)
{
    double *m = CALLOC(double, n);
    for (int i = 0; i < nz; i++) {
        int row = coo[i].i;
        m[row]++;
    }
    int r = 0.0;
    for (int i = 0; i < n; i++) {
        if (m[i] > r) r = m[i];
    }
    FREE(m);
    return r;
}

// CSR matrix
void csr_dmult(matrix_csr *mat, DOUBLE *x, DOUBLE *y)
{
	#pragma omp parallel for \
		shared(mat, x, y)
    for (int k = 0; k < mat->n; k++) {
        DOUBLE t = 0.0;
        for (int l = mat->i[k]; l < mat->i[k + 1]; l++)
            t += mat->A[l] * x[mat->j[l]];
        y[k] = t;
    }
}

void csr_smult(matrix_csr *mat, FLOAT *x, FLOAT *y)
{
	#pragma omp parallel for \
		shared(mat, x, y)
    for (int k = 0; k < mat->n; k++) {
        FLOAT2 t = 0.0;
        for (int l = mat->i[k]; l < mat->i[k + 1]; l++)
            t += mat->A[l] * x[mat->j[l]];
        y[k] = t;
    }
}

matrix *csr_create(int n, int nz, matrix_coo *coo)
{
    int *i = ALLOC(int, n + 1);
    int *j = ALLOC(int, nz);
    DOUBLE *A = ALLOC(DOUBLE, nz);

    i[0] = 0;
    int l = 0;
    for (int k = 0; k < n; k++) {
        while (l < nz && coo[l].i == k) {
            j[l] = coo[l].j;
            A[l] = coo[l].a;
            l++;
        }
        i[k + 1] = l;
    }

    matrix_csr *mat = new matrix_csr();
    mat->n = n;
    mat->i = i;
    mat->j = j;
    mat->A = A;
    mat->dmult = (void (*)(matrix *, DOUBLE *, DOUBLE *))csr_dmult;
    mat->smult = (void (*)(matrix *, FLOAT *, FLOAT *))csr_smult;

    return (matrix *)mat;
}

// dense matrix
void dense_dmult(matrix_dense *mat, DOUBLE *x, DOUBLE *y)
{
	#pragma omp parallel for \
		shared(mat, x, y)
    for (int i = 0; i < mat->n; i++) {
        DOUBLE t = 0.0;
        for (int j = 0; j < mat->n; j++)
            t += mat->A[i * mat->n + j] * x[j];
        y[i] = t;
    }
}

void dense_smult(uint8_t m, matrix_dense *mat, FLOAT *x, FLOAT *y)
{
	#pragma omp parallel for \
		shared(mat, x, y)
    for (int i = 0; i < mat->n; i++) {
        FLOAT2 t = 0.0;
        for (int j = 0; j < mat->n; j++)
            t += mat->A[i * mat->n + j] * x[j];
        y[i] = t;
    }
}

matrix *dense_create(int n, int nz, matrix_coo *coo)
{
    DOUBLE *A = CALLOC(DOUBLE, n * n);

    // row major format
    for (int k = 0; k < nz; k++) A[coo[k].i * n + coo[k].j] = coo[k].a;
    
    matrix_dense *mat = new matrix_dense();
    mat->n = n;
    mat->dmult = (void (*)(matrix *, DOUBLE *, DOUBLE *))dense_dmult;
    mat->smult = (void (*)(matrix *, FLOAT *, FLOAT *))dense_smult;
    mat->A = A;

    return (matrix *)mat;
}

// Jacobi preconditioner
void jacobi_smult(precond_jacobi *pre, FLOAT *x, FLOAT *y)
{
	#pragma omp parallel for \
		shared(y, x, pre)
    for (int k = 0; k < pre->n; k++) {
       y[k] = x[k] / pre->d[k];
    }
}

struct matrix *jacobi_create(int n, int nz, matrix_coo *coo)
{
    FLOAT *d = ALLOC(FLOAT, n);
	// ManSegArray *d = new ManSegArray(n);
    for (int k = 0; k < nz; k++)
        if (coo[k].i == coo[k].j)
		{
            d[coo[k].i] = coo[k].a;
		}

    precond_jacobi *pre = new precond_jacobi();
    pre->n = n;
    pre->dmult = NULL;
    pre->smult = (void (*)(matrix *, FLOAT *, FLOAT *))jacobi_smult;
    pre->d = d;

    return (matrix *)pre;
}

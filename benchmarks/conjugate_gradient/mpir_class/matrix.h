#ifndef __cg_MATRIX_H__
#define __cg_MATRIX_H__

// definitions for matrices
class matrix
{
public:
    int n;
    
    void (*dmult)(matrix *, DOUBLE *, DOUBLE *);
    void (*smult)(matrix *, FLOAT *, FLOAT *);
};

class matrix_coo : public matrix
{
public:
    int i, j;
    double a;
};

class matrix_csr : public matrix
{
public:
    int n;
    int *i, *j;
    DOUBLE *A;
};

class matrix_dense : public matrix
{
public:
    int n;
    DOUBLE *A;
};

class precond_jacobi : public matrix
{
public:
    int n;
    FLOAT *d;
};

static inline void matrix_mult(matrix *mat, DOUBLE *x, DOUBLE *y) {
    mat->dmult(mat, x, y);
}

static inline void floatm_mult(matrix *mat, FLOAT *x, FLOAT *y) {
    mat->smult(mat, x, y);
}

extern matrix_coo *coo_load(const char *fname, int *n, int *nz);
extern double coo_norm_inf(int n, int nz, matrix_coo *coo);
extern double coo_max_nz(int n, int nz, matrix_coo *coo);

extern matrix *csr_create(int n, int nz, matrix_coo *coo);
extern matrix *dense_create(int n, int nz, matrix_coo *coo);
extern matrix *jacobi_create(int n, int nz, matrix_coo *coo);


#endif

#ifndef __cg_MATRIX_H__
#define __cg_MATRIX_H__

using namespace ManSeg;

// definitions for matrices
class matrix
{
public:
    int n;
    bool useTail;
    
    void (*dmult)(matrix *, DOUBLE *, DOUBLE *);
    void (*smult)(matrix *, FLOAT *, FLOAT *);

    void (*precision_increase)(matrix *);
	void (*precision_reduce)(matrix *);
};

class matrix_coo : public matrix
{
public:
    int i, j;
    double a;
    // ManSegArray *a;
};

class matrix_csr : public matrix
{
public:
    int n;
    int *i, *j;
    // double *A;
    ManSegArray *A;
};

class matrix_dense : public matrix
{
public:
    int n;
    // double *A;
    ManSegArray *A;
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

static inline void mat_increase_precision(matrix* mat) {
    mat->useTail = true;
    mat->precision_increase(mat);
}

static inline void mat_reduce_precision(matrix* mat) {
    mat->useTail = false;
    mat->precision_reduce(mat);
}

extern matrix_coo *coo_load(const char *fname, int *n, int *nz);
extern double coo_norm_inf(int n, int nz, matrix_coo *coo);
extern double coo_max_nz(int n, int nz, matrix_coo *coo);

extern matrix *csr_create(int n, int nz, matrix_coo *coo);
extern matrix *dense_create(int n, int nz, matrix_coo *coo);
extern matrix *jacobi_create(int n, int nz, matrix_coo *coo);


#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter, bool *precision_increase);

void iterative_refinement(int n, matrix *A, matrix *M, DOUBLE *b, DOUBLE *b_dash, DOUBLE *x, int out_maxiter, DOUBLE out_tol, 
    int in_maxiter, DOUBLE in_tol, int step_check, int *out_iter, int *in_iter)
{
    DOUBLE *e = ALLOC(DOUBLE, n);
    FLOAT *r = ALLOC(FLOAT, n);
    FLOAT *d = ALLOC(FLOAT, n);

    mixed_copy(n, x, d); // d = (float)x
    vector_set(n, 0.0, x);
	// fill r with low precision version first
    mixed_copy(n, b_dash, r); // r = (float)b_dash

	double x_norm, x_norm_prev;
    DOUBLE residual;
	bool precision_increase = false;
	bool switched = false;
    do
    {
		conjugate_gradient(n, A, M, r, d, in_maxiter, in_tol, step_check, in_iter, &precision_increase);
        
		mixed_axpy(n, 1.0, d, x); // x = x + d
		matrix_mult(A, x, e);

		// i *think* this is helping things, but it requires more testing
		// it would maybe make sense, since b is the other component of Ax = b [or AMx = Mb or whatever]
		// so having it be more accurate would (when we switch) would help hone in on the result
		if(!precision_increase)
			vector_xpby(n, b_dash, -1.0, e); // r = b - Ax
		else
			vector_xpby(n, b, -1.0, e); // r = b - Ax

        residual = vector_norm2(n, e);
        mixed_copy(n, e, r);
        floatm_set(n, 0.0, d);

		/* if(!switched && precision_increase)
		{
			mat_increase_precision(A);
			switched = true;
		} */

		printf("%d: outer: residual = %e\n", *out_iter, residual);
        // *energy += iter * (bits + 12) / 8;
        // printf("%d %d %d %e %e\n", *in_iter, 0, bits, (double)residual, (double)residual);

        (*out_iter)++;		
    } while ((residual > out_tol) && (*out_iter < out_maxiter));

    if (residual <= out_tol)
        printf("\n==== residual less than out_tol ====\n");
	if(*out_iter >= out_maxiter)
		printf("out_iter >= out_maxiter\n");

    FREE(e);
    FREE(r);
    FREE(d);
}
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <float.h>
#include <stdbool.h>

#include "cg.h"
#include "vector.h"
#include "matrix.h"

void conjugate_gradient(int n, matrix *A, matrix *M, FLOAT *b, FLOAT *x, int maxiter, FLOAT umbral, int step_check, int *in_iter)
{
    int iter = 0;
    FLOAT2 alpha, beta;
	FLOAT2 rho; // rho = r(k-1)^Tz(k-1)
	FLOAT2 tau; // tau = r(k)^Tz(k)
	FLOAT2 tol; // tol = ||r(k)||v2

	// x(0) = [0, 0, .., 0]

    FLOAT *r = ALLOC(FLOAT, n);
    FLOAT *p = ALLOC(FLOAT, n); // this may be the actual "conjugate gradient" (conjugate vectors)
    FLOAT *z = ALLOC(FLOAT, n);

    FLOAT *tr = ALLOC(FLOAT, n);

	// r(0) = b - Ax(0)
    floatm_mult(A, x, r);
    floatm_xpby(n, b, -1.0, r);

    if (M) {
		// Mz(0) = r(0)
        floatm_mult(M, r, z);
		// calculate rho
        rho = floatm_dot(n, r, z);
		// calculate tol
		tol = floatm_norm2(n, r);
    } else { // small optimization
		// z(0) = r(0)
        floatm_copy(n, r, z);
		// rho = r(0)^Tr(0)
        rho = floatm_dot(n, r, r);
        // tol = sqrt(rho)
		tol = sqrt(rho);
    }
	// p(1) = z(0)
	floatm_copy(n, z, p);

    int step = 0;
    FLOAT2 residual = tol;

	FLOAT *x_prev = CALLOC(FLOAT, n);

	// tol = recurrence residual
	// residual = true residual
	// 10 = deviation (ie. over 10x) from the true residual
    while ((iter < maxiter) && (tol > umbral)) {
		// w = Ap(k+1)
		floatm_mult(A, p, z);

		if (step < step_check) {
			step++;
		}
		else {
			floatm_mult(A, x, tr);		 // tr = Ax
			floatm_xpby(n, b, -1.0, tr); // tr = b - Ax
			residual = floatm_norm2(n, tr);
			printf("# rescheck: total_cg_iter=%d current_iter=%d tol=%e resid=%e\n", *in_iter, iter, (double)tol, (double)residual);

			double explicit_residual_deviation = residual/tol;

			float max_diff = floatm_max_diff_and_copy(n, x, x_prev);
			printf("max diff = %e\n", max_diff);

			if (explicit_residual_deviation > 10) {
				printf("broke out : iter = %d\n", *in_iter);
				break;
			}

			step = 1;
		}

		// alpha(k+1) = rho / p(k+1)w
		alpha = rho / floatm_dot(n, z, p);

		// x(k+1) = x + alpha(k+1)p(k+1)
		floatm_axpy(n, alpha, p, x);

		// r(k+1) = r(k) - alpha(k+1)(w)
		floatm_axpy(n, -alpha, z, r);

		// apply preconditioner
		if (M) {
			// Solve Mz(k) = r(k)
		    floatm_mult(M, r, z);
			// calculate tau
		    tau = floatm_dot(n, r, z);
			// calculate tol
		    tol = floatm_norm2(n, r);
		} else {
			// tau = inner product(r)
		    tau = floatm_dot(n, r, r);
			// tol = sqrt(tau)
		    tol = sqrt(tau);
		}
		
		// beta = tau / rho
		beta =  tau / rho;
		// update r(k-1)^Tz(k-1) to r(k)^Tz(k)
		rho = tau;
		
		// p(k+1) = z(k) + beta(k) + p(k)
		if (M) {
			floatm_xpby(n, z, beta, p);
		}
		else {
			floatm_xpby(n, r, beta, p);
		}

		iter++;
		(*in_iter)++;
	}

	if(iter >= maxiter)
		printf("======= iter > maxiter =======\n");
	if(tol <= umbral)
		printf("======= tol <= umbral - > %e =======\n", tol);

    FREE(r);
    FREE(p);
    FREE(z);
    FREE(tr);
	FREE(x_prev);
}
// BLAS like functions
#ifndef __cg_VECTOR_H__
#define __cg_VECTOR_H__

static inline void vector_set(int n, DOUBLE a, DOUBLE *x) {
    for (int i = 0; i < n; i++) x[i] = a;
}

static inline void vector_rand(int n, DOUBLE *x) {
    for (int i = 0; i < n; i++) x[i] = rand() / (double)RAND_MAX;
}

static inline void vector_copy(int n, DOUBLE *x, DOUBLE *y) {
    for (int i = 0; i < n; i++) y[i] = x[i];
}

static inline void vector_axpy(int n, DOUBLE a, DOUBLE *x, DOUBLE *y) {
    for (int i = 0; i < n; i++) y[i] = y[i] + x[i] * a;
}

static inline void vector_xpby(int n, DOUBLE *x, DOUBLE b, DOUBLE *y) {
    for (int i = 0; i < n; i++) y[i] = x[i] + b * y[i];
}

// i.e. inner product
static inline DOUBLE vector_dot(int n, DOUBLE *x, DOUBLE *y) {
    DOUBLE r = 0.0;
    for (int i = 0; i < n; i++) r += y[i] * x[i];
    return r;
}

static inline DOUBLE vector_norm2(int n, DOUBLE *x) {
    return sqrt(vector_dot(n, x, x));
}

// limited precision version

static inline void floatm_set(int n, FLOAT a, FLOAT *x) {
    FLOAT b = a;
    for (int i = 0; i < n; i++) x[i] = b;
}

static inline void floatm_rand(int n, FLOAT *x) {
    for (int i = 0; i < n; i++) x[i] = rand() / (double)RAND_MAX;
}

// y = x
static inline void floatm_copy(int n, FLOAT *x, FLOAT *y) {
    for (int i = 0; i < n; i++) y[i] = x[i];
}

// y = y + x*a
static inline void floatm_axpy(int n, FLOAT a, FLOAT *x, FLOAT *y) {
    for (int i = 0; i < n; i++) y[i] = y[i] + x[i] * a;
}

// y = x + b*y
static inline void floatm_xpby(int n, FLOAT *x, FLOAT b, FLOAT *y) {
    for (int i = 0; i < n; i++) y[i] = x[i] + b * y[i];
}

static inline FLOAT floatm_max_diff_and_copy(int n, FLOAT *x, FLOAT *y) {
	float maxv = -__FLT_MIN__;
	#pragma omp parallel for reduction(max: maxv)
	for(int i = 0; i < n; i++) {
		float val = fabsf(x[i] - y[i]);
		y[i] = x[i];
		
		if(val > maxv)
			maxv = val;
	}

	return maxv;
}

// the dot product can be computed with more precision than FLOAT

// inner product
static inline FLOAT2 floatm_dot(int n, FLOAT *x, FLOAT *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) r += y[i] * x[i];
    return r;
}

// sqrt(inner product(x))
static inline FLOAT2 floatm_norm2(int n, FLOAT *x) {
    return sqrt(floatm_dot(n, x, x));
}

// sum(y*(xa))
static inline FLOAT2 floatm_axpy_dot(int n, FLOAT a, FLOAT *x, FLOAT *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) {
        FLOAT2 t = y[i] + x[i] * a;
        y[i] = t;
        r += t * t;
    }
    return r;
}

// sqrt( sum((x-y)^2) )
static inline FLOAT2 floatm_diff_norm2(int n, FLOAT *x, FLOAT *y) {
    FLOAT2 r = 0.0;
    for (int i = 0; i < n; i++) {
        FLOAT2 t = x[i] - y[i];
        r += t * t;
    }
    return sqrt(r);
}

// z = max(abs(x - y), tol)
static inline void floatm_max_abs_diff(int n, FLOAT *x, FLOAT *y, FLOAT2* z, FLOAT2 tol) {
	for(int i = 0; i < n; i++) {
		z[i] = std::max(fabs(((double)x[i] - (double)y[i])), tol);
	}
}

static inline FLOAT floatm_max_diff(int n, FLOAT *x, FLOAT *y) {
	float max = -__FLT_MIN__;
	for(int i = 0; i < n; i++) {
		float val = fabsf(x[i] - y[i]);
		
		if(val > max)
			max = val;
	}

	return max;
}

// z = y ./ x
static inline void floatm_ratio(int n, FLOAT *x, FLOAT *y, FLOAT2 *z) {
	for(int i = 0; i < n; i++) z[i] = (y[i] / x[i]);
}

// z = abs(x - y)
static inline void floatm_abs_diff(int n, FLOAT *x, FLOAT *y, FLOAT *z) {
	for(int i = 0; i < n; i++) z[i] = fabsf(x[i] - y[i]);
}

static inline DOUBLE double_norm_diff(int n, DOUBLE *x, DOUBLE *y) {
    DOUBLE d = 0.0;
    DOUBLE err = 0.0;
    for(int i = 0; i < n; i++) {
        DOUBLE temp = d;
        DOUBLE r = std::fabs(y[i] - x[i]) + err;
        d = temp + r;
        err = temp - d;
        err += r;
    }
    return d;
}

static inline FLOAT2 float_norm_diff(int n, FLOAT *x, FLOAT *y) {
    FLOAT2 d = 0.0;
    FLOAT2 err = 0.0;
    for(int i = 0; i < n; i++) {
        FLOAT2 temp = d;
        FLOAT2 r = std::fabs(y[i] - x[i]) + err;
        d = temp + r;
        err = temp - d;
        err += r;
    }
    return d;
}

// mixed precision routines

static inline void mixed_copy(int n, DOUBLE *x, FLOAT *y) {
    for (int i = 0; i < n; i++) y[i] = x[i];
}

static inline void mixed_axpy(int n, DOUBLE a, FLOAT *x, DOUBLE *y) {
    for (int i = 0; i < n; i++) y[i] += x[i] * a;
}

static inline void mixed_xpby(int n, DOUBLE *x, DOUBLE b, FLOAT *y) {
    for (int i = 0; i < n; i++) y[i] = x[i] + b * y[i];
}

#endif

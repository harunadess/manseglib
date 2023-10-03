#include "matrix.h"

void matr_::fill(double* x)
{
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            v[i].pairs[j] = x[i + j];
}

void matr_::mm(double* x, double* y)
{
    for(int i = 0; i < n; ++i)
    {
        double t = 0.;
        for(int j = 0; j < n; ++j)
            t += v[i].pairs[j] * x[j];
        y[i] = t;
    }
}

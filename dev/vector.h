#ifndef __VEC_H__
#define __VEC_H__

#include "manseglib.hpp"
#include "matrix.h"

double prod(matr_* x, int n)
{
    double p = 0.;

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            p += (x->v[i].pairs[j] * x->v[i].pairs[j]);

    return p;
}

#endif

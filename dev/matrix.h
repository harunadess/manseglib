#ifndef __MATR_H__
#define __MATR_H__

#include "manseglib.hpp"

class matr_
{
public:
    matr_(int n)
        : n(n)
    {
        v = new ManSeg::ManSegArray[n];
        for(int i = 0; i < n; ++i)
            v[i].alloc(n);
    }

    ManSeg::ManSegArray* v;

    void fill(double* x);
    void mm(double* x, double* y);

private:
    int n;
};

#endif

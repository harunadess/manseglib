#pragma once

template<typename T>
int partition(T* a, int lo, int hi)
{
    int pivot = a[(hi+lo) / 2];
    int i = lo-1;
    int j = hi+1;

    while(true)
    {
        while(a[++i] < pivot);
        while(a[--j] > pivot);
        if(i >= j) return j;
        std::swap(a[i], a[j]);
    }
}

template<typename T>
void quicksort(T* a, int lo, int hi)
{
    if(lo < hi)
    {
        int p = partition(a, lo, hi);
        quicksort(a, lo, p);
        quicksort(a, (p+1), hi);
    }
}

template<typename T>
int partition_pair(T* a, T* b, int lo, int hi)
{
    int pivot = b[(hi+lo) / 2];
    int i = lo-1;
    int j = hi+1;

    while(true)
    {
        while(b[++i] < pivot);
        while(b[--j] > pivot);
        if(i >= j) return j;
        std::swap(b[i], b[j]);
        std::swap(a[i], a[j]);
    }
}

// implementation of quicksort that sorts values based on primary array
// using C.A.R Hoare method
template<typename T>
void quicksort_pair(T* secondary, T* primary, int lo, int hi)
{
    if(lo < hi)
    {
        int p = partition_pair(secondary, primary, lo, hi);
        quicksort_pair(secondary, primary, lo, p);
        quicksort_pair(secondary, primary, (p+1), hi);
    }
}

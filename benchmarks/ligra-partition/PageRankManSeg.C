// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra-numa.h"
#include "math.h"
#include "../../manseglib.hpp"
using namespace ManSeg;
int MaxIter=100;
template<class vertex, class ArrayType>
struct PR_F
{
    ArrayType& p_curr;
    ArrayType& p_next;
    double damping;
    vertex* V;
    static const bool use_cache = true;

    struct cache_t
    {
        double p_next;
    };
    PR_F(ArrayType& _p_curr, ArrayType& _p_next, double _damping, vertex* _V) :
        p_curr(_p_curr), p_next(_p_next), damping(_damping), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next[d] += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
       /*  writeAdd(&p_next[d], (damping*p_curr[s]/V[s].getOutDegree()));
		
        does not work
		volatile unsigned int oh, ot, nh, nt;
        volatile double next;
        double c = p_curr[s]/V[s].getOutDegree();
        do
        {
            next = (p_next[d] + c);
            oh = p_next[d].head;
            ot = p_next[d].tail;
            volatile std::uint64_t nl = *reinterpret_cast<volatile std::uint64_t*>(&next);
            nh = (nl >> 32U);
            nt = (nl & 0xFFFFFFFF);
        } while((!CAS(&p_next[d].head, oh, nh) || !CAS(&p_next[d].tail, ot, nt)) && (p_next[d] != next)); */

        cerr << "updateAtomic not implemented for ManSeg type.. yet\n";
        abort();
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        // Cache is used only in sequential mode
        p_next[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

template <class vertex>
struct PR_F_d
{
    double* p_curr, *p_next;
    double damping;
    vertex* V;
    static const bool use_cache = true;

    struct cache_t
    {
        double p_next;
    };
    PR_F_d(double* _p_curr, double* _p_next, double _damping, vertex* _V) :
        p_curr(_p_curr), p_next(_p_next), damping(_damping), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next[d] += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        writeAdd(&p_next[d], damping*(p_curr[s]/V[s].getOutDegree()));
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        // Cache is used only in sequential mode
        p_next[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

template<class vertex>
struct PR_Interim_F
{
    HeadsArray& p_curr;
    double* p_next;
    double damping;
    vertex* V;
    static const bool use_cache = true;

    struct cache_t
    {
        double p_next;
    };
    PR_Interim_F(HeadsArray& _p_curr, double* _p_next, double _damping, vertex* _V) :
        p_curr(_p_curr), p_next(_p_next), damping(_damping), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next[d] += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        /* writeAdd(&p_next[d], (damping*p_curr[s]/V[s].getOutDegree()));
		
        does not work
		volatile unsigned int oh, ot, nh, nt;
        volatile double next;
        double c = p_curr[s]/V[s].getOutDegree();
        do
        {
            next = (p_next[d] + c);
            oh = p_next[d].head;
            ot = p_next[d].tail;
            volatile std::uint64_t nl = *reinterpret_cast<volatile std::uint64_t*>(&next);
            nh = (nl >> 32U);
            nt = (nl & 0xFFFFFFFF);
        } while((!CAS(&p_next[d].head, oh, nh) || !CAS(&p_next[d].tail, ot, nt)) && (p_next[d] != next)); */

        cerr << "updateAtomic not implemented for ManSeg type.. yet\n";
        // abort();
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += damping*(p_curr[s]/V[s].getOutDegree());
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        // Cache is used only in sequential mode
        p_next[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//resets p
template<class ArrayType>
struct PR_Vertex_Reset
{
    ArrayType& p_curr;
    PR_Vertex_Reset(ArrayType& _p_curr) :
        p_curr(_p_curr) {}
    inline bool operator () (intT i)
    {
        p_curr[i] = 0.0;
        return 1;
    }
};

// double version
struct PR_Vertex_Reset_d
{
    double* p_curr;
    PR_Vertex_Reset_d(double* _p_curr) :
        p_curr(_p_curr) {}
    inline bool operator () (intT i)
    {
        p_curr[i] = 0.0;
        return 1;
    }
};

template<class ArrayType>
double seqsum(ArrayType& a, intT s, intT e)
{
    double d = 0.;
    double err = 0.;
    for( intT i=s; i < e; ++i ) 
    {
	    //The code below achieves
	    // d += a[i];
        // but does so with high accuracy
        
        double tmp = d;
        double y = a[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}
template<class ArrayType>
double sumArray(const partitioner &part, ArrayType& a, intT n)
{
    double d = 0.;
	double err = 0.;
	double tmp, y;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqsum( a, s, e);
    });
    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err  += y;
    }
    delete [] psum;
    return d;
}

template<class ArrayType>
double seqnormdiff(ArrayType& a, ArrayType& b, intT s, intT e)
{
    double d = 0.;
    double err = 0.;
    for( intT i=s; i < e; ++i ) 
    {
        //The code does d += fabs(a[i]- b[i]);
        // but does so with high accuracy
        double tmp = d;
        double y = fabs(a[i] - b[i]) + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}
template<class ArrayType>
double normDiff(const partitioner &part, ArrayType& a, ArrayType& b, intT n)
{
    double d = 0.;
    int p= part.get_num_partitions();
    double *psum = new double [p];
    double err = 0.;
    double tmp, y;
    /*calculate sum per partitions*/
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqnormdiff( a, b, s, e);
    } );

    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}

double interim_seqnormdiff(PairsArray& a, double* b, intT s, intT e)
{
    double d = 0.;
    double err = 0.;
    for( intT i=s; i < e; ++i ) 
    {
        //The code does d += fabs(a[i]- b[i]);
        // but does so with high accuracy
        double tmp = d;
        double y = fabs(a[i] - b[i]) + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}
double interim_normDiff(const partitioner &part, PairsArray& a, double* b, intT n)
{
    double d = 0.;
    int p= part.get_num_partitions();
    double *psum = new double [p];
    double err = 0.;
    double tmp, y;
    /*calculate sum per partitions*/
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = interim_seqnormdiff( a, b, s, e);
    } );

    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}

// double version
double seqnormdiff(double* a, double* b, intT n)
{
    double d = 0.;
    double err = 0.;
    for( intT i=0; i < n; ++i ) 
    {
        //does d += fabs(a[i]- b[i]); but with high accuracy
        double tmp = d;
        double y = fabs(a[i] - b[i]) + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}
double normDiff(const partitioner &part, double* a, double* b, intT n)
{
    double d = 0.;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    double err = 0.;
    double tmp, y;
    /* calculate sum per partition */
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqnormdiff( &a[s], &b[s], e-s);
    } );

    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}
// double version
double seqsum(double* a, intT n)
{
    double d = 0.;
    double err = 0.;
    for( intT i=0; i < n; ++i ) 
    {
	    //does d += a[i]; but with high accuracy
        double tmp = d;
        double y = a[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    return d;
}
double sumArray(const partitioner &part, double* a, intT n)
{
    double d = 0.;
	double err = 0.;
	double tmp, y;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    // get sum of each partition
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = seqsum( &a[s], e-s);
    });
    for( int i=0; i < p; ++i ) 
    {
        // does d += psum[i]; but with high accuracy
        tmp = d;
        y = psum[i] + err;
        d = tmp + y;
        err = tmp - d;
        err += y;
    }
    delete [] psum;
    return d;
}

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphTyp
    const partitioner &part = GA.get_partitioner();
    graph<vertex> & WG = GA.get_partition();
    const int perNode = part.get_num_per_node_partitions();
    const intT partElem = part.get_num_elements();
    intT n = GA.n;
    intT m = GA.m;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    //Data Array
    //x and y to do special allocation
    //frontier also need special node allocation
    //blocksize equal to the szie of each partitioned
    double one_over_n = 1/(double)n;

    ManSegArray p_curr(partElem);
    ManSegArray p_next(partElem);
	// alloc full
    p_curr.full = new double[partElem];
    p_next.full = new double[partElem];

    double delta = 2.0;

    loop(j, part, perNode, p_curr.heads[j] = one_over_n);
    loop(j, part, perNode, p_next.heads[j] = 0.0);
	loop(j, part, perNode, p_next.full[j] = 0.0);
    cerr << setprecision(16);

    int count=0;
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    while(count<MaxIter) // heads only
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, HeadsArray>(p_curr.heads,p_next.heads,damping,WG.V),m/20);
       
        // find value to scale PR vals by to make vector add to 1
        double scaleAdditive = (1 - sumArray<HeadsArray>(part, p_next.heads, n))*one_over_n;
        {
			loop(j, part, perNode, p_next.heads[j] += scaleAdditive);
		}

        // delta = abs(p_curr - p_next)
        delta = normDiff<HeadsArray>(part, p_curr.heads, p_next.heads, n);

        // reset p_curr and swap vertices
        {
			loop(j, part, perNode, p_curr.heads[j] = 0);
		}
		// swap parts manually
		// not doing so causes segmentation fault
        swap(p_curr.heads,p_next.heads);
        swap(p_curr.pairs,p_next.pairs);
        swap(p_curr.full,p_next.full);
        // manage frontier stuff
        Frontier.del();
        Frontier = output;

        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray<HeadsArray>(part, p_curr.heads, n) << "\n";
        // ensure swap & reset happens *before* breaking from loop
        
		if(delta <= AdaptivePrecisionBound)
		{
			cerr << "switching precision at iter " << count << "\n";
			break;
		}
    }

    // Interim Iteration required for switching precisions.
    /*
        if delta is close to current working precision, the following
        steps are executed: 

        1) run the PageRank iteration while reading in current
           precision, and writing in the new, extended precision;
        2) normalize the vector in the new precision, so we can
           ensure that the norm of the PageRank-vector p^k stays
           at 1; and
        3) set the new precision as current working precision;
           The normalization is necessary because, while writing back,
           the floating point representation is “cut”, which leads to a
           rounding towards zero. This, in turn, leads to ||p^k|| < 1
    */
	if(count < MaxIter)
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_Interim_F<vertex>(p_curr.heads,p_next.full,damping,WG.V),m/20);
       
        // find value to scale PR vals by to make vector add to 1
        double scaleAdditive = (1 - sumArray(part, p_next.full, n))*one_over_n;
        {
			loop(j, part, perNode, p_next.full[j] += scaleAdditive);
		}
        // calculate delta value between current and new pageranks
        delta = interim_normDiff(part, p_curr.pairs, p_next.full, n);

        // reset p_curr and swap vertices
        {
			loop(j, part, perNode, p_curr.full[j] = 0);
		}
        swap(p_curr.heads,p_next.heads);
        swap(p_curr.pairs,p_next.pairs);
        swap(p_curr.full,p_next.full);
        // manage frontier stuff
        Frontier.del();
        Frontier = output;

        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray(part, p_curr.full, n) << "\n";
    }

    while(count<MaxIter) // full precision
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F_d<vertex>(p_curr.full,p_next.full,damping,WG.V),m/20);

        // find value to scale PR vals by to make vector add to 1
        double scaleAdditive = (1 - sumArray(part, p_next.full, n))*one_over_n;
        {
        	loop(j, part, perNode, p_next.full[j] += scaleAdditive);
		}

        // delta = abs(p_curr - p_next)
        delta = normDiff(part, p_curr.full, p_next.full, n);
		if(delta < epsilon)
        {
            cerr << count << ": delta = " << delta << "\n";
            cerr << "successfully converged in " << count << " iterations\n";
            break;
        }
        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray(part, p_curr.full, n) << "\n";

        //reset p_curr
        {
			loop(j, part, perNode, p_curr.full[j] = 0);
		}
        swap(p_curr.heads,p_next.heads);
        swap(p_curr.pairs,p_next.pairs);
        swap(p_curr.full,p_next.full);
        // manage frontier stuff
        Frontier.del();
        Frontier = output;
    }

    // clean up memory
    Frontier.del();
	p_curr.delSegments();
	p_curr.del();
	p_next.delSegments();
	p_next.del();
}

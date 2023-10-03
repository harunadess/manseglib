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
int MaxIter=100;
template <class vertex, typename F_in, typename F_out>
struct PR_F
{
    F_in* p_curr_in;
    F_out* p_next_out;
    double damping;
    vertex* V;
    static const bool use_cache = true;

    struct cache_t
    {
        F_out p_next;
    };
    PR_F(F_in* _p_curr_in, F_out* _p_next_out, double _damping, vertex* _V) :
        p_curr_in(_p_curr_in), p_next_out(_p_next_out), damping(_damping), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next_out[d] += (F_out)( damping*(p_curr_in[s]/V[s].getOutDegree()) );
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        writeAdd<F_out>(&p_next_out[d], (F_out)( damping*(p_curr_in[s]/V[s].getOutDegree()) ));
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next_out[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += (F_out)(F_out)( damping*(p_curr_in[s]/V[s].getOutDegree()) );
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        // Cache is used only in sequential mode
        p_next_out[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//resets p
template<typename F>
struct PR_Vertex_Reset
{
    F* p_curr;
    PR_Vertex_Reset(F* _p_curr) :
        p_curr(_p_curr) {}
    inline bool operator () (intT i)
    {
        p_curr[i] = 0.0;
        return 1;
    }
};

template<typename F>
double seqsum(F* a, intT n)
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

template<typename F>
double sumArray(const partitioner &part, F* a, intT n)
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
        psum[k] = seqsum<F>( &a[s], e-s);
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

template<typename F1, typename F2>
double seqnormdiff(F1* a, F2* b, intT n)
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

template<typename F1, typename F2>
double normDiff(const partitioner &part, F1* a, F2* b, intT n)
{
    double d = 0.;
    int p = part.get_num_partitions();
    double *psum = new double [p];
    double err = 0.;
    double tmp, y;

	// requires auto declaration of function because
	// map_partition macro doesn't like there being more than one type
	// in the function template parameters
	auto snd = seqnormdiff<F1, F2>;
    /* calculate sum per partition */
    map_partition( k, part, {
        intT s = part.start_of(k);
        intT e = part.start_of(k+1);
        psum[k] = snd( &a[s], &b[s], e-s);
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
    const double float_limit = 5e-5f;
    //Data Array
    //x and y to do special allocation
    //frontier also need special node allocation
    //blocksize equal to the szie of each partitioned
    double one_over_n = 1.0/(double)n;

    // float version
    mmap_ptr<float> p_curr_f;
    p_curr_f.part_allocate (part);
    mmap_ptr<float> p_next_f;
    p_next_f.part_allocate (part);

    // double version
    mmap_ptr<double> p_curr_d;
    p_curr_d.part_allocate (part);
    mmap_ptr<double> p_next_d;
    p_next_d.part_allocate (part);

    double delta = 2.0;
    loop(j, part, perNode, p_curr_f[j] = one_over_n);
    loop(j, part, perNode, p_next_f[j] = 0.0f);
	loop(j, part, perNode, p_next_d[j] = 0.0);
    cerr << setprecision(16);

    int count=0;
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    // floats
    while(count<MaxIter)
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, float, float>(p_curr_f,p_next_f,damping,WG.V),m/20);
        
        // find value to scale PR vals by to make vector add to 1
        float scaleAdditive = (1 - sumArray<float>(part, p_next_f, n))*one_over_n;
        {
			loop(j, part, perNode, p_next_f[j] += scaleAdditive);
		}

        // delta = abs(p_curr - p_next)
        delta = normDiff<float, float>(part, p_curr_f, p_next_f, n);
    
        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray<float>(part, p_next_f, n) << "\n";

        // reset p_curr and swap vertices
        {
			loop(j, part, perNode, p_curr_f[j] = 0);
		}
        swap(p_curr_f, p_next_f);
        // manage fronteir stuff
        Frontier.del();
        Frontier = output;

        if(delta <= float_limit)
        {
            cerr << "hit float limit at iter " << count << "\n";
            break;
        }
    }

    // interim (float -> double) step
	if(count < MaxIter)
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, float, double>(p_curr_f,p_next_d,damping,WG.V),m/20);
        
        // find value to scale PR vals by to make vector add to 1
        double scaleAdditive = (1 - sumArray<double>(part, p_next_d, n))*one_over_n;
        {
			loop(j, part, perNode, p_next_d[j] += scaleAdditive);
		}

        // delta = abs(p_curr - p_next)
        delta = normDiff<float, double>(part, p_curr_f, p_next_d, n);

        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray<double>(part, p_next_d, n) << "\n";

        // reset p_curr and swap vertices
		{
			loop(j, part, perNode, p_curr_d[j] = 0);
		}
        swap(p_curr_d, p_next_d);
        // manage fronteir stuff
        Frontier.del();
        Frontier = output;
    }

    // doubles
    while(count<MaxIter)
    {
        ++count;

        // p_next[d] += damping * (p_curr[s]/V[s].getOutDegree())
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, double, double>(p_curr_d,p_next_d,damping,WG.V),m/20);
        
        // find value to scale PR vals by to make vector add to 1
        double scaleAdditive = (1 - sumArray<double>(part, p_next_d, n))*one_over_n;
        {
			loop(j, part, perNode, p_next_d[j] += scaleAdditive);
		}

        // delta = abs(p_curr - p_next)
        delta = normDiff<double, double>(part, p_curr_d, p_next_d, n);
		if(delta < epsilon)
        {
            cerr << count << ": delta = " << delta << "\n";
            cerr << "successfully converged in " << count << " iterations\n";
            break;
        }
        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray<double>(part, p_next_d, n) << "\n";

        // reset p_curr and swap vertices
        {
			loop(j, part, perNode, p_curr_d[j] = 0);
		}
        swap(p_curr_d, p_next_d);
        // manage fronteir stuff
        Frontier.del();
        Frontier = output;
    }


    Frontier.del();
    p_curr_f.del();
    p_next_f.del();
    p_curr_d.del();
    p_next_d.del();
}

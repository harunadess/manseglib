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
#include "mantissaSegmentation.hpp"
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
        p_next[d] += (damping*p_curr[s]/V[s].getOutDegree());
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        // writeAdd(&p_next[d], (damping*p_curr[s]/V[s].getOutDegree()));
		
        // does not work
		// volatile unsigned int oh, ot, nh, nt;
        // volatile double next;
        // double c = p_curr[s]/V[s].getOutDegree();
        // do
        // {
        //     next = (p_next[d] + c);
        //     oh = p_next[d].head;
        //     ot = p_next[d].tail;
        //     volatile std::uint64_t nl = *reinterpret_cast<volatile std::uint64_t*>(&next);
        //     nh = (nl >> 32U);
        //     nt = (nl & 0xFFFFFFFF);
        // } while((!CAS(&p_next[d].head, oh, nh) || !CAS(&p_next[d].tail, ot, nt)) && (p_next[d] != next));

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
        cache.p_next += (damping*p_curr[s]/V[s].getOutDegree());
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

template<class ArrayType>
double sumArray(ArrayType& a, intT n)
{
    double d = 0.0;
    double err = 0.0;
    for(intT i = 0; i < n; ++i)
    {
        // does d += a[i] with high accuracy
        double temp = d;
        double y = a[i] + err;
        d = temp + y;
        err = temp - d;
        err += y;
    }
    return d;
}

template<class ArrayType>
double normDiff(ArrayType& a, ArrayType& b, intT n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
    {
        // does d += fabs(b[i] - a[i]) with high accuracy
        double temp = d;
        double y = std::fabs(b[i] - a[i]) + err;
        d = temp + y;
        err = temp - d;
        err += y;
    }
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
    double one_over_n = 1.0/(double)n;

    // define nicer types
    using HeadsArray = ManSeg::TwoSegArray<false>;
    using PairsArray = ManSeg::TwoSegArray<true>;
    // use combined conveniencetype
    ManSeg::ManSegArray p_curr(partElem);
    ManSeg::ManSegArray p_next(partElem);

    double delta = 2.0;
    int count=0;

    // for(intT i = 0; i < partElem; ++i)
    // {
    //     p_curr.heads[i] = one_over_n;
    //     p_next.heads[i] = 0.0;
    // }
    loop(j, part, perNode, p_curr.heads[j] = one_over_n);
    loop(j, part, perNode, p_next.heads[j] = 0.0);

    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);

    cerr << setprecision(16);
    while(count<MaxIter) // heads only
    {
        // power method step (main page rank step)
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, HeadsArray>(p_curr.heads,p_next.heads,damping,WG.V),m/20);
       
        // calculate current sum of new pageranks (will sum to < 1)
        double nextPrSum = sumArray(p_next.heads, partElem);
        // find scaling value to ensure new pageranks sum to 1
        double w = (1.0 - nextPrSum)*one_over_n;
       
        // for(intT i = 0; i < partElem; ++i)
        //     p_next.heads[i] += w;
        loop(j, part, perNode, p_next.heads[j] += w);
        
        // calculate delta value between current and new pageranks
        delta = normDiff(p_curr.heads, p_next.heads, partElem);
        ++count;

        // reset p_curr and swap vertices
        vertexMap(part, Frontier,PR_Vertex_Reset<HeadsArray>(p_curr.heads));
        swap(p_curr.heads,p_next.heads);
        // manage frontier stuff
        Frontier.del();
        Frontier = output;

        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray(p_curr.heads, partElem) << "\n";
        // ensure swap & reset happens *before* breaking from loop
        if(delta < ManSeg::MaxSingleSegmentPrecision)
        {   
            cerr << "hit limit of single segment\n\n";
            break;
        }
    }

    while(count<MaxIter) // full precision
    {
        // power method step (main page rank step)
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex, PairsArray>(p_curr.pairs,p_next.pairs,damping,WG.V),m/20);

        // calculate current sum of new pageranks (will sum to < 1)
        double nextPrSum = sumArray(p_next.pairs, partElem);
        // find scaling value to ensure new pageranks sum to 1
        double w = (1.0 - nextPrSum)*one_over_n;

        // scale values  < -- todo: make parallel via loop macro or something
        // for(intT i = 0; i < partElem; ++i)
        //     p_next.pairs[i] += w;
        loop(j, part, perNode, p_next.pairs[j] += w);

        // calculate delta value between current and new pageranks
        delta = normDiff(p_curr.pairs, p_next.pairs, partElem);
        ++count;

        //reset p_curr
        vertexMap(part, Frontier,PR_Vertex_Reset<PairsArray>(p_curr.pairs));
        swap(p_curr.pairs,p_next.pairs);
        // manage frontier stuff
        Frontier.del();
        Frontier = output;

        cerr << count << ": delta = " << delta << "  xnorm = " << sumArray(p_curr.pairs, partElem) << "\n";
        if(delta < epsilon)
        {
            cerr << "successfully converged\n";
            break;
        }
    }

    // clean up memory
    Frontier.del();
    p_curr.del();
    p_next.del();
}

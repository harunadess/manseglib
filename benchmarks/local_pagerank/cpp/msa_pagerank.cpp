#include <iostream>
#include <iomanip>
#include <chrono>

#include <fstream>
#include <string>
#include <regex>
#include <cmath>

#include "../../../manseglib.hpp"
#include "quicksort.h"

using namespace std;
using namespace ManSeg;

constexpr double d = 0.85;
constexpr double tol = 1e-7;
// constexpr double tol = 1e-10;
constexpr int maxIter = 100;
// constexpr int maxIter = 200;

enum MatrixType {COO, CSR, CSC, SNAP_COO, SNAP_CSR, SNAP_CSC};

template<int matrixType>
class SparseMatrix;

template<>
class SparseMatrix<COO>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in COO" << endl;
            exit(1);
        }
        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        source = new int[numEdges];
        destination = new int[numEdges];

        for(int i = 0; i < numEdges; ++i)
            f >> source[i] >> destination[i];

        f.close();
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            newPr[dest] += d*(prevPr[src]/outdeg[src]);
        }
    }

    void iterate(double d, double* _prevPr, double* _newPr, int outdeg[], double* contr)
    {    
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            _newPr[dest] += d*(_prevPr[src]/outdeg[src]);
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            double val = d*(prevPr[src]/outdeg[src]);
            newPr[dest] += val;
        }
    }

    int numVertices;
    int numEdges;
private:
    int* source;
    int* destination;
};

template<>
class SparseMatrix<CSR>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSR" << endl;
            exit(1);
        }

        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        index = new int[numVertices + 1];
        dest = new int[numEdges];

        index[0] = 0;
        int pos = 0;

        getline(f, line);
        for(int i = 0; i < numVertices; ++i)
        {
            getline(f, line);
            istringstream iss(line);
            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

            pos += (results.size() - 1);
            index[i + 1] = pos;

            for(int j = 1; j < results.size(); ++j)
                dest[index[i] + (j - 1)] = stoi(results[j]);
        }

        f.close();
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numVertices; ++i)
            outdeg[i] = (index[i + 1] - index[i]);
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    void iterate(double d, double* prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    int numVertices;
    int numEdges;
private:
    int* index;
    int* dest;
};

template<>
class SparseMatrix<CSC>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSC" << endl;
            exit(1);
        }

        string line;
        getline(f, line);

        f >> numVertices;
        f >> numEdges;

        index = new int[numVertices + 1];
        source = new int[numEdges];

        index[0] = 0;
        int pos = 0;

        getline(f, line);
        for(int i = 0; i < numVertices; ++i)
        {
            getline(f, line);
            istringstream iss(line);
            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

            pos += (results.size() - 1);
            index[i + 1] = pos;

            for(int j = 1; j < results.size(); ++j)
                source[index[i] + (j - 1)] = stoi(results[j]);
        }

        f.close();
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
        contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    void iterate(double d, double* prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    int numVertices;
    int numEdges;
private:
    int* index;
    int* source;
};


/*
    **********************************************************
                    SNAP input file format
    **********************************************************
*/
template<>
class SparseMatrix<SNAP_COO>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in COO" << endl;
            exit(1);
        }
        
        string skipLine, nodes, edges;
        getline(f, skipLine); // skip first line

        f >> nodes >> edges;
        getline(f, skipLine); // skip next descriptor

        string n = regex_replace(nodes, regex("[^0-9]*([0-9]+).*"), "$1");
        numVertices = atoi(n.c_str());
        
        n = regex_replace(edges, regex("[^0-9]*([0-9]+).*"), "$1");
        numEdges = atoi(n.c_str());

        cout << "numVertices=" << numVertices << "\n";
        cout << "numEdges=" << numEdges << "\n";

        source = new int[numEdges];
        destination = new int[numEdges];

        getline(f, skipLine); // skip next descriptor

        for(int i = 0; i < numEdges; ++i)
            f >> source[i] >> destination[i];

        f.close();
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            newPr[dest] += d*(prevPr[src]/outdeg[src]);
        }
    }

    void iterate(double d, double* _prevPr, double* _newPr, int outdeg[], double* contr)
    {    
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            _newPr[dest] += d*(_prevPr[src]/outdeg[src]);
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            double val = d*(prevPr[src]/outdeg[src]);
            newPr[dest] += val;
        }
    }

    int numVertices;
    int numEdges;
private:
    int* source;
    int* destination;
};

template<>
class SparseMatrix<SNAP_CSR>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSR" << endl;
            exit(1);
        }

       string line, nodes, edges;
        getline(f, line); // skip first line

        f >> nodes >> edges;
        getline(f, line); 

        string n = regex_replace(nodes, regex("[^0-9]*([0-9]+).*"), "$1");
        numVertices = atoi(n.c_str());
        
        n = regex_replace(edges, regex("[^0-9]*([0-9]+).*"), "$1");
        numEdges = atoi(n.c_str());

        cout << "numVertices=" << numVertices << "\n";
        cout << "numEdges=" << numEdges << "\n";

        index = new int[numVertices + 1];
        dest = new int[numEdges];

        index[0] = 0; // first index

        int indexPos = 0;
        int edgeCounter = 0;
        int s; // hold source for comparison

        getline(f, line); // skip next descriptor
        for(int i = 0; i < numEdges; ++i)
        {
            f >> s >> dest[i];

            while(indexPos < s)
                index[++indexPos] = (edgeCounter);

            ++edgeCounter;
        }

        while(indexPos < (numVertices + 1))
            index[++indexPos] = numEdges;

        f.close();
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numVertices; ++i)
            outdeg[i] = (index[i + 1] - index[i]);
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    void iterate(double d, double* prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[dest[j]] += contr[i];
        }
    }

    int numVertices;
    int numEdges;
private:
    int* index;
    int* dest;
};

template<>
class SparseMatrix<SNAP_CSC>
{
public:
    SparseMatrix(string file)
    {
        ifstream f;
        f.open(file, std::ios_base::in);
        if(!f.is_open())
        {
            cerr << "error opening file in CSR" << endl;
            exit(1);
        }

       string line, nodes, edges;
        getline(f, line); // skip first line

        f >> nodes >> edges;
        getline(f, line);

        string n = regex_replace(nodes, regex("[^0-9]*([0-9]+).*"), "$1");
        numVertices = atoi(n.c_str());
        
        n = regex_replace(edges, regex("[^0-9]*([0-9]+).*"), "$1");
        numEdges = atoi(n.c_str());

        cout << "numVertices=" << numVertices << "\n";
        cout << "numEdges=" << numEdges << "\n";

        index = new int[numVertices + 1];
        // source = new int[numEdges];

        getline(f, line);  // skip next line

        // read in edges as COO
        int* src = new int[numEdges];
        int* dst = new int[numEdges];

        for(int i = 0; i < numEdges; ++i)
            f >> src[i] >> dst[i];

        // perform sort on dst, so we can create a CSC representation
        quicksort_pair(src, dst, 0, (numEdges-1));

        index[0] = 0; // first index
        int indexPos = 0;
        int edgeCounter = 0;
        for(int i = 0; i < numEdges; ++i)
        {
            while(indexPos < dst[i])
                index[++indexPos] = (edgeCounter);
            ++edgeCounter;
        }

        while(indexPos < (numVertices + 1))
            index[++indexPos] = numEdges;

        source = src;

        // not super pretty but we need this for the same ordering as
        // you would get with CSR/COO, otherwise deltas/pagerank values
        // end up being off due to error
        for(int i = 0; i < (numVertices); ++i)
            quicksort(source, index[i], index[i+1]-1);

        f.close();

        delete[] dst;
    }

    void calculateOutDegree(int outdeg[])
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    template<class TwoSegArray>
    void iterate(double d, TwoSegArray& prevPr, TwoSegArray& newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    void iterate(double d, double* prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    void iterate_fullSet(double d, TwoSegArray<false>& prevPr, double* newPr, int outdeg[], double* contr)
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        int curr, next;
        for(int i = 0; i < numVertices; ++i)
        {
            curr = index[i];
            next = index[i + 1];
            for(int j = curr; j < next; ++j)
                newPr[i] += contr[source[j]];
        }
    }

    int numVertices;
    int numEdges;

private:
    int* index;
    int* source;
};

template<class TwoSegArray>
double sum(TwoSegArray& a, int& n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
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
double sum(double* a, int& n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
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

template<class TwoSegArray>
// x, y, n
double normDiff(TwoSegArray& a, TwoSegArray& b, int& n)
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
template<class TwoSegArray>
double normDiff(TwoSegArray& a, double* b, int& n)
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
// x, y, n
double normDiff(double* a, double* b, int& n)
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

template<class SparseMatrix>
void pr(SparseMatrix* matrix, std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds>& tmStart, string& inputFile)
{
    auto totalSt = tmStart;
    auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Reading input: " << tmInput << " seconds" << endl;
    tmStart = chrono::high_resolution_clock::now();

    int n = matrix->numVertices;
    ManSegArray x(n); // pagerank
	x.full = new double[n];
    // ManSegArray v(n); 
    ManSegArray y(n);
	y.full = new double[n];

    int* outdeg = new int[n];
    double* contr = new double[n];

	double prevDelta = 2.0, delta = 2.0;
	// int checkFrequency = 4;
	// double threshold = 3; 

	// int switch_point = 11;
	double thresh = 0.25;

    int iter = 0;

    double oneOverN = (1.0/static_cast<double>(n));
    for(int i = 0; i < n; ++i)
    {
        x.heads[i] = oneOverN;
        outdeg[i] = y.heads[i] = 0.0;
    }

    matrix->calculateOutDegree(outdeg);

    auto tmInit = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Initialisation: " << tmInit << " seconds" << endl;
    auto tmTotal = 0.0;
    tmStart = chrono::high_resolution_clock::now();

    auto iterateS = chrono::high_resolution_clock::now();
    
	int step = 0;
    while(iter < maxIter) // use heads only
    {
        matrix->iterate(d, x.heads, y.heads, outdeg, contr);

		double w = (1.0 - sum(y.heads, n))*oneOverN;
        for(int i = 0; i < n; ++i)
            y.heads[i] += w;

        delta = normDiff(x.heads, y.heads, n);
        ++iter;

        for(int i = 0; i < n; ++i)
        {
            x.heads[i] = y.heads[i];
            y.heads[i] = 0.0;
        }

        auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x.heads, n) 
            << " time=" << tmStep << " seconds" << endl;
		// cout << "iteration " << iter << ": delta=" << delta << "\n";
        
        tmStart = chrono::high_resolution_clock::now();

		++step;
		
		// if(step == checkFrequency)
		{
			// double ratio = prevDelta/delta;
			double change = (prevDelta-delta)/prevDelta;
			// cout  << "ratio = " << ratio << "\n";
			cout << "% change = " << change << "\n";
			if(change <= thresh)
				break;
			// step = 0;

			// if(delta < thresh)
				// break;
		}
		prevDelta = delta;
    }
    
    cout << "\n=========================\nInterim Step\n=========================" << endl;
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
    {
        matrix->iterate_fullSet(d, x.heads, y.full, outdeg, contr);
    
        double w = (1.0 - sum(y.full, n))*oneOverN;
        for(int i = 0; i < n; ++i)
            y.full[i] += w;

        delta = normDiff(x.pairs, y.full, n);
        ++iter;

        for(int i = 0; i < n; ++i)
        {
            x.full[i] = y.full[i];
            y.full[i] = 0.0;
        }

        auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x.full, n) 
            << " time=" << tmStep << " seconds" << endl;
		// cout << "iteration " << iter << ": delta=" << delta << "\n";

        tmStart = chrono::high_resolution_clock::now();

		// {
		// 	double ratio = prevDelta/delta;
		// 	cout  << "ratio = " << ratio << "\n";
		// 	cout << "% change = " << (prevDelta-delta)/prevDelta << "\n";
		// 	// if(ratio <= threshold)
		// 	// 	break;
		// 	// step = 0;
		// }
		// prevDelta = delta;
    }
    
    cout << "\n=========================\nIncreased Precision\n=========================" << endl;
    while(iter < maxIter && delta > tol) // use full precision
    {
        // matrix->iterate(d, &x.pairs, &y.pairs, outdeg, contr);
        matrix->iterate(d, x.full, y.full, outdeg, contr);

        double w = (1.0 - sum(y.full, n))*oneOverN;
        for(int i = 0; i < n; ++i)
            y.full[i] += w;

        delta = normDiff(x.full, y.full, n);
        ++iter;

        for(int i = 0; i < n; ++i)
        {
            x.full[i] = y.full[i];
            y.full[i] = 0.0;
        }
 
        auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x.full, n) 
            << " time=" << tmStep << " seconds" << endl;
		// cout << "iteration " << iter << ": delta=" << delta << "\n";

        tmStart = chrono::high_resolution_clock::now();

		// {
		// 	double ratio = prevDelta/delta;
		// 	cout  << "ratio = " << ratio << "\n";
		// 	cout << "% change = " << (prevDelta-delta)/prevDelta << "\n";
		// 	// if(ratio <= threshold)
		// 	// 	break;


		// 	// step = 0;
		// }
		// prevDelta = delta;
    }

    auto endT = chrono::high_resolution_clock::now();

    auto iterateT = chrono::duration_cast<chrono::nanoseconds>(endT - iterateS).count()*1e-9;
    auto totalT = chrono::duration_cast<chrono::nanoseconds>(endT - totalSt).count()*1e-9;

    cout << "\nTotal time:" << totalT << " seconds\n"
        << "Total iterate time:" << iterateT << " seconds\n"
        << "Total seq time:" << (totalT - iterateT) << " seconds" << endl;

    if(delta > tol)
        cerr << "error: solution has not converged" << endl;

    // write to file
    // string outPath = "";
    // for(int i = inputFile.length()-1; i >= 0; --i)
    // {
    //     char c = inputFile.at(i);
    //     if(c == '/' || c == '\\') break;
    //     outPath += c;
    // }
    // reverse(outPath.begin(), outPath.end());

    // ofstream of;
    // of.open(("./results/msa_" + outPath + ".prvals"));

    // of << setprecision(16);
    // for(int i = 0; i < n; ++i)
    //     of << i << " " << x.full[i] << "\n";

    // of.close();

    // x.del();
    // y.del();
    delete[] contr;
	x.delSegments();
	y.delSegments();
}

int main(int argc, char** argv)
{
    cout << setprecision(16);
    cerr << setprecision(16);

    if(argc < 4)
    {
        cerr << "usage: ./msa_pagerank <type> <format> <input_file>" << endl;
        return 1;
    }

    string type = argv[1];
    string format = argv[2];
    string inputFile = argv[3];

    cout << "\nType: " << type
    << "\nFormat: " << format 
    << "\nInput file: " << inputFile
    << endl;

    auto tmStart = chrono::high_resolution_clock::now();

    SparseMatrix<COO>* coo;
    SparseMatrix<CSR>* csr;
    SparseMatrix<CSC>* csc;
    SparseMatrix<SNAP_COO>* snap_coo;
    SparseMatrix<SNAP_CSR>* snap_csr;
    SparseMatrix<SNAP_CSC>* snap_csc;
    if(type.compare("SNAP") == 0) // type is SNAP
    {
        cout << "SNAP format" << endl;
        if(format.compare("CSR") == 0)
        {
            snap_csr = new SparseMatrix<SNAP_CSR>(inputFile);
            pr(snap_csr, tmStart, inputFile);
        }
        else if(format.compare("CSC") == 0)
        {
            snap_csc = new SparseMatrix<SNAP_CSC>(inputFile);
            pr(snap_csc, tmStart, inputFile);
        }
        else if(format.compare("COO") == 0)
        {
            snap_coo = new SparseMatrix<SNAP_COO>(inputFile);
            pr(snap_coo, tmStart, inputFile);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }
    else
    {
        cout << "GraphType format" << endl;
        if(format.compare("CSR") == 0)
        {
            csr = new SparseMatrix<CSR>(inputFile);
            pr(csr, tmStart, inputFile);
        }
        else if(format.compare("CSC") == 0)
        {
            csc = new SparseMatrix<CSC>(inputFile);
            pr(csc, tmStart, inputFile);
        }
        else if(format.compare("COO") == 0)
        {
            coo = new SparseMatrix<COO>(inputFile);
            pr(coo, tmStart, inputFile);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }
    

    return 0;
}

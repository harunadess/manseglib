#include <iostream>
#include <iomanip>
#include <chrono>

#include <fstream>
#include <string>
#include <regex>

#include "../../../mantissaSegmentation.hpp"

using namespace std;
using namespace ManSeg;

constexpr double d = 0.85;
constexpr double tol = 1e-7;
// constexpr double tol = 1e-10;
constexpr int maxIter = 100;
// constexpr int maxIter = 200;

// double* npr, *opr, *cnt;

// SNAP csc use:
// convert into coo
// then sort on dest
// then convert to csc
// initialisation will be slow, but should be O(n log n)

// todo: convert this from virtual nonsense to standalone objects/classes
// there isn't really any need for there to be inhertiance here,
// considering SparseMatrix is an abstract.
// just convert func pr() to be template
// and instead of SparseMatrix* use T*
// provided T has the functions it can be used..
class SparseMatrix
{
public:
    virtual void calculateOutDegree(int outdeg[]) = 0;
    virtual void iterate(double d, TwoSegArray<false>& prevPr, TwoSegArray<false>& newPr, int outdeg[], TwoSegArray<false>& contr) = 0;
    virtual void iterate(double d, TwoSegArray<true>& prevPr, TwoSegArray<true>& newPr, int outdeg[], TwoSegArray<true>& contr) = 0;

    int numVertices;
    int numEdges;
};

class SparseMatrixCOO : public SparseMatrix
{
public:
    SparseMatrixCOO(string file)
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
    }

    virtual void calculateOutDegree(int outdeg[]) override
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    virtual void iterate(double d, TwoSegArray<false>& prevPr, TwoSegArray<false>& newPr, int outdeg[], TwoSegArray<false>& contr) override
    {
        for(int i = 0; i < numEdges; ++i)
            newPr[destination[i]] += d*(prevPr[source[i]]/outdeg[source[i]]);

        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numEdges; ++i)
        //     npr[destination[i]] += d*(opr[source[i]]/outdeg[source[i]]);

        // newPr = npr;
    }

    virtual void iterate(double d, TwoSegArray<true>& prevPr, TwoSegArray<true>& newPr, int outdeg[], TwoSegArray<true>& contr) override
    {
        for(int i = 0; i < numEdges; ++i)
            newPr[destination[i]] += d*(prevPr[source[i]]/outdeg[source[i]]);

        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numEdges; ++i)
        //     npr[destination[i]] += d*(opr[source[i]]/outdeg[source[i]]);

        // newPr = npr;
    }


private:
    int* source;
    int* destination;
};

class SparseMatrixCSR : public SparseMatrix
{
public:
    SparseMatrixCSR(string file)
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
    }

    virtual void calculateOutDegree(int outdeg[]) override
    {
        for(int i = 0; i < numVertices; ++i)
            outdeg[i] = (index[i + 1] - index[i]);
    }

    virtual void iterate(double d, TwoSegArray<false>& prevPr, TwoSegArray<false>& newPr, int outdeg[], TwoSegArray<false>& contr) override
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        for(int i = 0; i < numVertices; ++i)
        {
            for(int j = index[i]; j < index[i + 1]; ++j)
                newPr[dest[j]] += contr[i];
        }

        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numVertices; ++i)
        //     cnt[i] = d*(opr[i]/outdeg[i]);

        // for(int i = 0; i < numVertices; ++i)
        //     for(int j = index[i]; j < index[i + 1]; ++j)
        //         npr[dest[j]] += cnt[i];

        // newPr = npr;
    }

    virtual void iterate(double d, TwoSegArray<true>& prevPr, TwoSegArray<true>& newPr, int outdeg[], TwoSegArray<true>& contr) override
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        for(int i = 0; i < numVertices; ++i)
        {
            for(int j = index[i]; j < index[i + 1]; ++j)
                newPr[dest[j]] += contr[i];
        }

        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numVertices; ++i)
        //     cnt[i] = d*(opr[i]/outdeg[i]);

        // for(int i = 0; i < numVertices; ++i)
        //     for(int j = index[i]; j < index[i + 1]; ++j)
        //         npr[dest[j]] += cnt[i];

        // newPr = npr;
    }

private:
    int* index;
    int* dest;
};

class SparseMatrixCSC : public SparseMatrix
{
public: 
    SparseMatrixCSC(string file)
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
    }

    virtual void calculateOutDegree(int outdeg[]) override
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    virtual void iterate(double d, TwoSegArray<false>& prevPr, TwoSegArray<false>& newPr, int outdeg[], TwoSegArray<false>& contr) override
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        for(int i = 0; i < numVertices; ++i)
        {
            for(int j = index[i]; j < index[i + 1]; ++j)
                newPr[i] += contr[source[j]];
        }
        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numVertices; ++i)
        //     cnt[i] = d*(opr[i]/outdeg[i]);

        // for(int i = 0; i < numVertices; ++i)
        //     for(int j = index[i]; j < index[i + 1]; ++j)
        //         npr[i] += cnt[source[j]];

        // newPr = npr;
    }

    virtual void iterate(double d, TwoSegArray<true>& prevPr, TwoSegArray<true>& newPr, int outdeg[], TwoSegArray<true>& contr) override
    {
        for(int i = 0; i < numVertices; ++i)
            contr[i] = d*(prevPr[i]/outdeg[i]);

        for(int i = 0; i < numVertices; ++i)
        {
            for(int j = index[i]; j < index[i + 1]; ++j)
                newPr[i] += contr[source[j]];
        }
        // prevPr.toNewDoubleArray(opr);
        // newPr.toNewDoubleArray(npr);

        // for(int i = 0; i < numVertices; ++i)
        //     cnt[i] = d*(opr[i]/outdeg[i]);

        // for(int i = 0; i < numVertices; ++i)
        //     for(int j = index[i]; j < index[i + 1]; ++j)
        //         npr[i] += cnt[source[j]];

        // newPr = npr;
    }

private:
    int* index;
    int* source;
};


// SNAP reading version of things
class SNAPSparseMatrixCOO : public SparseMatrix
{
public:
    SNAPSparseMatrixCOO(string file)
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
        getline(f, skipLine);

        string n = regex_replace(nodes, regex("[^0-9]*([0-9]+).*"), "$1");
        numVertices = atoi(n.c_str());
        
        n = regex_replace(edges, regex("[^0-9]*([0-9]+).*"), "$1");
        numEdges = atoi(n.c_str());

        cout << "numVertices=" << numVertices << "\n";
        cout << "numEdges=" << numEdges << "\n";

        source = new int[numEdges];
        destination = new int[numEdges];

        for(int i = 0; i < numEdges; ++i)
            f >> source[i] >> destination[i];
    }

    virtual void calculateOutDegree(int outdeg[]) override
    {
        for(int i = 0; i < numEdges; ++i)
            ++outdeg[source[i]];
    }

    virtual void iterate(double d, TwoSegArray<false>& prevPr, TwoSegArray<false>& newPr, int outdeg[], TwoSegArray<false>& contr) override
    {
        for(int i = 0; i < numEdges; ++i)
            newPr[destination[i]] += d*(prevPr[source[i]]/outdeg[destination[i]]);
    }

    virtual void iterate(double d, TwoSegArray<true>& prevPr, TwoSegArray<true>& newPr, int outdeg[], TwoSegArray<true>& contr) override
    {
        for(int i = 0; i < numEdges; ++i)
            newPr[destination[i]] += d*(prevPr[source[i]]/outdeg[destination[i]]);
    }


private:
    int* source;
    int* destination;
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

template<class TwoSegArray>
double normDiff(TwoSegArray& a, TwoSegArray& b, int& n)
{
    double d = 0.0;
    double err = 0.0;
    for(int i = 0; i < n; ++i)
    {
        // does d += abs(b[i] - a[i]) with high accuracy
        double temp = d;
        double y = abs(b[i] - a[i]) + err;
        d = temp + y;
        err = temp - d;
        err += y;
    }
    return d;
}

void pr(SparseMatrix* matrix, std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds>& tmStart)
{
    auto totalSt = tmStart;
    auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Reading input: " << tmInput << " seconds" << endl;
    tmStart = chrono::high_resolution_clock::now();

    int n = matrix->numVertices;
    TwoSegArray<false>x_f(n); // pagerank
    TwoSegArray<true>x_t = x_f.createFullPrecision();     // <--- should *more or less* equivalent to just creating true then assigning things
    TwoSegArray<false>v_f(n);                //      i.e. this should work fine and if it doesn't there's something wrong
    TwoSegArray<true>v_t = v_f.createFullPrecision();
    TwoSegArray<false>y_f(n); // new pagerank
    TwoSegArray<true>y_t = y_f.createFullPrecision();
    int* outdeg = new int[n];
    TwoSegArray<false>contr_f(n); // contribution for each vertex
    TwoSegArray<true>contr_t = contr_f.createFullPrecision();

    // npr = new double[n];
    // opr = new double[n];
    // cnt = new double[n];

    double delta = 2.0;
    int iter = 0;

    for(int i = 0; i < n; ++i)
    {
        x_f[i] = v_t[i] = (1.0 / (double)(n));
        outdeg[i] = y_f[i] = 0.0;
    }

    matrix->calculateOutDegree(outdeg);

    auto tmInit = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Initialisation: " << tmInit << " seconds" << endl;
    auto tmTotal = 0.0;
    tmStart = chrono::high_resolution_clock::now();

    auto iterateS = chrono::high_resolution_clock::now();
    
    
    // set this to true to only use full precision
    bool changePrecision = false;
    
    
    while(iter < maxIter && delta > tol)
    {
        if(!changePrecision)
        {
            matrix->iterate(d, x_f, y_f, outdeg, contr_f);

            double w = (1.0 - sum(y_f, n));
            for(int i = 0; i < n; ++i)
                y_f[i] += w * v_f[i];

            delta = normDiff(x_f, y_f, n);
            ++iter;
            
            if(delta <= 1e-5) 
            {
                changePrecision = true; // i think this is the limit of precision heads only can do
                cout << "==========\nincreased precision\n==========\n";
            }

            for(int i = 0; i < n; ++i)
            {
                x_f[i] = y_f[i];
                y_f[i] = 0.0;
            }

            auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
            cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x_f, n) 
                << " time=" << tmStep << " seconds" << endl;
        }
        else
        {
            matrix->iterate(d, x_t, y_t, outdeg, contr_t);

            double w = (1.0 - sum(y_t, n));
            for(int i = 0; i < n; ++i)
                y_t[i] += w * v_t[i];

            delta = normDiff(x_t, y_t, n);
            ++iter;

            for(int i = 0; i < n; ++i)
            {
                x_t[i] = y_t[i];
                y_t[i] = 0.0;
            }

            auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
            cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x_t, n) 
                << " time=" << tmStep << " seconds" << endl;
        }

        tmStart = chrono::high_resolution_clock::now();
    }
    auto endT = chrono::high_resolution_clock::now();

    auto iterateT = chrono::duration_cast<chrono::nanoseconds>(endT - iterateS).count()*1e-9;
    auto totalT = chrono::duration_cast<chrono::nanoseconds>(endT - totalSt).count()*1e-9;

    cout << "\nTotal time:" << totalT << " seconds\n"
        << "Total iterate time:" << iterateT << " seconds\n"
        << "Total seq time:" << (totalT - iterateT) << " seconds" << endl;

    if(delta > tol)
        cerr << "error: solution has not converged" << endl;

    // for(int i = 0; i < n; ++i)
        // cerr << i << " " << x_t[i] << endl;

    x_t.del();
    y_t.del();
    v_t.del();
    contr_t.del();
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

    SparseMatrix* matrix;
    if(type.compare("SNAP") == 0)
    {
        if(format.compare("CSR") == 0)
        {
            // matrix = new SNAPSparseMatrixCSR(inputFile);
            // pr(matrix, tmStart);
        }
        else if(format.compare("CSC") == 0)
        {
            // matrix = new SNAPSparseMatrixCSC(inputFile);
            // pr(matrix, tmStart);
        }
        else if(format.compare("COO") == 0)
        {
            matrix = new SNAPSparseMatrixCOO(inputFile);
            pr(matrix, tmStart);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }
    else
    {
        cout << "using default method...\n";
        if(format.compare("CSR") == 0)
        {
            matrix = new SparseMatrixCSR(inputFile);
            pr(matrix, tmStart);
        }
        else if(format.compare("CSC") == 0)
        {
            matrix = new SparseMatrixCSC(inputFile);
            pr(matrix, tmStart);
        }
        else if(format.compare("COO") == 0)
        {
            matrix = new SparseMatrixCOO(inputFile);
            pr(matrix, tmStart);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }
    
    

    return 0;
}

#include <iostream>
#include <iomanip>
#include <chrono>

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <regex>

#include "quicksort.h"

using namespace std;

constexpr double d = 0.85;     // this is already a real-ish world precision
constexpr double tol = 1e-7;
// constexpr double tol = 1e-10;  // a more real-ish world precision
constexpr int maxIter = 100;
// constexpr int maxIter = 200;   // double potential iterations to counter this..?

enum MatrixType {COO, CSR, CSC, SNAP_COO, SNAP_CSR, SNAP_CSC};

template<int matrixType>
class SparseMatrix
{
public:
    int numVertices;
    int numEdges;
};

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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            newPr[dest] += d*(prevPr[src]/outdeg[src]);
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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
    {
        int src, dest;
        for(int i = 0; i < numEdges; ++i)
        {
            src = source[i];
            dest = destination[i];
            newPr[dest] += d*(prevPr[src]/outdeg[src]);
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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
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

    void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[])
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

double normDiff(double* a, double* b, int& n)
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

void readPageRankValues(double* x, int n, string file)
{
    ifstream f;
    f.open(file, std::ios_base::in);
    if(!f.is_open())
    {
        cerr << "error opening pagerank values file" << endl;
        exit(1);
    }

    int v;
    for(int i = 0; i < n; ++i)
        f >> v >> x[i]; 

    f.close();
}

template<class SparseMatrix>
void pr(SparseMatrix* matrix, std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds>& tmStart, 
    string& inputFile, string& outputFile)
{
    auto totalSt = tmStart;
    auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Reading input: " << tmInput << " seconds" << endl;
    tmStart = chrono::high_resolution_clock::now();

    int n = matrix->numVertices;
    double* x = new double[n]; // pagerank
    double* v = new double[n];
    double* y = new double[n]; // new pagerank
    int* outdeg = new int[n];
    double* contr = new double[n]; // contribution for each vertex

    double delta = 2.0;
    int iter = 0;

    for(int i = 0; i < n; ++i)
    {
        x[i] = v[i] = 1.0 / (double)(n);
        y[i] = outdeg[i] = 0.0;
    }

    matrix->calculateOutDegree(outdeg);

    readPageRankValues(x, n, outputFile);

    auto tmInit = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
    cout << "Initialisation: " << tmInit << " seconds" << endl;
    auto tmTotal = 0.0;
    tmStart = chrono::high_resolution_clock::now();

    auto iterateS = chrono::high_resolution_clock::now();
    while(iter < maxIter && delta > tol)
    {
        matrix->iterate(d, x, y, outdeg, contr);

        double w = 1.0 - sum(y, n);
        for(int i = 0; i < n; ++i)
            y[i] += w * v[i];

        delta = normDiff(x, y, n);
        ++iter;

        for(int i = 0; i < n; ++i)
        {
            x[i] = y[i];
            y[i] = 0.0;
        }

        auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
        cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x, n) 
            << " time=" << tmStep << " seconds" << endl;

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

    // write to file
    string outPath = "";
    for(int i = inputFile.length()-1; i >= 0; --i)
    {
        char c = inputFile.at(i);
        if(c == '/' || c == '\\')
            break;
        outPath += c;
    }
    reverse(outPath.begin(), outPath.end());

    ofstream of;
    of.open(("./results/std_" + outPath + ".prvals"));

    of << setprecision(16);
    for(int i = 0; i < n; ++i)
        of << i << " " << x[i] << "\n";

    of.close();
}

int main(int argc, char** argv)
{
    cout << setprecision(16);
    cerr << setprecision(16);

    if(argc < 5)
    {
        cerr << "usage: ./pagerank <type> <format> <input_file> <output_file>" << endl;
        return 1;
    }

    string type = argv[1];
    string format = argv[2];
    string inputFile = argv[3];
    string outputFile = argv[4];

    cout << "Type: " << type
    << "\nFormat: " << format 
    << "\nInput file: " << inputFile
    << "\nOutput file: " << outputFile
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
        cout << "executing SNAP method.." << endl;
        if(format.compare("CSR") == 0)
        {
            snap_csr = new SparseMatrix<SNAP_CSR>(inputFile);
            pr(snap_csr, tmStart, inputFile, outputFile);
        }
        else if(format.compare("CSC") == 0)
        {
            snap_csc = new SparseMatrix<SNAP_CSC>(inputFile);
            pr(snap_csc, tmStart, inputFile, outputFile);
        }
        else if(format.compare("COO") == 0)
        {
            snap_coo = new SparseMatrix<SNAP_COO>(inputFile);
            pr(snap_coo, tmStart, inputFile, outputFile);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }
    else
    {
        cout << "executing default method.." << endl;
        if(format.compare("CSR") == 0)
        {
            csr = new SparseMatrix<CSR>(inputFile);
            pr(csr, tmStart, inputFile, outputFile);
        }
        else if(format.compare("CSC") == 0)
        {
            csc = new SparseMatrix<CSC>(inputFile);
            pr(csc, tmStart, inputFile, outputFile);
        }
        else if(format.compare("COO") == 0)
        {
            coo = new SparseMatrix<COO>(inputFile);
            pr(coo, tmStart, inputFile, outputFile);
        }
        else
        {
            cerr << "Unknown format: " << format << endl;
            return 1;
        }
    }

    return 0;
}

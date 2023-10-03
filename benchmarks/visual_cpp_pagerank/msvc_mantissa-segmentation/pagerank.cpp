//#include <iostream>
//#include <iomanip>
//#include <chrono>
//
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <vector>
//#include <iterator>
//
//using namespace std;
//
//constexpr double d = 0.85;
//constexpr double tol = 1e-7;
//constexpr int maxIter = 100;
//
//class SparseMatrix
//{
//public:
//    virtual void calculateOutDegree(int outdeg[]) = 0;
//    virtual void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[]) = 0;
//
//    int numVertices;
//    int numEdges;
//};
//
//class SparseMatrixCOO : public SparseMatrix
//{
//public:
//    SparseMatrixCOO(string file)
//    {
//        ifstream f;
//        f.open(file, std::ios_base::in);
//        if(!f.is_open())
//        {
//            cerr << "error opening file in COO" << endl;
//            exit(1);
//        }
//        string line;
//        getline(f, line);
//
//        f >> numVertices;
//        f >> numEdges;
//
//        source = new int[numEdges];
//        destination = new int[numEdges];
//
//        for(int i = 0; i < numEdges; ++i)
//            f >> source[i] >> destination[i];
//    }
//
//    virtual void calculateOutDegree(int outdeg[]) override
//    {
//        for(int i = 0; i < numEdges; ++i)
//            ++outdeg[source[i]];
//    }
//
//    virtual void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[]) override
//    {
//        int src, dest;
//        // pre calculate contribution
//        /*
//        for(int i = 0; i < numEdges; ++i)
//        {
//            s = source[i];
//            contr[s] = d*(x[s]/outdeg[s]);
//        }
//        */
//        
//        // pre calculating contribution doesn't help for coo
//        // does it?
//        for(int i = 0; i < numEdges; ++i)
//        {
//            src = source[i];
//            dest = destination[i];
//            newPr[dest] += d*(prevPr[src]/outdeg[src]);
//        }
//    }
//
//
//private:
//    int* source;
//    int* destination;
//};
//
//class SparseMatrixCSR : public SparseMatrix
//{
//public:
//    SparseMatrixCSR(string file)
//    {
//        ifstream f;
//        f.open(file, std::ios_base::in);
//        if(!f.is_open())
//        {
//            cerr << "error opening file in CSR" << endl;
//            exit(1);
//        }
//
//        string line;
//        getline(f, line);
//
//        f >> numVertices;
//        f >> numEdges;
//
//        index = new int[numVertices + 1];
//        dest = new int[numEdges];
//
//        index[0] = 0;
//        int pos = 0;
//
//        getline(f, line);
//        for(int i = 0; i < numVertices; ++i)
//        {
//            getline(f, line);
//            istringstream iss(line);
//            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
//
//            pos += (results.size() - 1);
//            index[i + 1] = pos;
//
//            for(size_t j = 1; j < results.size(); ++j)
//                dest[index[i] + (j - 1)] = stoi(results[j]);
//        }
//    }
//
//    virtual void calculateOutDegree(int outdeg[]) override
//    {
//        for(int i = 0; i < numVertices; ++i)
//            outdeg[i] = (index[i + 1] - index[i]);
//    }
//
//    virtual void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[]) override
//    {
//        for(int i = 0; i < numVertices; ++i)
//            contr[i] = d*(prevPr[i]/outdeg[i]);
//
//        int curr, next;
//        for(int i = 0; i < numVertices; ++i)
//        {
//            curr = index[i];
//            next = index[i + 1];
//            for(int j = curr; j < next; ++j)
//                newPr[dest[j]] += contr[i];
//        }
//    }
//
//private:
//    int* index;
//    int* dest;
//};
//
//class SparseMatrixCSC : public SparseMatrix
//{
//public:
//    SparseMatrixCSC(string file)
//    {
//        ifstream f;
//        f.open(file, std::ios_base::in);
//        if(!f.is_open())
//        {
//            cerr << "error opening file in CSC" << endl;
//            exit(1);
//        }
//
//        string line;
//        getline(f, line);
//
//        f >> numVertices;
//        f >> numEdges;
//
//        index = new int[numVertices + 1];
//        source = new int[numEdges];
//
//        index[0] = 0;
//        int pos = 0;
//
//        getline(f, line);
//        for(int i = 0; i < numVertices; ++i)
//        {
//            getline(f, line);
//            istringstream iss(line);
//            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
//
//            pos += (results.size() - 1);
//            index[i + 1] = pos;
//
//            for(size_t j = 1; j < results.size(); ++j)
//                source[index[i] + (j - 1)] = stoi(results[j]);
//        }
//    }
//
//    virtual void calculateOutDegree(int outdeg[]) override
//    {
//        for(int i = 0; i < numEdges; ++i)
//            ++outdeg[source[i]];
//    }
//
//    virtual void iterate(double d, double prevPr[], double newPr[], int outdeg[], double contr[]) override
//    {
//        for(int i = 0; i < numVertices; ++i)
//            contr[i] = d*(prevPr[i]/outdeg[i]);
//
//        int curr, next;
//        for(int i = 0; i < numVertices; ++i)
//        {
//            curr = index[i];
//            next = index[i + 1];
//            for(int j = curr; j < next; ++j)
//                newPr[i] += contr[source[j]];
//        }
//    }
//
//private:
//    int* index;
//    int* source;
//};
//
//
//double sum(double* a, int& n)
//{
//    double d = 0.0;
//    double err = 0.0;
//    for(int i = 0; i < n; ++i)
//    {
//        // does d += a[i] with high accuracy
//        double temp = d;
//        double y = a[i] + err;
//        d = temp + y;
//        err = temp - d;
//        err += y;
//    }
//    return d;
//}
//
//double normDiff(double* a, double* b, int& n)
//{
//    double d = 0.0;
//    double err = 0.0;
//    for(int i = 0; i < n; ++i)
//    {
//        // does d += abs(b[i] - a[i]) with high accuracy
//        double temp = d;
//        double y = abs(b[i] - a[i]) + err;
//        d = temp + y;
//        err = temp - d;
//        err += y;
//    }
//    return d;
//}
//
//int pageRank(string format, string inputFile)
//{
//	cout << "Format: " << format
//		<< "\nInput file: " << inputFile
//		<< endl;
//
//	auto tmStart = chrono::high_resolution_clock::now();
//	auto totalSt = tmStart;
//
//	SparseMatrix* matrix;
//	if (format.compare("CSR") == 0)
//	{
//		matrix = new SparseMatrixCSR(inputFile);
//	}
//	else if (format.compare("CSC") == 0)
//	{
//		matrix = new SparseMatrixCSC(inputFile);
//	}
//	else if (format.compare("COO") == 0)
//	{
//		matrix = new SparseMatrixCOO(inputFile);
//	}
//	else
//	{
//		cerr << "Unknown format: " << format << endl;
//		return 1;
//	}
//
//	auto tmInput = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
//	cout << "Reading input: " << tmInput << " seconds" << endl;
//	tmStart = chrono::high_resolution_clock::now();
//
//	int n = matrix->numVertices;
//	double* x = new double[n]; // pagerank
//	double* v = new double[n];
//	double* y = new double[n]; // new pagerank
//	int* outdeg = new int[n];
//	double* contr = new double[n]; // contribution for each vertex
//
//	double delta = 2.0;
//	int iter = 0;
//
//	for (int i = 0; i < n; ++i)
//	{
//		x[i] = v[i] = 1.0 / (double)(n);
//		y[i] = 0;
//		outdeg[i] = 0.0;
//	}
//
//	matrix->calculateOutDegree(outdeg);
//
//	auto tmInit = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
//	cout << "Initialisation: " << tmInit << " seconds" << endl;
//	auto tmTotal = 0.0;
//	tmStart = chrono::high_resolution_clock::now();
//
//	auto iterateS = chrono::high_resolution_clock::now();
//	while (iter < maxIter && delta > tol)
//	{
//		matrix->iterate(d, x, y, outdeg, contr);
//
//		double w = 1.0 - sum(y, n);
//		for (int i = 0; i < n; ++i)
//			y[i] += w * v[i];
//
//		delta = normDiff(x, y, n);
//		++iter;
//
//		for (int i = 0; i < n; ++i)
//		{
//			x[i] = y[i];
//			y[i] = 0.0;
//		}
//
//		auto tmStep = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tmStart).count()*1e-9;
//		cout << "iteration " << iter << ": delta=" << delta << " xnorm=" << sum(x, n)
//			<< " time=" << tmStep << " seconds" << endl;
//
//		tmStart = chrono::high_resolution_clock::now();
//	}
//	auto endT = chrono::high_resolution_clock::now();
//
//	auto iterateT = chrono::duration_cast<chrono::nanoseconds>(endT - iterateS).count()*1e-9;
//	auto totalT = chrono::duration_cast<chrono::nanoseconds>(endT - totalSt).count()*1e-9;
//
//	cout << "\nTotal time:" << totalT << " seconds\n"
//		<< "Total iterate time:" << iterateT << " seconds\n"
//		<< "Total seq time:" << (totalT - iterateT) << " seconds" << endl;
//
//	if (delta > tol)
//		cerr << "error: solution has not converged" << endl;
//
//	return 0;
//}
//
//int main(int argc, char** argv)
//{
//    cout << setprecision(16);
//	
//	int returnCode;
//	//returnCode = pageRank("COO", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\smallGraph.coo");
//	//returnCode = pageRank("COO", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\LiveJournal_dir.coo");
//	//returnCode = pageRank("COO", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\orkut_undir.coo");
//	//returnCode = pageRank("COO", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\rMatGraph_J_5_100.coo");
//
//	//returnCode = pageRank("CSR", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\smallGraph.csr");
//	returnCode = pageRank("CSR", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\LiveJournal_dir.csr");
//	//returnCode = pageRank("CSR", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\orkut_undir.csr");
//	//returnCode = pageRank("CSR", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\rMatGraph_J_5_100.csr");
//
//	//returnCode = pageRank("CSC", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\smallGraph.csc");
//	//returnCode = pageRank("CSC", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\LiveJournal_dir.csc");
//	//returnCode = pageRank("CSC", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\orkut_undir.csc");
//	//returnCode = pageRank("CSC", "D:\\Jorta\\Documents\\Uni\\4th Year\\CSC3021_Concurrent_Programming\\Assignments\\graphs\\rMatGraph_J_5_100.csc");
//
//
//    return 0;
//}
//
//

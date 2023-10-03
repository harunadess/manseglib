/*
 * Use command-line flag -ea for java VM to enable assertions.
 */

import java.lang.Long;
import java.lang.Integer;
import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.StringTokenizer;
import java.util.concurrent.*;

abstract class SparseMatrix {
    int num_vertices; // Number of vertices in the graph
    int num_edges;    // Number of edges in the graph

    // Auxiliary in preparation of PageRank iteration: pre-calculate the
    // out-degree (number of outgoing edges) for each vertex
    abstract void calculateOutDegree(int outdeg[]);

    // Perform one PageRank iteration.
    //    dampingFactor: damping factor
    //    prevPageRanks[]: previous PageRank values, read-only
    //    newPageRanks[]: new PageRank values, initialised to zero
    //    outdeg[]: values pre-calculated by calculateOutDegree()
    abstract void iterate(double dampingFactor, double[] prevPageRanks, double[] newPageRanks, int outdeg[]);
}

// This class represents the adjacency matrix of a graph as a sparse matrix
// in coordinate format (COO)
class SparseMatrixCOO extends SparseMatrix {
    private int source[];
    private int destination[];

    SparseMatrixCOO(String file) {
        try {
            InputStreamReader is
                    = new InputStreamReader(new FileInputStream(file), "UTF-8");
            BufferedReader rd = new BufferedReader(is);
            readFile(rd);
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + e);
            return;
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unsupported encoding exception: " + e);
            return;
        } catch (Exception e) {
            System.err.println("Exception: " + e);
            return;
        }
    }

    int getNext(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        return Integer.parseInt(line);
    }

    void getNextPair(BufferedReader rd, int pair[]) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        StringTokenizer st = new StringTokenizer(line);
        pair[0] = Integer.parseInt(st.nextToken());
        pair[1] = Integer.parseInt(st.nextToken());
    }

    void readFile(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        if (!line.equalsIgnoreCase("COO"))
            throw new Exception("file format error -- header");

        num_vertices = getNext(rd);
        num_edges = getNext(rd);

        source = new int[num_edges];
        destination = new int[num_edges];

        int edge[] = new int[2];
        for (int i = 0; i < num_edges; ++i) {
            getNextPair(rd, edge);
            source[i] = edge[0];
            destination[i] = edge[1];
        }
    }

    // Auxiliary function for PageRank calculation
    void calculateOutDegree(int outdeg[]) {
        for (int i = 0; i < num_edges; ++i) {
            ++outdeg[source[i]];
        }
    }

    void iterate(double dampingFactor, double[] prevPageRanks, double[] newPageRanks, int outdeg[]) {
        int src, dest;
        for (int i = 0; i < num_edges; ++i) {
            src = source[i];
            dest = destination[i];
            newPageRanks[dest] += dampingFactor * (prevPageRanks[src] / outdeg[src]);
        }
    }

    void iterateSinglePrecision(float dampingFactor, float[] prevPageRanks, float[] newPageRanks, int outdeg[]) {
    }
}

// This class represents the adjacency matrix of a graph as a sparse matrix
// in compressed sparse rows format (CSR), where a row index corresponds to
// a source vertex and a column index corresponds to a destination
class SparseMatrixCSR extends SparseMatrix {
    private int index[];
    private int destination[];

    SparseMatrixCSR(String file) {
        try {
            InputStreamReader is
                    = new InputStreamReader(new FileInputStream(file), "UTF-8");
            BufferedReader rd = new BufferedReader(is);
            readFile(rd);
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + e);
            return;
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unsupported encoding exception: " + e);
            return;
        } catch (Exception e) {
            System.err.println("Exception: " + e);
            return;
        }
    }

    int getNext(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        return Integer.parseInt(line);
    }

    void readFile(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        if (!line.equalsIgnoreCase("CSR") && !line.equalsIgnoreCase("CSC-CSR"))
            throw new Exception("file format error -- header");

        num_vertices = getNext(rd);
        num_edges = getNext(rd);

        index = new int[num_vertices + 1];
        destination = new int[num_edges];
        int pos = 0; // init index counter for destination array

        for (int i = 0; i < num_vertices; ++i) {
            line = rd.readLine();
            if (line == null)
                throw new Exception("premature end of file");
            String elm[] = line.split(" ");
            assert Integer.parseInt(elm[0]) == i : "Error in CSR file";

            pos += (elm.length - 1);
            index[i + 1] = pos; // assign next index in array for vertex edges start

            for (int j = 1; j < elm.length; ++j) {
                int dest = Integer.parseInt(elm[j]);
                destination[index[i] + (j - 1)] = dest;
            }
        }
    }

    // Auxiliary function for PageRank calculation
    void calculateOutDegree(int outdeg[]) {
        for(int i = 0; i < num_vertices; ++i) {
            outdeg[i] = (index[i + 1] - index[i]);
        }

    }

    void iterate(double dampingFactor, double[] prevPageRanks, double[] newPageRanks, int outdeg[]) {
        int currentIndex;
        int nextIndex;
        double c[] = new double[num_vertices];
        for(int i = 0; i < num_vertices; ++i) {
            c[i] = dampingFactor*(prevPageRanks[i]/outdeg[i]);
        }
        
        for (int i = 0; i < num_vertices; ++i) {
            currentIndex = index[i];
            nextIndex = index[i + 1];
            for (int j = currentIndex; j < nextIndex; ++j) {
                newPageRanks[destination[j]] += c[i];
            }
        }
    }
}

// This class represents the adjacency matrix of a graph as a sparse matrix
// in compressed sparse columns format (CSC). The incoming edges for each
// vertex are listed.
class SparseMatrixCSC extends SparseMatrix {
    private int index[];
    private int source[];

    SparseMatrixCSC(String file) {
        try {
            InputStreamReader is
                    = new InputStreamReader(new FileInputStream(file), "UTF-8");
            BufferedReader rd = new BufferedReader(is);
            readFile(rd);
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + e);
            return;
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unsupported encoding exception: " + e);
            return;
        } catch (Exception e) {
            System.err.println("Exception: " + e);
            return;
        }
    }

    int getNext(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        return Integer.parseInt(line);
    }

    void readFile(BufferedReader rd) throws Exception {
        String line = rd.readLine();
        if (line == null)
            throw new Exception("premature end of file");
        if (!line.equalsIgnoreCase("CSC") && !line.equalsIgnoreCase("CSC-CSR"))
            throw new Exception("file format error -- header");

        num_vertices = getNext(rd);
        num_edges = getNext(rd);

        index = new int[num_vertices + 1];
        source = new int[num_edges];

        int pos = 0; // init index counter for source array

        for (int i = 0; i < num_vertices; ++i) {
            line = rd.readLine();
            if (line == null)
                throw new Exception("premature end of file");
            String elm[] = line.split(" ");
            assert Integer.parseInt(elm[0]) == i : "Error in CSC file";

            pos += (elm.length - 1);
            index[i + 1] = pos;

            for (int j = 1; j < elm.length; ++j) {
                int src = Integer.parseInt(elm[j]);
                source[index[i] + (j - 1)] = src;
            }
        }
    }

    // Auxiliary function for PageRank calculation
    void calculateOutDegree(int outdeg[]) {
        for (int i = 0; i < num_edges; ++i) {
            ++outdeg[source[i]];
        }
    }

    void iterate(double dampingFactor, double[] prevPageRanks, double[] newPageRanks, int outdeg[]) {
        int currentIndex;
        int nextIndex;
        double c[] = new double[num_vertices];
        for(int i = 0; i < num_vertices; ++i)
            c[i] = dampingFactor*(prevPageRanks[i]/outdeg[i]);

        for (int i = 0; i < num_vertices; ++i) {
            currentIndex = index[i];
            nextIndex = index[i + 1];
            for (int j = currentIndex; j < nextIndex; ++j) {
                newPageRanks[i] += c[source[j]];
            }
        }
    }
}

// Main class with main() method. Performs the PageRank computation until
// convergence is reached.
class PageRank {
    public static void main(String args[]) {

        if (args.length < 2) {
            System.err.println("Usage: java pagerank format inputfile");
            return;
        }

        String format = args[0];
        String inputFile = args[1];

        // Tell us what you're doing
        System.err.println("Format: " + format);
        System.err.println("Input file: " + inputFile);

        long tm_start = System.nanoTime();
        long total_st = tm_start;

        SparseMatrix matrix;

        //find format, create specific matrix from that format
        if (format.equalsIgnoreCase("CSR")) {
            matrix = new SparseMatrixCSR(inputFile);
        } else if (format.equalsIgnoreCase("CSC")) {
            matrix = new SparseMatrixCSC(inputFile);
        } else if (format.equalsIgnoreCase("COO")) {
            matrix = new SparseMatrixCOO(inputFile);
        } else {
            System.err.println("Unknown format '" + format + "'");
            return;
        }

        double tm_input = (double) (System.nanoTime() - tm_start) * 1e-9;
        System.err.println("Reading input: " + tm_input + " seconds");
        tm_start = System.nanoTime();

        final int n = matrix.num_vertices; //total number of vertices
        double x[] = new double[n];        // pageRank
        double v[] = new double[n];
        double y[] = new double[n];        // newPageRank
        final double d = 0.85; // Leave this value as is    // dampening factor
        final double tol = 1e-7; // Leave this value as is  // failure tolerance
        final int max_iter = 100;                           // max iterations
        final boolean verbose = true;
        double delta = 2;                                   //default delta - Ll-norm?
        int iter = 0;

        //this loop is doing line 1 of the pseudocode algorithm:
        //initially fills page ranks with 1/n (number of vertices)
        for (int i = 0; i < n; ++i) {
            x[i] = v[i] = 1.0 / (double) n;
            y[i] = 0;
        }

        int outdeg[] = new int[n];
        matrix.calculateOutDegree(outdeg); //precalc outDegree of vertices in graph

        double tm_init = (double) (System.nanoTime() - tm_start) * 1e-9;
        System.err.println("Initialisation: " + tm_init + " seconds");
        tm_start = System.nanoTime();

        long iterate_s = System.nanoTime();
        while (iter < max_iter && delta > tol) {
            // Power iteration step.
            // 1. Transfering weight over out-going links (summation part)
            matrix.iterate(d, x, y, outdeg);
            // 2. Constants (1-d)v[i] added in separately.
            double w = 1.0 - sum(y, n); // ensure y[] will sum to 1
            for (int i = 0; i < n; ++i)
                y[i] += w * v[i];

            // Calculate residual error
            delta = normdiff(x, y, n);
            iter++;

            // Swap x[] and y[] and reset y[]
            for (int i = 0; i < n; ++i) {
                x[i] = y[i];
                y[i] = 0.;
            }

            double tm_step = (double) (System.nanoTime() - tm_start) * 1e-9;
            if (verbose)
                System.err.println("iteration " + iter + ": delta=" + delta
                        + " xnorm=" + sum(x, n)
                        + " time=" + tm_step + " seconds");
            tm_start = System.nanoTime();
        }
        double end_t = System.nanoTime();

        double iterate_t = (end_t - iterate_s)*1e-9;
        double total_t = (end_t - total_st)*1e-9;
        System.err.println("Total time:" + total_t + " seconds");
        System.err.println("Total iterate time:" + iterate_t + " seconds");
        System.err.println("Total seq time:" + (total_t - iterate_t) + " seconds");

        if (delta > tol)
            System.err.println("Error: solution has not converged.");

        // Dump PageRank values to file
        // for(int i = 0; i < n; ++i)
        //     System.out.println(i + " " + x[i]);
    }

    static double sum(double[] a, int n) {
        double d = 0.;
        double err = 0.;
        for (int i = 0; i < n; ++i) {
            // The code below achieves
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

    static double normdiff(double[] a, double[] b, int n) {
        double d = 0.;
        double err = 0.;
        for (int i = 0; i < n; ++i) {
            // The code below achieves
            // d += Math.abs(b[i] - a[i]);
            // but does so with high accuracy
            double tmp = d;
            double y = Math.abs(b[i] - a[i]) + err;
            d = tmp + y;
            err = tmp - d;
            err += y;
        }
        return d;
    }
}

/*// single-precision version of main() for COO, (Q3, part iv)
class PageRank {
    public static void main(String args[]) {

        if (args.length < 3) {
            System.err.println("Usage: java pagerank format inputfile outputfile");
            return;
        }

        String format = args[0];
        String inputFile = args[1];
        String outputFile = args[2];

        // Tell us what you're doing
        System.err.println("Format: " + format);
        System.err.println("Input file: " + inputFile);
        System.err.println("Output file: " + outputFile);

        long tm_start = System.nanoTime();

        SparseMatrixCOO matrix;

        //find format, create specific matrix from that format
        if (format.equalsIgnoreCase("COO")) {
            matrix = new SparseMatrixCOO(inputFile);
        } else {
            System.err.println("No single-precision version of '" + format + "'");
            return;
        }

        double tm_input = (double) (System.nanoTime() - tm_start) * 1e-9;
        System.err.println("Reading input: " + tm_input + " seconds");
        tm_start = System.nanoTime();

        final int n = matrix.num_vertices; //total number of vertices
        float x[] = new float[n];        // pageRank
        float v[] = new float[n];
        float y[] = new float[n];        // newPageRank
        final float d = 0.85f; // Leave this value as is    // dampening factor
        final float tol = 1e-7f; // Leave this value as is  // failure tolerance
        final int max_iter = 100;                           // max iterations
        final boolean verbose = true;
        float delta = 2;                                   //default delta - Ll-norm?
        int iter = 0;

        //this loop is doing line 1 of the pseudocode algorithm:
        //initially fills page ranks with 1/n (number of vertices)
        for (int i = 0; i < n; ++i) {
            x[i] = v[i] = 1.0f / (float) n;
            y[i] = 0;
        }

        int outdeg[] = new int[n];
        matrix.calculateOutDegree(outdeg); //precalc outDegree of vertices in graph

        double tm_init = (double) (System.nanoTime() - tm_start) * 1e-9;
        System.err.println("Initialisation: " + tm_init + " seconds");
        tm_start = System.nanoTime();

        while (iter < max_iter && delta > tol) {
            // Power iteration step.
            // 1. Transfering weight over out-going links (summation part)
            matrix.iterateSinglePrecision(d, x, y, outdeg);
            // 2. Constants (1-d)v[i] added in separately.
            double w = 1.0 - sum(y, n); // ensure y[] will sum to 1
            for (int i = 0; i < n; ++i)
                y[i] += w * v[i];

            // Calculate residual error
            delta = normdiff(x, y, n);
            iter++;

            // Swap x[] and y[] and reset y[]
            for (int i = 0; i < n; ++i) {
                x[i] = y[i];
                y[i] = 0.f;
            }

            double tm_step = (double) (System.nanoTime() - tm_start) * 1e-9;
            if (verbose)
                System.err.println("iteration " + iter + ": delta=" + delta
                        + " xnorm=" + sum(x, n)
                        + " time=" + tm_step + " seconds");
            tm_start = System.nanoTime();
        }

        if (delta > tol)
            System.err.println("Error: solution has not converged.");

        // Dump PageRank values to file
        writeToFile(outputFile, x, n);
    }

    static double sum(float[] a, int n) {
        double d = 0.;
        double err = 0.;
        for (int i = 0; i < n; ++i) {
            // The code below achieves
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

    static float normdiff(float[] a, float[] b, int n) {
        float d = 0.f;
        float err = 0.f;
        for (int i = 0; i < n; ++i) {
            // The code below achieves
            // d += Math.abs(b[i] - a[i]);
            // but does so with high accuracy
            float tmp = d;
            float y = Math.abs(b[i] - a[i]) + err;
            d = tmp + y;
            err = tmp - d;
            err += y;
        }
        return d;
    }

    static void writeToFile(String file, float[] v, int n) {
        try {
            OutputStreamWriter os
                    = new OutputStreamWriter(new FileOutputStream(file), "UTF-8");
            BufferedWriter wr = new BufferedWriter(os);
            writeToBuffer(wr, v, n);
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + e);
            return;
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unsupported encoding exception: " + e);
            return;
        } catch (Exception e) {
            System.err.println("Exception: " + e);
            return;
        }
    }

    static void writeToBuffer(BufferedWriter buf, float[] v, int n) {
        PrintWriter out = new PrintWriter(buf);
        for (int i = 0; i < n; ++i)
            out.println(i + " " + v[i]);
        out.close();
    }
}*/

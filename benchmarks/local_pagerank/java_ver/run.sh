#!/bin/bash

javac pagerank.java

# COO
# java -Xmx8g PageRank "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo"
# java -Xmx8g PageRank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo"
java -Xmx8g PageRank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo"
# java -Xmx8g PageRank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo"

echo ""
# CSR
# java -Xmx8g PageRank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr"
# java -Xmx8g PageRank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr"
java -Xmx8g PageRank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr"
# java -Xmx8g PageRank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr"

echo ""
# CSC
# java -Xmx8g PageRank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc"
# java -Xmx8g PageRank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc"
java -Xmx8g PageRank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc"
# java -Xmx8g PageRank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc"
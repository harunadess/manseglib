#!/bin/bash

# rm pagerank pagerank.o msa_pagerank msa_pagerank.o
# rm msa_pagerank msa_pagerank.o

# make pagerank msa_pagerank
# make msa_pagerank

## default small graphs ##
# std
# echo "std version"
# COO
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo" "./results/std_smallGraph.coo.prvals"
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo" "./results/std_orkut_undir.coo.prvals"
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo" "./results/std_LiveJournal_dir.coo.prvals"
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo" "./results/std_rMatGraph_J_5_100.coo.prvals"

# CSR
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr" "./results/std_smallGraph.csr.prvals"
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr" "./results/std_orkut_undir.csr.prvals"
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr" "./results/std_LiveJournal_dir.csr.prvals"
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr" "./results/std_rMatGraph_J_5_100.csr.prvals"

# CSC
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc" "./results/std_smallGraph.csc.prvals"
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc" "./results/std_orkut_undir.csc.prvals"
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc" "./results/std_LiveJournal_dir.csc.prvals"
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc" "./results/std_rMatGraph_J_5_100.csc.prvals"

# echo ""
echo "msa version"
# msa version
# COO
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo" "./results/msa_smallGraph.coo.prvals"
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo" "./results/msa_orkut_undir.coo.prvals"
./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo" "./results/msa_LiveJournal_dir.coo.prvals"
# ./pagerank_check default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo" "./results/msa_rMatGraph_J_5_100.coo.prvals"

# CSR
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr" "./results/msa_smallGraph.csr.prvals"
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr" "./results/msa_orkut_undir.csr.prvals"
./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr" "./results/msa_LiveJournal_dir.csr.prvals"
# ./pagerank_check default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr" "./results/msa_rMatGraph_J_5_100.csr.prvals"

# CSC
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc" "./results/msa_smallGraph.csc.prvals"
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc" "./results/msa_orkut_undir.csc.prvals"
./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc" "./results/msa_LiveJournal_dir.csc.prvals"
# ./pagerank_check default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc" "./results/msa_rMatGraph_J_5_100.csc.prvals"



## SNAP ##
# echo "std version"
# std
# COO
# ./pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# CSR
# ./pagerank SNAP CSR "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP CSR "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# CSC
# ./pagerank SNAP CSC "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP CSC "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# echo ""
# echo "msa version"
# msa
# COO
# ./msa_pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"

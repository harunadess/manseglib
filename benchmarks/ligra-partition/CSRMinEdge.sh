#!/bin/bash
# setup environment
module load compilers/gcc-4.9.0
export LD_LIBRARY_PATH="/home/jsun/swanModel/jacob/cilk-swan/lib/" 
cd /home/jsun/sjw/ligra-partition/JS-CSR
#cd /home/jsun/sjw/ligra-partition/hilbert-relabel
#cd /home/jsun/sjw/ligra-partition/coo-csr
    threads="48"
    partitions="384"
    data="/local2/jsun/relabel"
    graph="twitter_dir_MinEdges_b RMAT27_dir_MinEdges_b friendster_dir_MinEdges_b LJ_dir_MinEdges_b orkut_undir_MinEdges_b Yahoo_mem_undir_MinEdges_b USAroad_undir_MinEdges_b power_dir_MinEdges_b"
    app="BFS PageRankBit PageRank PageRankDelta BC Components BP" 
    OUTDIR=/home/jsun/paper_results/NOV/clang/NOCSC
    mkdir -p $OUTDIR
    
     for g in $graph;do
      outfile=${OUTDIR}/${g}_CSR.txt
      if [ "$g" == "Yahoo_mem_undir_MinEdges_b" ]; then
    	v="1472832"
      elif [ "$g" == "orkut_undir_MinEdges_b" ]; then
    	v="715576"
      elif [ "$g" == "twitter_dir_MinEdges_b" ]; then
    	v="38054692"
      elif [ "$g" == "friendster_dir_MinEdges_b" ]; then
    	v="103509270"
      elif [ "$g" == "RMAT27_dir_MinEdges_b" ]; then
    	v="64349498"
      elif [ "$g" == "LJ_dir_MinEdges_b" ]; then
    	v="518065"
      elif [ "$g" == "USAroad_undir_MinEdges_b" ]; then
    	v="6577110"
      elif [ "$g" == "power_dir_MinEdges_b" ]; then
    	v="62507412"
      else
    	v="100"
      fi
      for f in $app; do
       if [ "$f" == "BFS" ]; then
       	 ver="vertex"
       elif [ "$f" == "BC" ]; then
    	 ver="vertex"
       else
    	 ver="edge"
       fi
       echo Algorithm:${f} >>"$outfile" 2>&1
       	   for t in $threads;do
             export CILK_NWORKERS=${t}
             for p in $partitions;do
               echo Partitions:${p}>>"$outfile" 2>&1
               echo method:${ver}>>"$outfile" 2>&1
               LD_PRELOAD="./bin/interposer_cilk.so" numactl ./${f} -b -o -c ${p} -rounds 20 -v ${ver} -r ${v} ${data}/${g} >> "$outfile" 2>&1
             done
           done
        done
     done 

     wghdata="/local2/jsun/relabel"
     wghgraph="twitter_dir_MinEdges_wgh_b RMAT27_dir_MinEdges_wgh_b friendster_dir_MinEdges_wgh_b LJ_dir_MinEdges_wgh_b orkut_undir_MinEdges_wgh_b Yahoo_mem_undir_MinEdges_wgh_b USAroad_undir_MinEdges_wgh_b power_dir_MinEdges_wgh_b"
     appWgh="SPMV BellmanFord" 
    
     for wg in $wghgraph;do
     outfile=${OUTDIR}/${wg}_CSR.txt
     for wf in $appWgh; do
      if [ "$wg" == "Yahoo_mem_undir_MinEdges_wgh_b" ]; then
    	wv="1472832"
      elif [ "$wg" == "orkut_undir_MinEdges_wgh_b" ]; then
    	wv="715576"
      elif [ "$wg" == "twitter_dir_MinEdges_wgh_b" ]; then
    	wv="38054692"
      elif [ "$wg" == "friendster_dir_MinEdges_wgh_b" ]; then
    	wv="103509270"
      elif [ "$wg" == "RMAT27_dir_MinEdges_wgh_b" ]; then
    	wv="64349498"
      elif [ "$wg" == "LJ_dir_MinEdges_wgh_b" ]; then
    	wv="518065"
      elif [ "$wg" == "USAroad_undir_MinEdges_wgh_b" ]; then
    	wv="6577110"
      elif [ "$wg" == "power_dir_MinEdges_wgh_b" ]; then
    	wv="62507412"
      else
    	wv="100"
      fi
       if [ "$wf" == "BellmanFord" ]; then
       	 wver="vertex"
       else
    	 wver="edge"
       fi
        echo Algorithm:${wf} >>"$outfile" 2>&1
            for wt in $threads; do
            export CILK_NWORKERS=${wt}
              echo Vertex:$wv>>"$outfile" 2>&1
              echo Method:$wver>>"$outfile" 2>&1
               LD_PRELOAD="./bin/interposer_cilk.so" numactl ./${wf} -b -o -c 384 -rounds 20 -v ${wver} -r ${wv} ${wghdata}/${wg} >> "$outfile" 2>&1
         done
      done
   done



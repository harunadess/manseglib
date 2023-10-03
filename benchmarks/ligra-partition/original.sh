#!/bin/bash
# setup environment
module load compilers/gcc-4.9.0
export LD_LIBRARY_PATH="/home/jsun/swanModel/jacob/cilk-swan/lib/" 
#. /var/shared/intel/bin/iccvars.sh intel64
cd /home/jsun/sjw/ligra-partition/JS-hilbert
    threads="48"
    partitions="384"
    data="/local2/jsun/Galois"
    graph="twitter_dir_b RMAT27_dir_b friendster_dir_b LiveJournal_dir_b orkut_undir_b Yahoo_mem_undir_b USAroad_undir_b power_dir_b"
    app="BFS PageRankBit PageRank PageRankDelta BC Components BP" 
    OUTDIR=/home/jsun/paper_results/NOV/clang/
    mkdir -p $OUTDIR
    
    for g in $graph;do
      outfile=${OUTDIR}/${g}_Hilbert.txt
      if [ "$g" == "friendster_dir_b" ]; then
    	v="5000"
      elif [ "$g" == "RMAT27_dir_b" ]; then
    	v="1000"
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
               LD_PRELOAD="./bin/interposer_cilk.so" numactl ./${f} -b -c ${p} -rounds 20 -v ${ver} -r ${v} ${data}/${g} >> "$outfile" 2>&1
             done
           done
        done
     done 

     wghdata="/local2/jsun/Galois"
     wghgraph="twitter_dir_wgh_b RMAT27_dir_wgh_b friendster_dir_wgh_b LiveJournal_dir_wgh_b orkut_undir_wgh_b Yahoo_mem_undir_wgh_b USAroad_undir_wgh_b power_dir_wgh_b"
     appWgh="SPMV BellmanFord" 
    
     for wg in $wghgraph;do
     outfile=${OUTDIR}/${wg}_Hilbert.txt
     for wf in $appWgh; do
      if [ "$g" == "friendster_dir_wgh_b" ]; then
    	wv="5000"
      elif [ "$g" == "RMAT27_dir_wgh_b" ]; then
    	wv="1000"
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
               LD_PRELOAD="./bin/interposer_cilk.so" numactl ./${wf} -b -c 384 -rounds 20 -v ${wver} -r ${wv} ${wghdata}/${wg} >> "$outfile" 2>&1
         done
      done
   done



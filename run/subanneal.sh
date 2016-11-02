#!/bin/bash

DATETIME=$(date +%d-%m_%H%M%S)
RUNDIR=$HOME/Runs/N27/rho0.34_anneal
JOBNAME=rho0.34_anneal

if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi
#for taskid in $(seq 1 10);do
#sleep 8;
#echo $taskid
  qsub                                     \
    -S /bin/bash                           \
    -cwd                                   \
    -q dodeca.q@compute-0-9.local          \
    -N ${JOBNAME}                          \
	-V                                     \
    -m e                                   \
    -e ${RUNDIR}/${JOBNAME}.e    \
    -o ${RUNDIR}/${JOBNAME}.o    \
    -p -12                                 \
    -b y                                   \
	./runanneal.sh ${RUNDIR}
#done



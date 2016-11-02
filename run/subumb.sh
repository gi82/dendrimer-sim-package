#!/bin/bash

RUNDIR=$HOME/Runs/umbrella/G6/D7
left=$1
right=$2
histsize=$3
task=$4
JOBNAME="D7G6_w$1_$2_r$4"
if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi
echo "RUNDIR:$RUNDIR"
  qsub                                     \
    -S /bin/bash                           \
    -cwd                                   \
    -N ${JOBNAME}                          \
    -q dodeca.q@compute-0-8.local          \
    -V                                     \
    -m e                                   \
    -e ${RUNDIR}/${JOBNAME}.e    \
    -o ${RUNDIR}/${JOBNAME}.o    \
    -p 0                                \
    -b y                                   \
	./runumb.sh ${RUNDIR} ${left} ${right} ${histsize} $task	


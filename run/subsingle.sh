#!/bin/bash
#specify working directory and real number
ARGS=2 
#input arguments
#1: run directory
if [ $# -ne "$ARGS" ]
then
	echo "Please pass: 1:generation number 2:real number"
exit
fi
GEN=$1
REAL=$2
RUNDIR=$HOME/Runs/corrsingle/D7/G${GEN}
TASKID=$(printf "%03d" $REAL)
if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi

if [ "$HOSTNAME" == 'hal.smt' ];then
	JOBNAME=G${GEN}D7_${TASKID}
  qsub                                 \
    -S /bin/bash                       \
    -cwd                               \
    -q dodeca.q@compute-0-4.local      \
    -N ${JOBNAME}                      \
    -M ioannis.georgiou@tuwien.ac.at   \
    -m e                               \
    -e ${RUNDIR}/${JOBNAME}.e          \
    -o ${RUNDIR}/${JOBNAME}.o          \
    -p -10                             \
    -b y                               \
    ./runsingle.sh ${GEN} ${RUNDIR} ${TASKID}

else
	echo "$GEN $RUNDIR $TASKID"
	echo "Wrong host names"
fi

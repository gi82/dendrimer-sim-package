#!/bin/bash
ARGS=3
if [ $# -ne "$ARGS" ]
then
	echo "Please pass 1:generation 2:dendrimer type 3:realization num"
	exit
fi
G=$1
DENDTYPE=$2
real=$3

RUNDIR=$HOME/Runs/dom_prl_widom/G${G}/D${DENDTYPE}
if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi
if [ "$HOSTNAME" == 'liquid16' ] ;then
    	./runwidom.sh "$var" 
elif [ "$HOSTNAME" == 'hal.smt' ];then
	JOBNAME=dprl_G${G}_D${DENDTYPE}_${real}
  qsub                                \
    -S /bin/bash                      \
    -cwd                              \
    -N ${JOBNAME}                     \
    -M ioannis.georgiou@tuwien.ac.at  \
    -q dodeca.q@compute-0-6.local  \
    -m e                              \
    -e ${RUNDIR}/${JOBNAME}.e		  \
    -o ${RUNDIR}/${JOBNAME}.o    	  \
    -p -10                            \
    -b y                              \
    ./runwidom.sh ${RUNDIR} ${G} ${DENDTYPE} ${real}

else
	echo "Wrong host names"
fi

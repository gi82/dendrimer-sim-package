#!/bin/bash
#specify working directory
DENDTYPE=$1
RUNDIR=$HOME/Runs/eff_pot/fixCM/G4/longD7_5
if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi
DIST=$(seq 3.0 0.125 3.9)
TASKID=0
for var in ${DIST}
do
d=$(printf "%5.3f" $var)
echo $d
if [ "$HOSTNAME" == 'liquid16' ] ;then
    	./rundist.sh $var ${RUNDIR}
elif [ "$HOSTNAME" == 'hal.smt' ];then
	JOBNAME=D${DENDTYPE}_G4d${d}
  qsub                               \
    -S /bin/bash                     \
    -cwd                             \
    -q dodeca.q@compute-0-7.local          \
    -N ${JOBNAME}                    \
    -M ioannis.georgiou@tuwien.ac.at \
    -m e                             \
    -e ${RUNDIR}/${JOBNAME}.e        \
    -o ${RUNDIR}/${JOBNAME}.o        \
    -p -10                           \
    -b y                             \
    ./rundist.sh $d ${RUNDIR} ${DENDTYPE} ${TASKID}
sleep 10
else
	echo "Wrong host names"
fi
done

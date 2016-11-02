#!/bin/bash

DATETIME=$(date +%d-%m_%H%M%S)
NUMDENDRIMERS=220
BOXSIZE=43
RESTART=0
NUMMONOMERS=$(echo ${NUMDENDRIMERS}*62|bc)
factor=72.25657
rho=$(echo "scale=3;${NUMDENDRIMERS}*${factor}/${BOXSIZE}^3" |bc)
rho=$(printf "%5.3f" $rho)
echo $rho
RUNDIR=$HOME/Runs/N_220_D12/N_220_${BOXSIZE}
echo "Starting subscript for box size:" ${BOXSIZE}
sleep 5;

if [[ ! (-d $RUNDIR) ]];then mkdir -p "$RUNDIR";fi
for taskid in $(seq 1 1 1000);do
	echo $taskid
	INICONFIGFILE=$PWD/input.no
	sleep 2;
	./subintervsc.sh ${NUMDENDRIMERS} ${BOXSIZE} ${RESTART} ${INICONFIGFILE} ${RUNDIR} ${taskid}
	sleep 3;
done




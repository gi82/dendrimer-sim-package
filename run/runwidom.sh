#!/bin/bash 
#script requires [ARGS] parameters
ARGS=4 
#input is initial distance between dendrimers
if [ $# -ne "$ARGS" ]
then
	echo "Please pass 1:run directory 2:generation 3:dendrimer type as input 4:realnum"
exit
fi
RUNDIR=$1
G=$2
DENDTYPE=$3
REAL=$4
F=3
EXEC=widom 

NSAMPLE1=80000000
NSAMPLE2=2000
NDECOR=280
HISTSIZE=140
RMIN=0.0
RMAX=14.000

DATETIME=$(date +%d-%m_%H%M%S)

INPUTFILENAME=inputwidom.dat

SIMUL1="$RUNDIR/widomD${DENDTYPE}ns${NSAMPLE1}G${G}F${F}rmin${RMIN}_rmax${RMAX}_${REAL}"
INPUTFILE="$SIMUL1/$INPUTFILENAME"

if [[ ! (-d $SIMUL1) ]];then mkdir -p "$SIMUL1";fi
echo "input file:"$INPUTFILE
echo "simulation directory:"$SIMUL1

#***********************************#
# Write parameters to input file    #
#***********************************#
echo $DENDTYPE    >  "$INPUTFILE" 
echo $G           >> "$INPUTFILE"
echo $F           >> "$INPUTFILE"
echo $NSAMPLE1    >> "$INPUTFILE"
echo $NSAMPLE2    >> "$INPUTFILE"
echo $NDECOR      >> "$INPUTFILE"
echo $HISTSIZE    >> "$INPUTFILE"
echo $RMIN        >> "$INPUTFILE"
echo $RMAX        >> "$INPUTFILE"


# COPY EXECUTABLE 
#cd $RUNDIR/../src/
cp $EXEC $SIMUL1
# GO TO SIMULATION DIRECTORY
cd $SIMUL1
time ./$EXEC $INPUTFILE 

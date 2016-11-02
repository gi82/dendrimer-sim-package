#!/bin/bash
#script requires [ARGS] parameters
ARGS=3 
#input arguments
#1: run directory
if [ $# -ne "$ARGS" ]
then
	echo "Please pass: 1:generation number 2:rundir as arguments 3:TASKID"
exit
fi
INIDIST=0.0
RESTART=0
G=$1
RUNDIR=$2
TASKID=$3
NUMINICONFFILES=0
NUMXYZFILES=1
POTTYPE=7
SEED=0

DATETIME=$(date +%d-%m_%H%M%S)
#RUNDIR=$(pwd)

EXEC=singlenocell
INPUTFILENAME=input.dat
POTENTIALFILENAME=potential.dat

F=3
NUMDENDRIMERS=1

INICONFIGFILE=$HOME/noinput

TEMP=1.000

XBOX=200
YBOX=200
ZBOX=200

CELLLISTSIZE=1.6

MAXSTEP=0.3
DECORCYCLES=1000
PRODCYCLES=100000000
EQUILCYCLES=100000
ACCRATIO=0.5;
FRAMESFREQ=50000
SAMPLINGFREQ=50000

CC=0 
CS=1
BB=2
SS=2
# FENE type/K/L0/R
FENECC="$CC 15.0 0.0 1.5"
FENECS="$CS 15.0 0.0 1.5"
FENEBB="$BB 15.0 0.0 1.5"
# MORSE type/epsilon (in KbT)/a/d/min/max
MORSECC="$CC 1.0 20.00 0.80 0.2 1.25"
MORSECS="$CS 1.0 20.00 0.80 0.2 1.25"
MORSESS="$SS 1.0 20.00 0.80 0.2 1.25"

SIMUL1="$RUNDIR/${TASKID}"
INPUTFILE="$SIMUL1/$INPUTFILENAME"
POTENTIALFILE="$SIMUL1/$POTENTIALFILENAME"

if [[ ! (-d $SIMUL1) ]];then mkdir -p "$SIMUL1";fi
echo "input file:"$INPUTFILE
echo "simulation directory:"$SIMUL1

INICONFIGFILENAME="${INICONFIGFILE##*/}"
echo ${DATETIME}
# copy everything to simulation directory
cp $EXEC $SIMUL1
if [ $NUMINICONFFILES -eq 1 ];then 
	cp $INICONFIGFILE $SIMUL1/$INICONFIGFILENAME
fi
#***********************************#
# Write parameters to input file    #
#***********************************#
echo $TASKID                                                       >  "$INPUTFILE"
echo $G $F $NUMDENDRIMERS                                          >> "$INPUTFILE" 
echo $INIDIST                                                      >> "$INPUTFILE"
echo $SEED $NUMINICONFFILES $INICONFIGFILENAME  ${RESTART}         >> "$INPUTFILE"
echo $NUMXYZFILES												   >> "$INPUTFILE"
echo $TEMP $XBOX $YBOX $ZBOX                             		   >> "$INPUTFILE" 
echo $CELLLISTSIZE                                       		   >> "$INPUTFILE"
echo $DECORCYCLES $EQUILCYCLES $PRODCYCLES $ACCRATIO     		   >> "$INPUTFILE"
echo $SAMPLINGFREQ $FRAMESFREQ                           		   >> "$INPUTFILE"
echo $MAXSTEP                                            		   >> "$INPUTFILE"
echo $POTTYPE $POTENTIALFILENAME                         		   >> "$INPUTFILE"
if   (($POTTYPE<=0));then
echo $FENECC                                             		   >  "$POTENTIALFILE"
echo $FENECS                                             		   >> "$POTENTIALFILE"
echo $FENEBB                                             		   >> "$POTENTIALFILE"
echo $MORSECC                                            		   >> "$POTENTIALFILE"
echo $MORSECS                                            		   >> "$POTENTIALFILE"
echo $MORSESS                                            		   >> "$POTENTIALFILE"
fi

# COPY EXECUTABLE 
#cd $RUNDIR/../src/
cp $EXEC $SIMUL1
# GO TO SIMULATION DIRECTORY
cd $SIMUL1
time ./$EXEC $INPUTFILE 

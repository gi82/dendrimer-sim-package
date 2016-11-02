#!/bin/bash
#script requires [ARGS] parameters
ARGS=4 
#input arguments
#1: initial distance between dendrimers
#2: run directory
if [ $# -ne "$ARGS" ]
then
	echo "Please pass 1:initial distance 2:directory 3:potential type 4:TaskId as argument"
exit
fi
sleepparam=20
INIDIST=$1
RUNDIR=$2
POTTYPE=$3
TASKID=$4
NUMINICONFFILES=0
NUMXYZFILES=1
SEED=0
RESTART=0
DATETIME=$(date +%d-%m_%H%M%S)

EXEC=eff
INPUTFILENAME=input.dat
POTENTIALFILENAME=potential.dat

G=4
F=3
NUMDENDRIMERS=2

INICONFIGFILE="ptclsG${G}D7_"

TEMP=1.000

XBOX=70
YBOX=70
ZBOX=70

CELLLISTSIZE=1.3

MAXSTEP=0.3
DECORCYCLES=2000
EQUILCYCLES=200000
PRODCYCLES=50000000
ACCRATIO=0.5;
FRAMESFREQ=500
SAMPLINGFREQ=500

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

SIMUL1="$RUNDIR/D${POTTYPE}_G${G}F${F}/Rini${INIDIST}"
INPUTFILE="$SIMUL1/$INPUTFILENAME"
POTENTIALFILE="$SIMUL1/$POTENTIALFILENAME"

if [[ ! (-d $SIMUL1) ]];then mkdir -p "$SIMUL1";fi
echo "input file:"$INPUTFILE
echo "simulation directory:"$SIMUL1

#***********************************#
# Write parameters to input file    #
#***********************************#
echo $TASKID                                             	>  "$INPUTFILE"
echo $G $F $NUMDENDRIMERS                                	>> "$INPUTFILE" 
echo $INIDIST                                            	>> "$INPUTFILE"
echo $SEED $NUMINICONFFILES $RUNDIR/$INICONFIGFILE $RESTART >> "$INPUTFILE"
echo $NUMXYZFILES                                        	>> "$INPUTFILE"
echo $TEMP $XBOX $YBOX $ZBOX                             	>> "$INPUTFILE" 
echo $CELLLISTSIZE                                       	>> "$INPUTFILE"
echo $DECORCYCLES $EQUILCYCLES $PRODCYCLES $ACCRATIO     	>> "$INPUTFILE"
echo $SAMPLINGFREQ $FRAMESFREQ                           	>> "$INPUTFILE"
echo $MAXSTEP                                            	>> "$INPUTFILE"
echo $POTTYPE $POTENTIALFILENAME                         	>> "$INPUTFILE"

if   (($POTTYPE<=0));then
echo $FENECC                                             	>  "$POTENTIALFILE"
echo $FENECS                                             	>> "$POTENTIALFILE"
echo $FENEBB                                             	>> "$POTENTIALFILE"
echo $MORSECC                                            	>> "$POTENTIALFILE"
echo $MORSECS                                            	>> "$POTENTIALFILE"
echo $MORSESS                                            	>> "$POTENTIALFILE"
fi

# COPY EXECUTABLE 
#cd $RUNDIR/../src/
cp $EXEC $SIMUL1
# GO TO SIMULATION DIRECTORY
cd $SIMUL1
sleep $sleepparam
time ./$EXEC $INPUTFILE 

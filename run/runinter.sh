#!/bin/bash
#script requires [ARGS] parameters
ARGS=6 
#input arguments
#1: initial distance between dendrimers
#2: run directory
if [ $# -ne "$ARGS" ]
then
	echo "Please pass 1:num of dendrimers 2: boxsize 3: restart (0 or 1) 4:initial configuration file 5:full path of run directory 6:simulation number"
	exit;
fi
#initial distance set to zero
INIDIST=0.0
NUMDENDRIMERS=$1
NUMMONOMERS=$(echo ${NUMDENDRIMERS}*62|bc)
BOXSIZE=$2
# restart recommended when having an initial configuration file
# in the correct box size
RESTART=$3
INICONFIGFILE=$4
RUNDIR=$5
TASKID=$6
TID=$(printf "%03d" $TASKID)
echo $TID
EXEC=inter
BOOLINICONFFILES=0
# if restart the ini configuration files have to be provided
if [ $RESTART -eq 1 ]; then
BOOLINICONFFILES=1
fi

NUMXYZFILES=1
POTTYPE=7
SEED=0
#calculate rho
factor=72.25657
rho=$(echo "scale=3;${NUMDENDRIMERS}*${factor}/${BOXSIZE}^3" |bc)
rho=$(printf "%5.3f" $rho)
echo rho=$rho
DATETIME=$(date +%d-%m_%H%M%S)

INPUTFILENAME=input.dat
POTENTIALFILENAME=potential.dat

G=4
F=3

TEMP=1.00

XBOX=${BOXSIZE}
YBOX=${BOXSIZE}
ZBOX=${BOXSIZE}

CELLLISTSIZE=1.0

MAXSTEP=0.30
DECORCYCLES=1000
EQUILCYCLES=110000
PRODCYCLES=10000000
ACCRATIO=0.5;
FRAMESFREQ=100000
SAMPLINGFREQ=10000

CC=0 
CS=1
BB=2
SS=2

#FENECC="$CC 40.0 1.8750 0.3750"
#FENECS="$CS 20.0 2.8125 0.5625"
#FENEBB="$BB 40.0 1.8750 0.3750"
FENECC="$CC 20.0 0.7 0.3"
FENECS="$CS 20.0 0.7 0.3"
FENEBB="$BB 20.0 0.7 0.3"

MORSECC="$CC 1.0 24.00 0.80 0.2 1.25"
MORSECS="$CS 1.0 24.00 0.80 0.2 1.25"
MORSESS="$SS 1.0 24.00 0.80 0.2 1.25"

SIMUL1="$RUNDIR/N_${NUMDENDRIMERS}_rho_${rho}_real_${TID}"
INPUTFILE="$SIMUL1/$INPUTFILENAME"
POTENTIALFILE="$SIMUL1/$POTENTIALFILENAME"

if [[ ! (-d $SIMUL1) ]];then mkdir -p "$SIMUL1";fi
#extract cinfig filename from given path
INICONFIGFILENAME="${INICONFIGFILE##*/}"
echo ${DATETIME}
# copy everything to simulation directory
cp $EXEC $SIMUL1
if [ $BOOLINICONFFILES -eq 1 ];then 
	cp $INICONFIGFILE $SIMUL1/$INICONFIGFILENAME
fi
# GO TO SIMULATION DIRECTORY
cd $SIMUL1
echo "input file:"$INPUTFILE
echo "simulation directory:"$SIMUL1

#***********************************#
# Write parameters to input file    #
#***********************************#
echo $TASKID                                                    >  "$INPUTFILE"
echo $G $F $NUMDENDRIMERS                                       >> "$INPUTFILE" 
echo $INIDIST                                                   >> "$INPUTFILE"
echo $SEED $BOOLINICONFFILES $INICONFIGFILENAME ${RESTART}      >> "$INPUTFILE"
echo $NUMXYZFILES                                               >> "$INPUTFILE"
echo $TEMP $XBOX $YBOX $ZBOX                                    >> "$INPUTFILE" 
echo $CELLLISTSIZE                                              >> "$INPUTFILE"
echo $DECORCYCLES $EQUILCYCLES $PRODCYCLES $ACCRATIO            >> "$INPUTFILE"
echo $SAMPLINGFREQ $FRAMESFREQ                                  >> "$INPUTFILE"
echo $MAXSTEP                                                   >> "$INPUTFILE"
echo $POTTYPE $POTENTIALFILENAME                                >> "$INPUTFILE"
if   (($POTTYPE<=0));then
echo $FENECC                                                    >  "$POTENTIALFILE"
echo $FENECS                                                    >> "$POTENTIALFILE"
echo $FENEBB                                                    >> "$POTENTIALFILE"
echo $MORSECC                                                   >> "$POTENTIALFILE"
echo $MORSECS                                                   >> "$POTENTIALFILE"
echo $MORSESS                                                   >> "$POTENTIALFILE"
fi

time ./$EXEC $INPUTFILE ${RESTART}

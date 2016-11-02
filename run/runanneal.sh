#!/bin/bash
#script requires [ARGS] parameters
ARGS=1 
#input arguments
#1: initial distance between dendrimers
#2: run directory
#sleep 5
if [ $# -ne "$ARGS" ]
then
	echo "Please pass 1:full path of run directory "
exit
fi
#initial distance set to zero
NUMDENDRIMERS=200
BOXSIZE=110
# restart recommended when having an initial configuration file
# in the correct box size
RESTART=0
INICONFIGFILE="input.0"
RUNDIR=$1
TEMPHIGH=3.0
TEMPLOW=1.0
TEMPSTEP=0.05
TEMP=$TEMPHIGH
EXEC=anneal
BOOLINICONFFILES=0
# if restart the ini configuration files have to be provided
if [ $RESTART -eq 1 ]; then
BOOLINICONFFILES=1
fi

POTTYPE=7
SEED=0

DATETIME=$(date +%d-%m_%H%M%S)

INPUTFILENAME=input.dat
POTENTIALFILENAME=potential.dat

G=4
F=3

XBOX=${BOXSIZE}
YBOX=${BOXSIZE}
ZBOX=${BOXSIZE}

CELLLISTSIZE=1.35

MAXSTEP=0.30
DECORCYCLES=1000
EQUILCYCLES=200000
PRODCYCLES=1000000
ACCRATIO=0.5;
FRAMESFREQ=1000
SAMPLINGFREQ=1000

CC=0 
CS=1
BB=2
SS=2

FENECC="$CC 20.0 0.7 0.3"
FENECS="$CS 20.0 0.7 0.3"
FENEBB="$BB 20.0 0.7 0.3"

MORSECC="$CC 1.0 24.00 0.80 0.2 1.25"
MORSECS="$CS 1.0 24.00 0.80 0.2 1.25"
MORSESS="$SS 1.0 24.00 0.80 0.2 1.25"

SIMUL1="$RUNDIR/anneal${NUMDENDRIMERS}_${XBOX}_D${POTTYPE}_G${G}F${F}_${DATETIME}"
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
echo $G $F $NUMDENDRIMERS                                       >  "$INPUTFILE" 
echo $SEED $BOOLINICONFFILES $INICONFIGFILENAME ${RESTART}      >> "$INPUTFILE"
echo $TEMPHIGH $TEMPLOW $TEMPSTEP                               >> "$INPUTFILE" 
echo $TEMP                                                      >> "$INPUTFILE"
echo $XBOX $YBOX $ZBOX                                          >> "$INPUTFILE"
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

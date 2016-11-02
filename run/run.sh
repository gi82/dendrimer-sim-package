#!/bin/bash

INIDIST=$1;
DATETIME=$(date +%d-%m_%H%M%S)
RUNDIR=$(pwd)
EXEC=inter
INPUTFILENAME=input.dat

G=2
F=3
#INIDIST=6.0
#INICONFIGFILE="MC_G${G}F${F}ini"
INICONFIGFILE="ptclsG${G}F${F}D7_"
NUMDENDRIMERS=1
NUMINICONFFILES=0
SEED=0
TEMP=1.00

XBOX=1000
YBOX=1000
ZBOX=1000

CC=0 
CS=1
BB=2
SS=2

FENECC="$CC 40.0 1.8750 0.3750"
FENECS="$CS 20.0 2.8125 0.5625" 
FENEBB="$BB 40.0 1.8750 0.3750"
  
MORSECC="$CC 0.714 6.400 1.00 0.4 2.45"
MORSECS="$CS 0.014 19.20 1.25 0.85 2.41"
MORSESS="$SS 0.014 19.20 1.50 1.00 2.42"

MAXSTEP=0.3
RLISTPERCENT=0.35
PRODCYCLES=50000000
EQUILCYCLES=100000
ACCRATIO=0.5;
FRAMESFREQ=1000
SAMPLINGFREQ=1000

#SIMUL1="$RUNDIR/try_free/dend_n${NUMDENDRIMERS}G${G}F${F}_${DATETIME}"
SIMUL1="$RUNDIR/try_freeG${G}F${F}/Rini${INIDIST}"
INPUTFILE="$SIMUL1/$INPUTFILENAME"
#setenv SIMUL1  "$RUNDIR"
if [[ ! (-d $SIMUL1) ]];then mkdir -p "$SIMUL1";fi
echo "input file:"$INPUTFILE
echo "simulation directory:"$SIMUL1

#***********************************#
# Write parameters to input file    #
#***********************************#
echo $G $F $NUMDENDRIMERS                >  "$INPUTFILE" 
echo $INIDIST                            >> "$INPUTFILE"
echo $SEED $NUMINICONFFILES              >> "$INPUTFILE"
echo $TEMP $XBOX $YBOX $ZBOX             >> "$INPUTFILE" 
echo $EQUILCYCLES $PRODCYCLES $ACCRATIO  >> "$INPUTFILE"
echo $SAMPLINGFREQ $FRAMESFREQ           >> "$INPUTFILE"
echo $MAXSTEP $RLISTPERCENT              >> "$INPUTFILE"
echo $FENECC                             >> "$INPUTFILE"
echo $FENECS                             >> "$INPUTFILE"
echo $FENEBB                             >> "$INPUTFILE"
echo $MORSECC                            >> "$INPUTFILE"
echo $MORSECS                            >> "$INPUTFILE"
echo $MORSESS                            >> "$INPUTFILE"
echo $RUNDIR/$INICONFIGFILE              >> "$INPUTFILE"

# COPY EXECUTABLE 
#cd $RUNDIR/../src/
cp $EXEC $SIMUL1
# GO TO SIMULATION DIRECTORY
cd $SIMUL1
time ./$EXEC $INPUTFILE

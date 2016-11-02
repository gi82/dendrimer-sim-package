#!/bin/bash
#$ -S /bin/bash                     
#$ -d $HOME/Runs/aaa
#$ -cwd                             
#$ -q dodeca.q@compute-0-6.local    
#$ -t 1-3:1
#$ -N name.$TASK_ID                    
#$ -V                               
#$ -M ioannis.georgiou@tuwien.ac.at 
#$ -m e                             
#$ -b y                             

#./runinter.sh ${NUMDENDRIMERS} ${BOXSIZE} ${RESTART} ${RUNDIR} ${SIMNUM}

pwd
ls
sleep10


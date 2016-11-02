#!/bin/bash
args=5
if [ $# -ne "$args" ]
then
	echo "Please pass 1:full path of run directory 2:left boundary 3:right boundary 4:histsize 5:real"
exit
fi
rundir=$1
left=$2	
right=$3
histsize=$4
real=$(echo "$5" | awk '{printf "%02s\n",$0}') 
outfile="w${left}_${right}_${histsize}.out"
box=200
pottype=7
g=4
f=3
decorrelation=20000
equilibration=200000
production=50000000000
seed=0
maxstep=0.5
outfreq=10000
checkfreq=10000000
samplefreq=500000

datetime=$(date +%d-%m_%H%M%S)

exec=umbrella

simuldir="$rundir/${real}_umb_D${pottype}_G${g}F${f}_w${left}_${right}"
inputfilename="$simuldir/in_umbrella.in"

if [[ ! (-d $simuldir) ]];then mkdir -p "$simuldir";fi
echo "input file:"$inputfilename
echo "simulation directory:"$simuldir
# copy everything to simulation directory
cp $exec $simuldir
#***********************************#
# Write parameters to input file    #
#***********************************#
echo $box                         >  "$inputfilename"
echo $left                        >> "$inputfilename"
echo $right                       >> "$inputfilename" 
echo $histsize                    >> "$inputfilename"
echo $outfile                     >> "$inputfilename"
echo $pottype					  >> "$inputfilename"
echo $g             	   	  	  >> "$inputfilename" 
echo $f             	   	  	  >> "$inputfilename"
echo $decorrelation 	   	  	  >> "$inputfilename"
echo $equilibration 	   	  	  >> "$inputfilename"
echo $production    	   	  	  >> "$inputfilename"
echo $seed          	   	  	  >> "$inputfilename"
echo $maxstep       	   	  	  >> "$inputfilename" 
echo $outfreq       	   	  	  >> "$inputfilename"
echo $checkfreq     	   	  	  >> "$inputfilename"
echo $samplefreq    	   	  	  >> "$inputfilename"

# GO TO SIMULATION DIRECTORY
cd $simuldir
time ./$exec $inputfilename

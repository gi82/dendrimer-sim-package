#!/bin/bash
for i in $(seq 0.0 0.2 1.2)
do 
	j=$(echo $i+0.3|bc)
	b=$(printf "%4.1f" $i)
	e=$(printf "%4.1f" $j)
	./subumb.sh $b $e 100 0
	sleep 20;
done

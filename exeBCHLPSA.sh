#!/bin/bash
for i in $(find ${k}*"txt"); do
	./BCHLPSA ${i} R_${i}_${heur} ${heur} >Logs/O_${i}_${heur} &
done

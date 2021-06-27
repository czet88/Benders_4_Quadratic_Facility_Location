#!/bin/bash
for i in $(find ./Inputs -name "in_"${k}*"txt"); do
	file=$(basename ${i})
	./BCHLPSA ${file} R_${file} >./Logs/O_${file} &
done

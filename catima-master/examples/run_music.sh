#!/bin/bash

for ((j=0; j<20; j++)); do 
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo "Begin Anode ${j}" 
	echo " %%%%%%%%%%%%%%%%%%%%%%%%%"
	./music_dedx ${j} 2184
done

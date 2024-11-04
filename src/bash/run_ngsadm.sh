#!/bin/bash

# set env

WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen

# load module
module purge
module load tools
module load ngsadmix/32

# enter WDIR
cd $WDIR'/7-ngsadm'
IFS=$'\n'

if ! test -d rep; then
	mkdir rep
fi
cd rep


for i in {2..6};
do
	
	for rep in {1..20}
	do

		RUN=0	
		while (($RUN == 0));
		do
			PRUN=$(ps -ef | grep xinsun | grep NGSadmix | wc -l)
			
			if (($PRUN <= 4)); then
				NGSadmix -likes $WDIR/2-gl/gl_tv_maf05_mis50.qc.rel.per10.beagle.gz -outfiles gl_tv_maf10_mis50.qc.rel.per10.${i}.${rep} -minMaf 0.10 -minInd 30 -P 10 -K $i &
				RUN=1
			else
				sleep 10
			fi
		done
	done
done 
wait


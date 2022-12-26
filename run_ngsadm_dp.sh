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
	
	for rep in {1..100}
	do

		RUN=0	
		while (($RUN == 0));
		do
			PRUN=$(ps -ef | grep xinsun | grep NGSadmix | wc -l)
			
			if (($PRUN <= 4)); then
				NGSadmix -likes $WDIR/2-gl/gl_tv_maf05_mis50.dp.per10.beagle.gz -outfiles gl_tv_maf05_mis80.dp.per10.${i}.${rep} -minMaf 0.05 -minInd 20 -P 10 -K $i &
				RUN=1
			else
				sleep 10
			fi
		done
	done
done 
wait


# get best
WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/7-ngsadm

cd $WDIR

if ! test -d best; then
        mkdir best
fi


for i in {2..6}
do
        NAME=gl_tv_maf05_mis80.dp.per10
        K=$i

        best_l=$(for rep in {1..100};do L=$(tail -1 ./rep/${NAME}.${K}.${rep}.log | awk '{split($2,a,"="); print a[2]}'); echo $rep $L; done | sort -k 2 -r -n | head -1 )
        best=$(echo $best_l | awk '{print $1}')
        best_like=$(echo $best_l | awk '{print $2}')

        cp ./rep/${NAME}.${K}.${best}.qopt ./best/${NAME}.${K}.Q
        echo $K $best $best_like >> ./best/best_summary.${NAME}

        # output detail
        for rep in {1..100};do L=$(tail -1 ./rep/${NAME}.${K}.${rep}.log | awk '{split($2,a,"="); print a[2]}'); echo $rep $L; done | sort -k 2 -r -n > ./best/best_detail.${NAME}.${K}

done


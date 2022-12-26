#!/bin/bash

# set env

WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen
LIST_NAME_BAM=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_name_bam
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta

# load module
module purge
module load tools
module load htslib/1.8
module load angsd/0.931


# enter WDIR
cd $WDIR/8-het
if ! test -d sfs_by_sample;then
	mkdir sfs_by_sample
fi
cd sfs_by_sample
	

IFS=$'\n'

# # check bam index
# # check input bam file index
# for i in $(less $LIST_NAME_BAM)
# do
#         ID=$(echo $i | awk '{print $1}')
#         BAM=$(echo $i | awk '{print $2}')
#         BAM_INDEX_1=$BAM".bai"
#         BAM_INDEX_2=$(echo "$BAM" | sed -s 's/\.bam/\.bai/')
# 
#         if test -f $BAM_INDEX_1; then
#                 echo "bam index OK"
#         else
#                 if test -f $BAM_INDEX_2; then
#                         echo "bam index OK"
#                 else
#                         module load samtools/1.10
#                         samtools index $BAM &
#                 fi
#         fi
# done
# wait
# 
# 
# # gen saf
# for i in $(less $LIST_NAME_BAM);
# do
# 	ID=$(echo $i | awk '{print $1}')
#         BAM=$(echo $i | awk '{print $2}')
# 
# 	if ! test -d $ID; then
# 		mkdir $ID
# 	fi
# 	cd $ID	
# 
# 	RUN=0	
# 	while (($RUN == 0));
# 	do
# 		PRUN=$(ps -ef | grep xinsun | grep angsd | wc -l)
# 		
# 		if (($PRUN <= 15)); then
# 			angsd -GL 2 -C 50 -ref $REF -anc $REF -out $ID -nThreads 5 -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -i $BAM -dosaf 1 -fold 1  &
# 			RUN=1
# 		else
# 			sleep 1
# 		fi
# 	done
# 	cd ..
# done 
# wait


# get sfs
# switch to fatnode

for i in $(less $LIST_NAME_BAM);
do
	ID=$(echo $i | awk '{print $1}')
        BAM=$(echo $i | awk '{print $2}')

	if ! test -d $ID; then
		mkdir $ID
	fi
	cd $ID	

	RUN=0	
	while (($RUN == 0));
	do
		PRUN=$(ps -ef | grep xinsun | grep realSFS | wc -l)
		
		if (($PRUN <= 10)); then
			{
			realSFS -P 4 ${ID}.saf.idx > ${ID}.ml
			awk '{print $2/($1+$2)}' ${ID}.ml > ${ID}.het
			} &
			RUN=1
		else
			sleep 1
		fi
	done
	cd ..
done 
wait




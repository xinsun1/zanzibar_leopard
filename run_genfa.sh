#!/bin/bash

# set env

WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen

# load module
module purge
module load tools
module load samtools/1.13
module load bcftools/1.12



# enter WDIR
cd $WDIR'/6-psmc'
IFS=$'\n'


#  template check dir
# 	if test -f $MSA_FA ; then
# 		rm $MSA_FA
# 	fi

META_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/6-psmc/list_psmc_meta
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta
PSMC=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/program/psmc-master


for i in $(less $META_LIST)
do
	ID=$(echo $i | awk '{print $1}')
	DP=$(echo $i | awk '{print $2}')
	BAM=$(echo $i | awk '{print $3}')


	RUN=0	
	while (($RUN == 0));
	do
		PRUN=$(ps -ef | grep xinsun | grep samtools | wc -l)
		echo $PRUN
		
		if (($PRUN <= 20)); then

			
			cd $WDIR/6-psmc
			#### 0. gen_fa
			
			if ! test -d $ID; then
			        mkdir $ID
			fi
			cd $ID
			
			
			# change DP filter
			samtools mpileup -C 50 -q 30 -Q 30 -uf $REF $BAM | bcftools call -c - | vcfutils.pl vcf2fq -d $(awk '{a='$DP'*0.5;if(a<=5){print 5}else{print a}}') -D $(awk '{print 2*'$DP'}') > $ID.fq &
			RUN=1
		else
			sleep 10
		fi
	done

done
wait
	



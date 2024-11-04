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


if ! test -d bootstrap; then
	mkdir bootstrap
fi

if ! test -d bs_plot; then
	mkdir bs_plot
fi


#  template check dir
# 	if test -f $MSA_FA ; then
# 		rm $MSA_FA
# 	fi

META_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/6-psmc/list_psmc_meta
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta
PSMC=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/program/psmc-master


for i in $(less $META_LIST);
do
	ID=$(echo $i | awk '{print $1}')
	DP=$(echo $i | awk '{print $2}')
	BAM=$(echo $i | awk '{print $3}')


	RUN=0	
	while (($RUN == 0));
	do
		PRUN=$(ps -ef | grep xinsun | grep psmc | wc -l)
		
		if (($PRUN <= 60)); then

			{
			cd $WDIR/6-psmc
			#### 0. gen_fa
			
			# split by chr next time!!!
			
			if ! test -d $ID; then
			        mkdir $ID
			fi
			cd $ID
			
			if ! test -f ${ID}.fq.fai; then
		
				DP_LOW=$(echo $DP | awk '{a='$DP'*0.5;if(a<=5.0){print 5.0}else{print a}}')
				DP_HIGH=$(echo $DP | awk '{print 2.0*'$DP'}')
				
				# change DP filter
				# samtools mpileup -C 50 -q 30 -Q 30 -uf $REF $BAM | bcftools call -c - | vcfutils.pl vcf2fq -d $DP_LOW -D $DP_HIGH > $ID.fq
	
				# gzip -c $ID.fq > $ID.fq.gz &
				
				CHR_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_chr_felcat9_noMT_X
				# samtools fqidx -r $CHR_LIST $ID.fq | gzip > $ID.noMTX.fq.gz
			fi

		
			#### 1. run psmc
			RID=${ID}.noMTX
			
			# generate PSMC fasta
			# $PSMC/utils/fq2psmcfa -q 20 -v ${RID}.fq > ${RID}.tv.psmcfa
			# $PSMC/utils/splitfa ${RID}.tv.psmcfa > ${RID}_split.tv.psmcfa
			
			# run PSMC
			# $PSMC/psmc -N25 -t20 -p "4+25*2+4+6" -o ${RID}.tv.psmc ${RID}.tv.psmcfa
			 
			# #### 2. run bootstrap
			# for i in {1..100}
			# do
			#         RUN=0
			#         while  (($RUN == 0));
			#         do
			#                 PRUN=$(ps -ef |grep xins | grep psmc | wc -l)
			# 
			#                 if (($PRUN <= 40)); then
			#                         $PSMC/psmc -N25 -t20 -p "4+25*2+4+6" -b -o ../bootstrap/${RID}_${i}.psmc ${RID}_split.psmcfa &
			#                         sleep 0.1
			#                         RUN=1
			#                 else
			#                         sleep 10
			#                 fi
			#         done
			# done
			# wait
			# 
			
			#### 3. prep plot
			
			cat $RID.tv.psmc ../bootstrap/${RID}_*.tv.psmc > ./$RID.withbp.tv.psmc
			
			cd $WDIR/6-psmc/bs_plot
			
			$PSMC/utils/psmc_plot.pl -g7.5 -u3e-9 -R -x100 plot_${ID}_withbp.tv $WDIR/6-psmc/$ID/$RID.withbp.tv.psmc
			
			awk '{print $1,$2,"m", 0, "'$ID'"}' OFS='\t' plot_${ID}_withbp.tv.0.txt > $WDIR/6-psmc/$ID/plot_m_bp_${ID}.tv.txt
			
			for i in {1..10}
			do
			        awk '{print $1,$2,"bp", "'$i'", "'$ID'"}' OFS='\t' plot_${ID}_withbp.tv.$i.txt >> $WDIR/6-psmc/$ID/plot_m_bp_${ID}.tv.txt
			done
			} &
			
			RUN=1
		else
			sleep 10
		fi
	done
done 
wait


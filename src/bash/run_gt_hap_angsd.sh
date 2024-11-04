#!/bin/bash

# set working folder
WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/12-eig
ID_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_name
BAM_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_bam
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta
HAP_FOLDER=/home/projects/ku-cbd/data/zanzibar/pop_gen/0-angsd_fa
REF_GT=/home/projects/ku-cbd/data/zanzibar/pop_gen/0-angsd_fa/felcat9
HAP_CHR=/home/projects/ku-cbd/data/zanzibar/pop_gen/12-eig/hap_chr
CHR_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_chr_felcat9_noMT_X


#### check and generate consensus fasta
# load angsd before use
module purge
module load tools
module load htslib/1.8 angsd/0.931

IFS=$'\n'
HAP_SUFFIX='.minD3.q30.cons'

cd $HAP_FOLDER

for i in $(paste $ID_LIST $BAM_LIST)
do
	ID=$(echo $i | awk '{print $1}')
	BAM=$(echo $i | awk '{print $2}')
	
	 # test if file exist
        if test -d $ID; then
                cd $ID
        else
                mkdir $ID
                cd $ID
        fi

        # test if old fa exist
        if ! test -f ${ID}${HAP_SUFFIX}.arg; then
		angsd -nThreads 40 -i $BAM -dofasta 2 -doCounts 1 -minQ 20 -minmapq 30 \
			-setminDepthInd 3 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
			-ref $REF -out ${ID}${HAP_SUFFIX} &
	fi
	cd ..
done
wait



#### get gt per sample per base
# load samtools
module purge
module load tools
module load samtools/1.14

# check fa index
for i in $(less $ID_LIST)
do
        if ! test -f $HAP_FOLDER/$i/${i}${HAP_SUFFIX}.fa.gz.fai; then
                samtools faidx $HAP_FOLDER/$i/${i}${HAP_SUFFIX}.fa.gz &
        fi
done

wait

# get gt
cd $WDIR

if ! test -d gt_angsd_all; then
	mkdir gt_angsd_all
fi
cd gt_angsd_all

for i in $(less $ID_LIST)
do
	# check if subfolder exist
	if ! test -d $i; then
		mkdir $i
	fi
	cd $i

	for j in $(less $CHR_LIST)
	do
		RUN=0
		while (($RUN == 0));
		do
			PRUN=$(ps -ef | grep xinsun | grep samtools | wc -l)

			if (($PRUN <=40)); then
				samtools faidx $HAP_FOLDER/$i/${i}${HAP_SUFFIX}.fa.gz $j | tail -n +2 | awk -F "" '{for(i=1;i<=NF;i++){print $i}}' > ${i}_$j &
				RUN=1
			else
				sleep 0.5
			fi
		done
	done
	cd ..
done

wait

# prep ref gt
cd $HAP_FOLDER
if ! test -d felcat9; then
	mkdir felcat9
fi
cd felcat9

for chr in $(less $CHR_LIST)
do
	samtools faidx $REF $chr | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' | tail -n +2 | awk -F "" '{for(i=1;i<=NF;i++){print $i}}' > felcat9_$chr &
done
wait



# merge by chr
# change folder name and position before use
cd $WDIR
cd gt_angsd_all

BATCH='hap.angsd.all'

for C in $(less $CHR_LIST)
do
	for i in $(less $ID_LIST);do printf $i/${i}_${C}" " ;done | xargs paste -d "" $REF_GT/felcat9_$C | awk -F "" '{delete a;n_allele=0;alt=""; for(i=1;i<=NF;i++){a[$i]++}; for(j in a){if(j != "N"){n_allele++;if(j != $1){alt=j}}}; if(n_allele==2){print "'${C}':"NR,"'${C}'","0.0",NR,$1,alt >> "'$HAP_CHR'/'$C'.snp"; for(k=2;k<=NF;k++){if($k=="N"){printf 9}else{if($k==$1){printf 2}else{if($k==alt){printf 0}else{printf 9}}}}; printf "\n"} }' OFS='\t' > $HAP_CHR/${C}.geno &

done
wait

cd $HAP_CHR
for i in $(less $CHR_LIST);do cat $i.geno >> $WDIR/${BATCH}.geno; cat $i.snp >> $WDIR/${BATCH}.snp ;done

# clean and filter
cd $WDIR

awk '{if(($5=="A" && $6=="G") || ($5=="T" && $6=="C") || ($5=="C" && $6=="T") || ($5=="G" && $6=="A")){$7="1"}else{$7="0"}; print $7}' ${BATCH}.snp > ${BATCH}.snp_is_ts

paste ${BATCH}.geno ${BATCH}.snp_is_ts | awk '{if($2==0){print $1}}' > ${BATCH}_tv.geno &
paste ${BATCH}.snp ${BATCH}.snp_is_ts | awk '{if($7==0){print $1,$2,$3,$4,$5,$6}}' OFS='\t' > ${BATCH}_tv.snp &

wait

awk -F "" '{N=0;M=0;for(i=1;i<=NF;i++){if($i=="9"){N+=1}else{M+=$i}};if(NF==N){MAF=0}else{MAF=M/(NF*2-N*2)};print N/NF, MAF}' ${BATCH}_tv.geno > ${BATCH}_tv.maf_mis

# change maf nad mis filter before use, default maf05 mis20
paste ${BATCH}_tv.snp ${BATCH}_tv.maf_mis | awk '{if($7<=0.5 && $8>=0.05 && $8<=0.95){print $1,$2,$3,$4,$5,$6}}' OFS='\t' > ${BATCH}_tv_maf05_mis50.snp &
paste ${BATCH}_tv.geno ${BATCH}_tv.maf_mis | awk '{if($2<=0.5 && $3>=0.05 && $3<=0.95){print $1}}' > ${BATCH}_tv_maf05_mis50.geno &
wait


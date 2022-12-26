#!/bin/bash

#### set working directory
WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl
GL_FOLDER=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl/gl_chr
BAM_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_bam
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta

cd $WDIR


#### GL call by chr
# load angsd before use
module purge
module load tools
module load htslib/1.8 angsd/0.931

# use checkBamHeaders 0 to skip header check due to some files lacks MT
LIST_CHR=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_chr_felcat9_noMT_angsd

if ! test -d gl_chr; then
	mkdir gl_chr
fi

for i in $(less $LIST_CHR)
do
	angsd -GL 2 -out ./gl_chr/gl_chr$i -ref $REF -nThreads 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -bam $BAM_LIST -r $i &
done

wait


#### filter and merge
# change sample size NSAMPLE
# change filter options, default maf05, mis50

cd $WDIR
cd $GL_FOLDER

NSAMPLE=80

for i in $(less $LIST_CHR)
do
	# check maf file if -ref use $6,$8
	{
	zcat gl_chr${i}.mafs.gz | awk '{if(NR==1){print 1}else{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A") || $6 <0.05 || $6 >0.95 || $8/'$NSAMPLE' < 0.5){print 0}else{print 1} }}' > gl_chr${i}.is_tv_maf05_mis50
	zcat gl_chr${i}.beagle.gz | awk 'NR==FNR {a[FNR]=$1;next}{if(a[FNR]==1){print $0}}' gl_chr${i}.is_tv_maf05_mis50 - > gl_chr${i}_tv_maf05_mis50.beagle 
	zcat gl_chr${i}.mafs.gz | awk 'NR==FNR {a[FNR]=$1;next}{if(a[FNR]==1){print $0}}' gl_chr${i}.is_tv_maf05_mis50 - > gl_chr${i}_tv_maf05_mis50.mafs 
	} &
done
wait


LIST_CHR_NOX=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_chr_felcat9_noMT_X_angsd

head -1 gl_chrchrA1:_tv_maf05_mis50.beagle >> $WDIR/gl_tv_maf05_mis50.beagle
head -1 gl_chrchrA1:_tv_maf05_mis50.mafs >> $WDIR/gl_tv_maf05_mis50.mafs
for i in $(less $LIST_CHR_NOX)
do
	tail -n +2 gl_chr${i}_tv_maf05_mis50.beagle >> $WDIR/gl_tv_maf05_mis50.beagle
	tail -n +2 gl_chr${i}_tv_maf05_mis50.mafs >> $WDIR/gl_tv_maf05_mis50.mafs
done

gzip $WDIR/gl_tv_maf05_mis50.beagle &
gzip $WDIR/gl_tv_maf05_mis50.mafs &
wait





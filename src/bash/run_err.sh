#!/bin/bash


WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen
REF=/home/projects/ku-cbd/data/zanzibar/ref/felcat9.fasta
BAM_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_bam

CHR_LIST=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_chr_felcat9_noMT_angsd

# load angsd before use
module purge
module load tools
module load htslib/1.8 angsd/0.931


# gen anc
cd $WDIR/0-angsd_fa
BAM=/home/projects/ku-cbd/data/zanzibar/mapping/Amurleopard_PPO1.Felis_catus_9.realigned.bam
ID=Amurleopard_PPO1

# angsd -nThreads 4 -i $BAM -dofasta 2 -doCounts 1 -minQ 30 -minmapq 20 \
#                   -setminDepthInd 3 -remove_bads 1 -uniqueOnly 1 \
#                   -ref $REF -out $ID.minD3.q30.f2.cons
# 
# 
# 
ANC=$WDIR/0-angsd_fa/$ID.minD3.q30.f2.cons.fa.gz

module load samtools/1.13
# samtools faidx $ANC

# err estimate
##### run per chr next time
cd $WDIR/1-err
for i in $(less $CHR_LIST)
do
	angsd -doAncError 1 -anc $ANC -ref $REF -out err_chr.$i -bam $BAM_LIST -C 50 -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -checkBamHeaders 0 -r $i &
done
wait



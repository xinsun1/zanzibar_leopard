#!/bin/bash

module load tools
module load gcc/10.2.0 intel/perflibs/64/2020_update2  R/4.0.3
module load gsl/2.1
module load samtools/1.13
module load htslib/1.8
module load perl/5.20.1
module load ngstools/20190624
module load fastme/2.1.5


WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/5-wg-nj
cd $WDIR

BEAGLE=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl/gl_tv_maf05_mis50.qc.rel.beagle.gz
NAME=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_name_keep_rel_order
N_IND=59
N_SITE=2878628
OUT=gl_tv_maf05_mis50.qc.rel



ngsDist --pairwise_del --geno $BEAGLE --n_threads 10 --out $OUT --n_ind $N_IND --n_sites $N_SITE --probs --labels $NAME --n_boot_rep 100 --boot_block_size 10000


fastme -i $OUT -s -D 101 -o ${OUT}.fastme.tre -T 10


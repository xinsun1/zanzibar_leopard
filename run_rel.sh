#!/bin/bash

# set WDIR
WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl

MAF_FILE=gl_tv_maf05_mis50.mafs.gz
B_FILE=gl_tv_maf05_mis50.beagle.gz


# load module
module load tools
module load htslib/1.8


# prep freq
cd $WDIR
# if -ref cut f6
zcat $MAF_FILE | cut -f6 |sed 1d > freq

/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/program/ngsRelate/ngsRelate -f freq -G $B_FILE -O rel_maf05_mis50 -n 80 -p 40 



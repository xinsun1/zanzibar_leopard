#!/bin/bash

#### load modules
module purge
module load tools
module load anaconda3/2020.07

#### run
WDIR=/home/projects/ku-cbd/data/zanzibar/pop_gen/3-pca

cd $WDIR

IN_BEAGLE=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl/gl_tv_maf05_mis50.qc.rel.beagle.gz
OUT_NAME=pca_gl_tv_maf05_mis50.qc.rel

python /home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/15-gl/cinnamon/r1/pcangsd/pcangsd.py -e 5 -beagle $IN_BEAGLE -o $OUT_NAME -threads 40



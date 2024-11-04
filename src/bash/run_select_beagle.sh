#!/bin/bash

LIST_NAME_ALL=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_name
LIST_KEEP=/home/projects/ku-cbd/data/zanzibar/pop_gen/list_name_keep_rel
BEAGLE=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl/gl_tv_maf05_mis50.beagle
OUT=/home/projects/ku-cbd/data/zanzibar/pop_gen/2-gl/gl_tv_maf05_mis50.qc.rel.beagle

awk 'NR==FNR {a[$1]=1; next}{$2=0;if(a[$1]==1){$2=1};print FNR*3+1,$2; print FNR*3+2,$2; print FNR*3+3,$2}END{print 1,1; print 2,1;print 3,1}' $LIST_KEEP $LIST_NAME_ALL | awk 'NR==FNR {a[$1]=$2;next}{for(i=1;i<=NF;i++){if(a[i]==1){printf $i"\t"}}; print ""}' - $BEAGLE > $OUT

gzip -c $OUT > ${OUT}.gz
# rm $OUT

# awk 'NR==FNR {a[$1]=1;next}{if(a[$1]==1)print $1}' list_name_keep_rel list_name_qc > list_name_keep_rel_order



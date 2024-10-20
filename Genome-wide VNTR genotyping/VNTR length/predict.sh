#!/bin/bash
### danbing-tk pipeline(single sample)
### edit by sjzhang  -- 20230207
###



echo Pipeline danbingtk predict Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)

OUTDIR=$1
KMERS=$2
LSB=$3
nonTR_COV=$4



python3 /home2/sjzhang/git/danbing-tk/script/kmc2length.py --outdir $1 --ksize 21 --kmers $2 --trbed /home2/sjzhang/git/danbing-tk/resources/TR_loci.bed --LSB $3 --cov $4 --covbed /home2/sjzhang/git/danbing-tk/resources/ctrl.bed


end_time=$(date +%s)
echo "#######-- danbingtk predict End Time:`date +%Y/%m/%d--%H:%M`  --#######"



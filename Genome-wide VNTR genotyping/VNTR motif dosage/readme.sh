#!/bin/bash

# VNTR motif dosage genotyping
# select 38685 VNTR loci
less -S /home2/sjzhang/git/danbing-tk/resources/vntr.statistics.tsv|awk '{if($9=="True")print $0}'|grep -v chrX|sort -k2,2V -k3,3n -k4,4n  > 38685.stat.tsv 

#get motifs of 38685 VNTR loci
python match_38685motif.py

#get kmer counts of NyuWa cohort，calculate motif dosage
less -S 4129_sample.lst|cut -f1|grep -v Sample|split -l 40 /dev/stdin Nyuwa.sam
qsub getMotifDos_NyuBatch*.pbs


#combine motif dosage files
qsub comb_motifBatch.pbs

#get kmer counts of 1KGP and HGDP cohort，calculate motif dosage
less -S new_HGDP_1KGP_sampleinfo.txt|cut -f1|grep -v Sample|split -l 40 /dev/stdin world.sam
qsub getMotifDos_worldBatch*.pbs

#combine motif dosage files
qsub comb_motifBatch.pbs


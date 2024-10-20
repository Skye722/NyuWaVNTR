#!/bin/bash

# VNTR length genotyping
# a demo PBS file of each sample
demo.pbs

# running PBS file of each sample, get kmer file
qsub demo.pbs

# get loci specific bias(LSB) 
qsub getCovByLocus.pbs 
python /home2/sjzhang/git/danbing-tk/LSB/get_LSB.py 

# combine kmer files of all samples
ls *.tr.kmers > combine_kmers.lst

# combine LSB file of all samples. Each line represents a VNTR locus and each column represents a sample

# combine nonTR coverage file of all samples. Each line represents a sample and each column represents the coverage of each VNTR locus.

# predict VNTR length
qsub predict.pbs
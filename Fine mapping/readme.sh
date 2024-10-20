#!/bin/bash


# fine-mapping of eVNTRs by susieR
for i in `seq 1 22`;do Rscript causalVNTR.identification.R chr$i;done

# fine-mapping of eMotifs by susieR
for i in `seq 1 22`;do Rscript causalMotif.identification.R chr$i;done
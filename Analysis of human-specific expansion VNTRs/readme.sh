#!/bin/bash

# The file of VNTR loci where 38685 ref loci and hse loci in Course M M et al., 2021 overlap
hse_olp_col.loci

# calculation of fold change
for pop in AMR CSA EAS EUR MID OCE;do python get_VNTR_TL_FC.py -a ../AFR.lst -b ../$pop.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o AFR$pop.fc;done
for pop in  CSA EAS EUR MID OCE;do python get_VNTR_TL_FC.py -a ../AMR.lst -b ../$pop.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o AMR$pop.fc;done
for pop in  CSA EAS MID OCE;do python get_VNTR_TL_FC.py -a ../EUR.lst -b ../$pop.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o EUR$pop.fc;done
for pop in  CSA MID OCE;do python get_VNTR_TL_FC.py -a ../EAS.lst -b ../$pop.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o EAS$pop.fc;done
for pop in  MID OCE;do python get_VNTR_TL_FC.py -a ../CSA.lst -b ../$pop.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o CSA$pop.fc;done
python get_VNTR_TL_FC.py -a ../MID.lst -b ../OCE.lst -v 8222_estimated_good_TR_len_bc_final.tsv -o MIDOCE.fc


# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 16:27:48 2023

@author: skye
"""

with open("/home2/sjzhang/work/VNTRv2/4-eVNTR/Geuvadis_441.VNTR_gene_pair.csv", 'r') as f, open("Geuvadis_441_MotDos_rmOutliter_rmMotif.tsv", 'r') as fh, open("Geuvadis_441.motif_gene_pair.csv", 'w') as fout:
    fout.write("site_motif,gene_id,distance")
    _ = f.readline()
    _ = fh.readline()

    fd = {}
    for line in f:
        line = line.strip().split(",")
        fd.setdefault(line[0], []).append(line[1:])

    for line in fh:
        site_mot = line.strip().split("\t")[0]
        site = "_".join(site_mot.split("_")[0:3])
        motif = site_mot.split("_")[-1]
        #print(site,motif)

        if site in fd:
            gene_dis_s = fd[site]
            for gene, dis in gene_dis_s:
                fout.write(f"\n{site}_{motif},{gene},{dis}")

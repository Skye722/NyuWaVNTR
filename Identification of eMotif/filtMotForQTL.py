#codings:utf-8

import pandas as pd
import numpy as np
from scipy.stats import zscore


with open("Geuvadis_441_MotDos.tsv", 'r') as f,open("Geuvadis_441_MotDos_rmOutliter.csv", "w") as fout:
    head = f.readline()
    fout.write("\t".join([str(x) for x in head.strip().split("\t")])+"\n")
    for line in f:
        line = line.strip().split('\t')
        dos = [float(x) for x in line[2:]]
        row_mean = np.mean(dos)
        row_std = np.std(dos)
        dos_filt = [np.nan if abs(x-row_mean) > 2 * row_std else x for x in dos]
        fout.write(line[0] + "\t" + line[1] + "\t" + "\t".join([str(x) for x in dos_filt])+"\n")



# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 17:49:02 2022

@author: skye
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib
import matplotlib.pyplot as plt


def NonparametricOutlierRemoval(x, k=2):
    """
    input will be flattened. Return: x_non_outlier, mask_non_outlier
    """
    q0, q1 = np.quantile(x, [0.25,0.75])
    kiqr = (q1 - q0) * k
    mask = np.logical_and(x >= q0 - kiqr, x <= q1 + kiqr)
    return x[mask], mask

def getLSB(covm, csize):
    """
    Return: LSB, mask_non_outlier_loci
    """
    ncov = covm / (covm @ csize / np.sum(csize))[:,None]
    ms = np.mean(ncov, axis=0)
    vs = np.var(ncov, axis=0)
    
    mask = np.logical_and(NonparametricOutlierRemoval(ms)[1], NonparametricOutlierRemoval(vs)[1])
    return ncov, mask


covfile="/home2/sjzhang/git/danbing-tk/sample.nonTRandTR.cov"
ctrlbed="/home2/sjzhang/git/danbing-tk/resources/nonTR_TR_loci.bed"
covmat = np.loadtxt(covfile, dtype=object,delimiter='\t')[:,1:]
ctrlbed = np.loadtxt(ctrlbed, dtype=object)
csize = ctrlbed[:,2].astype(int) - ctrlbed[:,1].astype(int) #locus end - locus start

    
# covm = covmat[:,1:].astype(float) 

LSB, Lm = getLSB(covmat[:,1:].astype(float) , csize)

sample = np.insert(covmat[:,:1], 0,"Loci")
#print(len(sample))
pos = np.array(["_".join(i) for i in ctrlbed])
LSBOut = np.concatenate(([pos], LSB ),axis = 0)
LSBOut = np.concatenate((sample.reshape((len(sample),1)), LSBOut ),axis = 1)

np.savetxt("/home2/sjzhang/git/danbing-tk/sample.LSB.nonTRandTR.txt",LSBOut.T,fmt='%s',delimiter ='\t')




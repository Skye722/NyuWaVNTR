#codings:utf-8

import pandas as pd
import numpy as np
import sys
from scipy.stats import zscore

def rmOutliers(filename):
    df=pd.read_csv(filename, header=0, index_col=0)
    row_means = df.mean(axis=1)
    row_stds = df.std(axis=1)
    for row in df.index:
        for col in df.columns:
            if abs(df.loc[row, col] - row_means[row]) > 3 * row_stds[row]:
                df.loc[row, col] = np.nan
                #print(row, col)
    
    return df

def VNTRGenePair(win100filename):
    str_gene_pair = {}
    with open(win100filename) as fin:
        for line in fin:
            line = line.strip().split('\t')
            site = line[6] + "_" + line[7] + "_" + line[8]
            # distance
            d = str(abs((int(line[2])-int(line[1])) - (int(line[8])-int(line[7]))))
            gene = line[3]
            if site not in str_gene_pair:
                str_gene_pair[site] = [(gene, d)]
            else:
                str_gene_pair[site].append((gene, d))
    return str_gene_pair



def main():
    df = rmOutliers(sys.argv[1])
    df.to_csv(sys.argv[2],na_rep='NA')
    str_gene_pair = VNTRGenePair(sys.argv[3])
    with open("GD462.VNTR_gene_pair.csv",'wt') as fout:
        fout.write('site,gene_id,distance\n')
        num_str = len(str_gene_pair)
        genes = {}
        num_pairs = 0
        for i in str_gene_pair:
            num_pairs += len(str_gene_pair[i])
            for g,d in str_gene_pair[i]:
                genes[g] = 0
                fout.write(i + ',' + g + ',' + d + '\n')
    num_genes = len(genes)
    print("%s VNTRs." %(num_str))
    print("%s genes." %(num_genes))
    print("%s VNTR-gene pairs." %(num_pairs))


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python Geuvadis_441_VNTRlen.csv output Geuvadis_445.genes_window100k_vntr.txt")
    else:
        main()

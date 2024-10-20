#!/usr/bin/env python

'''
purpose: get length LFC of each STR locus for two populations

Usage:
pip install vcfpy
pip install scipy
python3 get_STR_FC.py -a pop1.lst -b pop2.lst -v input.vcf -o pop1_vs_pop2.FC.txt
python get_STR_FC.py -a lst.CHN.old -b lst.CHS.old -v /Parastor300s_G30S/shiyr/v7/v7-dump-top100-sort.vcf -o test.txt

output:
site Rst
chr17:32123599:tttc:4 0.01
'''

import sys
import os
import argparse
import logging
import numpy as np
from scipy.stats import ranksums
from statistics import mean
from math import log2


# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """
    description = "%(prog)s -- Get length LFC of each STR locus between two populations."
    epilog = "For command line options, type: %(prog)s -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description, epilog=epilog) #, usage = usage )
    argparser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    argparser.add_argument('-a', required=True, type=str, help='Samples of population 1 (one sample per line, no header).', metavar='', dest="pop1")
    argparser.add_argument('-b', required=True, type=str, help='Samples of population 2 (one sample per line, no header).', metavar='', dest="pop2")
    argparser.add_argument('-v', '--vcf', required=True, type=str, help='Path to STR VCF file.', metavar='', dest="vcf")
    argparser.add_argument('-o', '--output', required=True, type=str, help='Path to output file.', metavar='', dest="output")

    return argparser


def get_lfc(vntrLenfile, sam_pop1, sam_pop2):
    '''get log2 fold change, and p-value
    '''
    # read pop list
    with open(sam_pop1, 'r') as pop1:
        pop1s = [line.strip() for line in pop1]
    with open(sam_pop2, 'r') as pop2:
        pop2s = [line.strip() for line in pop2]

    # vntr length
    with open(vntrLenfile) as f:
        samples = f.readline().strip().split(",")
        pop1s_idx = [samples.index(item) for item in pop1s]
        pop2s_idx = [samples.index(item) for item in pop2s]
        # print(pop1s_idx[0:10])
        # print(pop2s_idx[0:10])

        for line in f:
            line = line.strip().split(",")
            locus = line[0]
            # get length of 2 alleles
            len_pop1 = [float(line[idx]) if line[idx] else np.nan for idx in pop1s_idx]
            len_pop2 = [float(line[idx]) if line[idx] else np.nan for idx in pop2s_idx]      
            # print(len_pop1[0:10])
            # print(len_pop2[0:10])

            # wilcox.test
            res = ranksums(len_pop1, len_pop2)
            p = res.pvalue
            # fold change
            m1 = log2(mean(len_pop1) + 0.0)
            m2 = log2(mean(len_pop2) + 0.0)
            LFC = m2 - m1

            yield locus, str(p), str(LFC), str(len(len_pop1)), str(len(len_pop2))


def main():
    """The main function
    """
    argparser = prepare_argparser()
    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    pop1 = args.pop1
    pop2 = args.pop2
    vntrLenfile = args.vcf

    # sample info
    sam_pop1 = {}
    try:
        logger.info("Read %s" %(pop1))
        fin = open(pop1, 'rt')
    except:
        logger.error("Read input file %s error." %(pop1))
        sys.exit(1)
    else:
        for line in fin:
            line = line.strip()
            sam_pop1[line] = 0
        fin.close()
        logger.info("There are %s samples in %s." %(len(sam_pop1), pop1))
    sam_pop2 = {}
    try:
        logger.info("Read %s" %(pop2))
        fin = open(pop2, 'rt')
    except:
        logger.error("Read input file %s error." %(pop2))
        sys.exit(1)
    else:
        for line in fin:
            line = line.strip()
            sam_pop2[line] = 0
        fin.close()
        logger.info("There are %s samples in %s." %(len(sam_pop2), pop2))
    # check olp of sam_pop1 and sam_pop2
    tmp_len = len(set(sam_pop1) & set(sam_pop2))
    if tmp_len > 0:
        logger.error("Some %s samples are in both %s and %s." %(tmp_len, pop1, pop2))
        sys.exit(1)


    # output
    output = args.output
    output_path = os.path.abspath(os.path.dirname(output))
    if not os.path.exists(output_path): os.mkdir(output_path)
    fout = open(output, 'wt')
    fout.write('site\talleleNum_pop1\talleleNum_po2\tLFC\tranksums.p\n')


    # parse
    # get site-level-stat
    for (locus, p, LFC, len_pop1, len_pop2) in get_lfc(vntrLenfile, pop1, pop2):
        to_w = [locus, len_pop1, len_pop2, LFC, p]
        # write
        fout.write('\t'.join(to_w) + '\n')
    logger.info("Complete!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)



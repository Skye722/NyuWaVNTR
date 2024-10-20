# coding:utf-8

import sys
import os
import gzip
import shutil

def read_f(fileName):
    if fileName.endswith('.gz'):
        f = gzip.open(fileName, "rt")
        for line in f:
            line = line.rstrip("\n").split("\t")
            line = line[:2] + [item.split(":")[-1] for item in line[3:]]
            yield line
        f.close()
    else:
        f = open(fileName, "rt")
        for line in f:
            line = line.rstrip("\n").split("\t")
            yield line
        f.close()
        


def combine_f(d1, d2, tmpFileName):
    with open(tmpFileName,"wt") as fout:
        for i1, i2 in zip(d1, d2):
            [i1.append(i2i) for i2i in i2[2:]]
            fout.write("\t".join(i1))
            fout.write("\n")


def main():
    combined = None
    outFileNames = [outFileName+".tmp1", outFileName+".tmp2"]
    idx = 0
    for fileName in fileNames:
        print(fileName)
        file_iter = read_f(fileName)
        if not combined:
            combined = file_iter
        else:
            tmpFileName = outFileNames[idx%2]
            idx+=1
            combine_f(combined, file_iter, tmpFileName)
            combined = read_f(tmpFileName)
    with gzip.open(outFileName,"wt") as fout:
        for line in combined:
            fout.write("\t".join(line))
            fout.write("\n")
    os.remove(outFileNames[0])
    if os.path.isfile(outFileNames[1]):
        os.remove(outFileNames[1])
    



if __name__ == '__main__':
    # python test.py f* out
    if len(sys.argv) < 3:
        pass
    else:
        outFileName = sys.argv[-1]
        fileNames = sys.argv[1:-1]
        main()

#coding:utf-8

# extract corresponding motifs from 38685 VNTR loci
with open('38685.stat_sort.tsv', 'r') as f:
    f.readline()
    nums_file1 = [int(line.strip().split()[0]) for line in f]


with open('/home2/sjzhang/git/danbing-tk/resources/danbing_1.3.1/metadata/TR_locus.4456881_motif.tsv', 'r') as f:
    with open('38685TR_locus_motif.tsv', 'w') as out:
        f.readline()
        for line in f:
            num = int(line.strip("\t").split()[1])
            if num in nums_file1:
                out.write(line)

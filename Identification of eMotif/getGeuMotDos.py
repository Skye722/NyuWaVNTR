#codings:utf-8
import select
import sys
import statistics

def read100kVNTR(filename):
    with open(filename, 'r') as f:
        site = [line.strip().split("\t")[-1] for line in f]
    return site

def getGeuMotDos(filename,samFile):
    with open(samFile, 'r') as samf:
        sam = [line.strip() for line in samf]
        sam.remove('NA20530')
        sam.remove('NA20506')
        sam.remove('HG00154')
        sam.remove('HG00361')
    with open(filename, 'r') as fHandle:
        dos=dict()
        head=fHandle.readline().strip().split("\t")
        idx = [head.index(sample) for sample in sam]
        for line in fHandle:
            line = line.strip().split("\t")
            if not line[1] in dos:
                dos[line[1]]=[[line[2]]+[line[x] for x in idx]]
            else:
                dos[line[1]].append([line[2]]+[line[x] for x in idx])
    return sam,dos

def main(): 
    site = read100kVNTR(sys.argv[1])
    sam,dos = getGeuMotDos(sys.argv[2],sys.argv[3])
    select_dos = {key:dos[key] for key in site if key in dos.keys()}
    with open(sys.argv[4], 'w') as outfile:
        outfile.write("site\tmotif\t"+"\t".join(sam)+"\n")
        for key,value in select_dos.items():
            for v in value:
                outfile.write(key+"\t"+"\t".join(v)+"\n")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python getGeuMotDos.py tmp /home2/sjzhang/work/VNTRv2/2-QC/callRate/World_and_Nyuwa_motifDos_final.tsv /home2/niuyw/project/Geuvadis/geno_pca/Geuvadis_olp_1KG2504.lst outfilename")
    else:
        main()    
            



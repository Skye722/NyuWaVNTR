# coding:utf-8

import os
import sys
import numpy
import gzip


def read_motif(filename):
    with open(filename, 'r') as f:
        header = f.readline()
        yield header.strip().split()
        for line in f:
            line = line.strip().split()
            yield line


def read_kmer(filename):
    with open(filename, 'r') as f:
        kmer_seqs = []
        mark = ""
        for line in f:
            if not line.startswith(">"):
                kmer_seqs.append(line.strip())
            else:
                if mark != "" and len(kmer_seqs) != 0:
                    yield mark, kmer_seqs
                kmer_seqs = []
                mark = line[1:].strip()
        yield mark, kmer_seqs


def read_sample_filenames(filename, dirPath="./"):
    with open(filename, 'r') as f:
        return [os.path.join(dirPath, "{}.tr.kmers".format(l.strip())) for l in f]


def main(args):
    if len(args) < 5:
        return

    motif_file, kmer_file, sample_file, sample_dir, out_file = args

    sample_filenames = read_sample_filenames(sample_file, sample_dir)

    out_f = gzip.open(out_file + ".gz", 'wt')

    kmer_iter = read_kmer(kmer_file)
    sample_iters = [read_kmer(sample_filename)
                    for sample_filename in sample_filenames]

    motif_iter = read_motif(motif_file)
    motif_header = next(motif_iter)

    out_f.write("\t".join(motif_header + ["kmers"] + [sample_filename.split(
        "/")[-1].replace(".tr.kmers", "") for sample_filename in sample_filenames])+"\n")

    motif_last_locus = "-1"
    kmer_info = None
    sample_infos = []
    cur_kmer_block_size = 1
    for motif_info in motif_iter:
        _, motif_locus, motif_seq = motif_info

        # 当读取到新的locus时
        if motif_locus != motif_last_locus:
            # 读取 一块 sample
            sample_infos = []
            for sample_filename_idx, sample_filename in enumerate(sample_filenames):
                sample_locus = -1
                # 判断一下与当前的的locus是不是一致的
                while int(sample_locus) < int(motif_locus):
                    sample_locus, sample_kmers = next(
                        sample_iters[sample_filename_idx])

                if sample_locus != motif_locus:
                    print("sample_locus != motif_locus")
                    raise Exception("sample_locus != motif_locus")

                sample_infos.append(
                    [sample_filename_idx, sample_kmers, sample_locus])
            # 读取 一块 kmer
            kmer_locus = -1
            while int(kmer_locus) < int(motif_locus):
                kmer_info = next(kmer_iter)
                kmer_locus, kmer_seqs = kmer_info
                cur_kmer_block_size += len(kmer_seqs)+1

            if kmer_locus != motif_locus:
                print("kmer_locus != motif_locus")
                raise Exception("kmer_locus != motif_locus")
            motif_last_locus = motif_locus
            cur_kmer_block_size -= (len(kmer_seqs)+1)

        if kmer_info and sample_infos:
            # 查找
            _, kmer_seqs = kmer_info
            _, tmp_motif_locus, tmp_motif_seq = motif_info

            # 获取哪些kmer在当前的motif
            hit_kmers = [[i, seq] for i, seq in enumerate(
                kmer_seqs) if seq in tmp_motif_seq]

            # 用于存储不同样本的kmer的那个数字
            hit_kmer_counts = [[] for sample_filename in sample_filenames]

            # 遍历所有样本的kmer数据
            for sample_filename_idx, sample_kmers, sample_locus in sample_infos:
                # 将 样本kmer后面的数字添加到对应的hit_kmer_counts列中
                [hit_kmer_counts[sample_filename_idx].append(
                    sample_kmers[hit_kmer[0]].strip().split()[1]) for hit_kmer in hit_kmers]
                
            # 将结果拼成一行数据，并输出
            out_info = [tmp_motif_locus, tmp_motif_seq] + [",".join([str("{}:{}").format(kmerfileidx+cur_kmer_block_size, kmerseq) for kmerfileidx, kmerseq in hit_kmers])] + [
                str("{}:{}").format(",".join(item), numpy.mean(list(map(float, item)))) for item in hit_kmer_counts]
            out_f.write("\t".join(out_info)+"\n")

    out_f.close()


if __name__ == "__main__":
    # motif_file, kmer_file, sample_file, sample_dir, out_file
    main(sys.argv[1:])

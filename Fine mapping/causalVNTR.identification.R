#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/eSTR.identification.R > analysis/eSTR.identification.log &

DOCNAME <- "fineMap"


# Tidyverse
library(tidyverse)
library(broom)

# plot
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
library(ggrastr)
library(EnvStats)
theme_set(theme_cowplot(
  font_size = 12, # axis.title
  rel_small = 10 / 12, # axis.text
  rel_tiny = 9 / 12,
  rel_large = 12 / 12
))

# heatmap
library(pheatmap)
library(corrplot)

# color
library(ggsci)
library(RColorBrewer)

#SusieR
library(susieR)

#output
library(data.table)

# Load exp
rna_exp <- readRDS(here::here("output/eVNTR","Geuvadis.exp_adjusted.rds"))
rna_exp[1:3, 1:3]
colnames(rna_exp) <- sub("\\..*", "", colnames(rna_exp))
rows_to_remove <- c("NA20530", "NA20506", "HG00154","HG00361")
rna_exp <- rna_exp[!(rownames(rna_exp) %in% rows_to_remove), ]
rna_exp[1:3, 1:3]
# write.csv(rna_exp, file = here::here("output/eVNTR", "Geuvadis.exp_adjusted.csv"),quote = FALSE)
# Load chr
argv <- commandArgs(T)

# Load VNTR-SNP geno
file <- paste(argv[1],"-tss1MBSNP.txt.gz",sep = "")
vntr_snp_geno <- read.csv(here::here("data/fineMap", file) ,header = TRUE,sep = ",")
vntr_snp_geno[1:7, 1:7]

# Load 543 VNTR-gene pairs
df.vntr_gene_pair <- read_csv(here::here("output/eVNTR",  "df.eVNTR.FDR0.05.csv"))
head(df.vntr_gene_pair)

# chr gene-vntr lst
chr_df <- subset(df.vntr_gene_pair, grepl(paste("^",argv[1],"_",sep = ""), site))
gene_vntr_lst <- split(chr_df$site, chr_df$gene)

## test
# ofn <- paste( "/home/zhangsj/vntr/output/fineMap/test", paste(".pip.", "L5", ".txt", sep=""), sep="")
# ff <- FALSE
# geno <- subset(vntr_snp_geno, geneId== "ENSG00000134463")
# geno_matrix <- t(geno[,6:ncol(geno)])
# colnames(geno_matrix) <- paste(geno$chrom,geno$snpPos,sep = "_")
# geno_matrix <- as.matrix(geno_matrix)
# gene <- "ENSG00000134463"
# g_exp <- rna_exp[,gene]
# names(g_exp) <- rownames(rna_exp)
# res <- susie(geno_matrix,g_exp,L=5,na.rm=FALSE)
# pip <- matrix(NA, 1, dim(geno_matrix)[2],dimnames = list(gene, colnames(geno_matrix)))
# pip <- res$pip
# cat(".")
# if (is.null(res$sets$cs$L1)) {
#   fwrite(list(NA), ofn, append=ff, sep="\t", col.names=FALSE, na="nan")
# } else {
#   inds <- c(1:dim(geno_matrix)[2])
#   L1 <- t(inds[res$sets$cs$L1])
#   fwrite(as.data.table(gene), ofn, append=ff, sep="\t", col.names=FALSE, na="nan")
# }
# ff <- TRUE
# output_data <- data.frame(colnames = names(pip), values = unname(pip))
# fwrite(output_data, ofn, append=ff, sep="\t", col.names=FALSE, na="nan",quote = FALSE)

# genes
genes <- names(gene_vntr_lst)
genes <- sub("\\..*", "", genes)

#func
replace_with_column_means <- function(col_values) {
  mean_value <- mean(col_values, na.rm = TRUE)
  col_values[is.na(col_values)] <- mean_value
  return(col_values)
}

# fit
ofn <- paste( "/home/zhangsj/vntr/output/fineMap/",argv[1], paste(".pip.", "L5", ".txt", sep=""), sep="")
fit_res <- lapply(1:length(genes), function(i) {
  print(i)
  g <- genes[i]
    # prep
  if (g %in% geneId) {
    geno <- subset(vntr_snp_geno, geneId==g)
    print(geno[1:7,1:7])
    geno_matrix <- t(geno[,6:ncol(geno)])
    colnames(geno_matrix) <- paste(geno$chrom,geno$snpPos,sep = "_")
    geno_matrix <- as.matrix(geno_matrix)
    print(geno_matrix[1:7,1:7])
    geno_matrix <- as.data.frame(apply(geno_matrix, 2, replace_with_column_means))
    geno_matrix <- as.matrix(geno_matrix)
    print(geno_matrix[1:7,1:7])
    g_exp <- rna_exp[, g]
    names(g_exp) <- rownames(rna_exp)
    # susieR
    res <- susie(geno_matrix,g_exp,L=5,na.rm=FALSE)
    # tidy
    pip <- matrix(NA, 1, dim(geno_matrix)[2])
    pip <- res$pip
    cat(".")
    if (is.null(res$sets$cs$L1)) {
      fwrite(list(NA), ofn, append=TRUE, sep="\t", col.names=FALSE, na="nan")
    } else {
      inds <- c(1:dim(geno_matrix)[2])
      L1 <- t(inds[res$sets$cs$L1])
      fwrite(as.data.table(g), ofn, append=TRUE, sep="\t", col.names=FALSE, na="nan")
    }
    output_data <- data.frame(colnames = names(pip), values = unname(pip))
    fwrite(output_data, ofn, append=TRUE, sep="\t", col.names=FALSE, na="nan",quote = FALSE)
  }else{
    print(g)
  }
})

# Save
# saveRDS(fit_res, file = here::here("output", DOCNAME,"eVNTR.exp_vntr.glm_fit-GD462.rds"))


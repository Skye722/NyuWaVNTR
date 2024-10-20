#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/eSTR.identification.R > analysis/eSTR.identification.log &

DOCNAME <- "eVNTR"


# Tidyverse
library(tidyverse)
library(broom)
library(readr)
library(dplyr)
# plot
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
# library(ggrastr)
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

#load FPKM
log_FPKM = read_csv(here::here("data/eVNTR", "Geuvadis.x_tmm.445.log2FPKM.csv"))

#Geuvadis 462
#cols_to_remove <- c("HG00104","HG00124","HG00134","HG00135","HG00152","HG00154","HG00156","HG00247","HG00249","HG00312","HG00359","HG00361","HG00377","NA07346","NA11993","NA18487","NA19150","NA20506","NA20530","NA20537","NA20816")

#Geuvadis 445
cols_to_remove <- c('NA20530','NA20506','HG00154','HG00361')
log_FPKM <- select(log_FPKM,-cols_to_remove)
dim(log_FPKM)

# load exp
rna_exp <- readRDS(here::here("output", DOCNAME,"Geuvadis.exp_adjusted.rds"))
rna_exp[1:3, 1:3]
colnames(rna_exp) <- sub("\\..*", "", colnames(rna_exp))
rna_exp[1:3, 1:3]
rows_to_remove <- c("NA20530", "NA20506", "HG00154","HG00361")
rna_exp <- rna_exp[!(rownames(rna_exp) %in% rows_to_remove), ]

# Load motif dosage
vntr_geno <- read_tsv(here::here("data/eVNTR",  "Geuvadis_441_MotDos_rmOutliter.csv"))
vntr_geno <- transform(vntr_geno, index_col = paste(vntr_geno$site, vntr_geno$motif, sep = "_"))
rownames(vntr_geno) <- vntr_geno$index_col
vntr_geno <- vntr_geno[, -(0:2)]
vntr_geno <- vntr_geno[, -ncol(vntr_geno)]
vntr_geno[0:3, 0:3]

#Neat motifs with non-0.0/nan > 50
row_counts <- rowSums(vntr_geno != 0 & !is.na(vntr_geno))
vntr_geno_filt <- vntr_geno[row_counts >50, ]
dim(vntr_geno_filt)
#write.table(vntr_geno_filt, file = here::here("output", DOCNAME, "Geuvadis_441_MotDos_rmOutliter_rmMotif.tsv"),row.names = TRUE,quote = FALSE,sep="\t")
#saveRDS(vntr_geno_filt, file = here::here("output", DOCNAME,"Geuvadis_441_MotDos_rmOutliter_rmMotif.rds"))

# VNTR-gene pairs
df.vntr_gene_pair <- read_csv(here::here("data/eVNTR",  "Geuvadis_441.motif_gene_pair.csv"))
head(df.vntr_gene_pair)

# gene-vntr lst
gene_vntr_lst <- split(df.vntr_gene_pair$site_motif, df.vntr_gene_pair$gene_id)

# test
# s <- "chr1_806322_806403_ACTTCAAACTCGGGCAGTGCATGCCTCCCC"
# g <- "ENSG00000228327"
# s_geno <- na.omit(unlist(vntr_geno[s, ]))
# g_exp <- rna_exp[, g]
# names(g_exp) <- rownames(rna_exp)
# g_exp <- g_exp[names(s_geno)]
# dat <- data.frame(sample = names(g_exp), exp = scale(g_exp), geno = scale(s_geno))
# x <- glm(formula = exp ~ geno, data = dat)

# genes
genes <- names(gene_vntr_lst)
genes <- sub("\\..*", "", genes)

# fit
fit_res <- lapply(1:length(gene_vntr_lst), function(i) {
  g <- genes[i]
  fit_res_lst <- lapply(gene_vntr_lst[[i]], function(s) {
    # prep
    s_geno <- na.omit(unlist(vntr_geno_filt[s, ]))
    g_exp <- rna_exp[, g]
    names(g_exp) <- rownames(rna_exp)
    g_exp <- g_exp[names(s_geno)]
    dat <- data.frame(sample = names(g_exp), exp = scale(g_exp), geno = scale(s_geno))
    print(g)
    print(s)
    # glm
    fit <- glm(formula = exp ~ geno, data = dat)
    # tidy
    fit_df <- broom::tidy(fit) %>%
      filter(term == "geno")
    fit_df[1, 1] <- s
    # neat
    fit_df %>%
      mutate(gene = g) %>%
      dplyr::select(gene, site = term, estimate, std.error, statistic, p.value)
  })
  do.call(rbind, fit_res_lst) %>%
    mutate(q.value = p.adjust(p.value, method = "bonferroni")) %>%
    arrange(-abs(estimate / std.error), p.value)
})
names(fit_res) <- genes

# Save
saveRDS(fit_res, file = here::here("output", DOCNAME,"eVNTR.exp_motif.glm_fit.rds"))


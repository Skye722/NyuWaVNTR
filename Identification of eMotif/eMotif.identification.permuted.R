#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/STR/str_rproj
# nohup Rscript analysis/eSTR.identification.permuted.R > analysis/eSTR.identification.permuted.log &

DOCNAME <- "eVNTR"


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

# load exp
rna_exp <- readRDS(here::here("output", DOCNAME,"Geuvadis.exp_adjusted.rds"))
rna_exp[1:3, 1:3]
colnames(rna_exp) <- sub("\\..*", "", colnames(rna_exp))
rna_exp[1:3, 1:3]
rows_to_remove <- c("NA20530", "NA20506", "HG00154","HG00361")
rna_exp <- rna_exp[!(rownames(rna_exp) %in% rows_to_remove), ]

# permutate exp (shuffle sample labels)
set.seed(1234)
rownames(rna_exp) <- sample(rownames(rna_exp), length(rownames(rna_exp)))

# Load motif dosage
vntr_geno <- read.table(here::here("output/eVNTR",  "Geuvadis_441_MotDos_rmOutliter_rmMotif.tsv"),header = TRUE, sep = "\t", row.names = 1)
vntr_geno[1:3, 1:3]

# motif-gene pairs
df.vntr_gene_pair <- read_csv(here::here("data/eVNTR",  "Geuvadis_441.motif_gene_pair.csv"))
head(df.vntr_gene_pair)

# gene-vntr lst
gene_vntr_lst <- split(df.vntr_gene_pair$site_motif, df.vntr_gene_pair$gene_id)

# # test
# s <- "chr1_806322_806403_ACTTCAAACTCGGGCAGTGCATGCCTCCCC"
# g <- "ENSG00000228327"
# s_geno <- na.omit(unlist(vntr_geno[s, ]))
# g_exp <- rna_exp[, g]
# print(g)
# print(s)
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
    s_geno <- na.omit(unlist(vntr_geno[s, ]))
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
saveRDS(fit_res, file = here::here("output", DOCNAME, "eVNTR.exp_motif.glm_fit.permuted.rds"))

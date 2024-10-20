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
rna_exp <- readRDS(here::here("output", DOCNAME, "Geuvadis.exp_adjusted-GD462.rds"))
rna_exp[1:3, 1:3]
colnames(rna_exp) <- sub("\\..*", "", colnames(rna_exp))
rna_exp[1:3, 1:3]

# permutate exp (shuffle sample labels)
set.seed(1234)
rownames(rna_exp) <- sample(rownames(rna_exp), length(rownames(rna_exp)))

# Load VNTR geno
vntr_geno <- read_csv(here::here("data/eVNTR",  "GD462_VNTRlen_rmOutliter.csv")) %>%
  column_to_rownames("site")
vntr_geno[1:3, 1:3]

# VNTR-gene pairs
df.vntr_gene_pair <- read_csv(here::here("data/eVNTR",  "GD462.VNTR_gene_pair.csv"))
head(df.vntr_gene_pair)

# gene-vntr lst
gene_vntr_lst <- split(df.vntr_gene_pair$site, df.vntr_gene_pair$gene_id)

# # test
# s <- "chr10_101985147_101985326"
# g <- "ENSG00000120029.13"
# s_geno <- na.omit(unlist(vntr_geno[s, ]))
# g_exp <- rna_exp[, g]
# names(g_exp) <- rownames(rna_exp)
# g_exp <- g_exp[names(s_geno)]
# dat <- data.frame(sample = names(g_exp), exp = scale(g_exp), geno = scale(s_geno))
# x <- glm(formula = exp ~ geno, data = dat)

# genes
genes <- names(gene_vntr_lst)

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
saveRDS(fit_res, file = here::here("output", DOCNAME, "eVNTR.exp_vntr.glm_fit-GD462.permuted.rds"))

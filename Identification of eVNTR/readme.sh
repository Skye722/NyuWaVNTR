#!/bin/bash

# 1. Calculation of PCs by PCAForQTL
#Load log2 FPKM
```{r load-exp, cache=FALSE, warning=FALSE, message=FALSE}
log_FPKM = read.csv(here::here("data/eVNTR", "Geuvadis.x_tmm.445.log2FPKM.csv"))
log_FPKM[1:3, 1:3]
```

#Run PCA
```{r fig.width=6, fig.height=2.5}
bb<-t(log_FPKM[,-1]) #row-samplename(445),col-genename(15627)
dim(bb)
prcompResult<-prcomp(bb,center=TRUE,scale.=TRUE) #calculate PC
PCs<-prcompResult$x
dim(PCs) #[1] 445 445
importanceTable<-summary(prcompResult)$importance
PVEs<-importanceTable[2,]
sum(PVEs)
```

#Choose K pcs
```{r fig.width=3.5, fig.height=3.5}
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow) #[1] 25
K_elbow<-resultRunElbow
p2<-PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow"),values=c(K_elbow),titleText="Colon - Transverse")
plot(p2,xlab="PC index",ylab="PVE")
PCsTop<-PCs[,1:25]
```

#Add known covariates for PCs calculation
```{r}
knownCovariates<-read.csv(here::here("data/eVNTR", "sex2numeric.combine.cov.nocomma.hasgenename.noheader.csv"),header=F,row.names=1)
colnames(knownCovariates)[1:11]<-paste0("SexPopPC",1:11)
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,PCsTop,unadjustedR2_cutoff=0.9)
PCsTop<-scale(PCsTop)
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)
```

#Plot hidden covariates correlation
```{r}
corrplot(cor(covariatesToUse), type = "lower",
         method = "square", tl.cex = 0.8)
```

# 2. Adjust expression matrix with selected PCs
#Adjust expression matrix
```{r}
# exp
dat <- as.data.frame(t(log_FPKM))
rownames(dat) <- dat$gene_id
colnames(dat) <- dat[1, ]
dat <- dat[-1, ]

# sex+Pop10+pc25
dat <- cbind(covariatesToUse, dat)

# genes
genes <- log_FPKM$gene_id

# formula
univ_formulas <- sapply(genes, function(y) {
  as.formula(paste0(y, " ~ SexPopPC1+SexPopPC2+SexPopPC3+SexPopPC4+SexPopPC5+SexPopPC6+SexPopPC7+SexPopPC8+SexPopPC9+SexPopPC10+SexPopPC11+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25"))
})

# adjust exp
dat_adj <- sapply(univ_formulas, function(f) {
  # print(x)
  lm.res <- lm(f, data = dat)
  lm.res$residuals
})
```

#Save
```{r}
saveRDS(dat_adj, file = here::here("output", "Geuvadis.exp_adjusted.rds"))
```

# 3. eVNTR identification
# get VNTRs within 100kb of genes expressed in LCL data set
awk 'BEGIN{FS=OFS="\t"}NR==FNR{A[$1]=0}NR>FNR{if($4 in A){print $0}}' /home2/niuyw/project/STR/eSTR/Geuvadis_445.expressed_genes.txt /home2/niuyw/RefData/Homo_sapiens/GENCODE_v34/gencode.v34.genes.bed6 | sort -k1,1 -k2,2n > Geuvadis_445.expressed_genes.bed6
sed '1d' ../2-QC/callRate/8222_estimated_good_TR_len_bc_final.tsv |cut -d',' -f1|awk -F'_' '{OFS="\t";print $1,$2,$3}'|sort -k1,1 -k2,2n|bedtools window -a /home2/niuyw/project/STR/eSTR/Geuvadis_445.expressed_genes.bed6 -b stdin -w 100000 > Geuvadis_445.genes_window100k_vntr.txt

#filter VNTRs and normalize VNTR length
less -S Geuvadis_445.genes_window100k_vntr.txt |awk -F'\t' '{OFS="_";print $7,$8,$9}'|sort|uniq > gene_VNTR.loci
>>> import pandas as pd
>>> df=pd.read_csv("../callRate/8222_estimated_good_TR_len_bc_final.tsv",index_col=0)
>>> with open("/home2/niuyw/project/Geuvadis/geno_pca/Geuvadis_olp_1KG2504.lst")as f:
...     sams=[i.strip() for i in f]
>>> sams.remove('NA20530')
>>> sams.remove('NA20506')
>>> sams.remove('HG00154')
>>> sams.remove('HG00361') #445中有4个样本不在最终集合里
>>> df1=df[sams]
>>> with open("gene_VNTR.loci") as f:
...     loci=[i.strip() for i in f]
>>> df2=df1.loc[loci]

>>> df2.to_csv("Geuvadis_441_VNTRlen.csv")

# filter VNTRs
python filtVNTRforQTL.py Geuvadis_441_VNTRlen.csv Geuvadis_441_VNTRlen_rmOutliter.csv Geuvadis_445.genes_window100k_vntr.txt

# eVNTR identification with glm
Rscript eVNTR.identification.R
Rscript eVNTR.identification.permuted.R
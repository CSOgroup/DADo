## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
rm(list=ls())
if(!require(DADo))
  devtools::install_github("marzuf/MANUSCRIPT_FIGURES", subdir="DADo")
  # alternatively: 
  # install.packages("DADo_0.0.0.1.tar.gz", repos = NULL, type ="source")
 # data("norm_ID")
library(DADo)
library(doMC)
library(foreach)
nCpu <- 2
registerDoMC(nCpu)

## ----prep_pipeline_data--------------------------------------------------
# list of genes used in the pipeline
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_pipeline_geneList") # this loads pipeline_geneList
# sample IDs
data("luad_ID") # this loads cond2_ID
data("norm_ID") # this loads cond1_ID

# table with gene-to-domain assignment
gene2tad_dt <- read.delim(system.file("extdata", "ENCSR489OCU_NCI-H460_all_genes_positions.txt", package = "DADo"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)

# take only the genes used in the pipeline
pip_g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]

data("ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_fpkmDT") # this loads fpkmDT
fpkm_dt <- fpkmDT
head_sq(fpkm_dt)
# emulate FPKM
# normalize sample-wise
fpkm_dt2 <- apply(fpkm_dt, 2, function(x)x/sum(x))
# stopifnot(colSums(fpkm_dt2) == 1)
stopifnot(abs(colSums(fpkm_dt2) - 1) <= 10^-4)
# and then multiply by 10^6 to have FPKM
fpkm_dt2 <- fpkm_dt2*10^6
fpkm_dt2 <- data.frame(fpkm_dt2, check.names = FALSE)
stopifnot(dim(fpkm_dt) == dim(fpkm_dt2))
stopifnot(rownames(fpkm_dt) == rownames(fpkm_dt2))
stopifnot(colnames(fpkm_dt) == colnames(fpkm_dt2))
fpkm_dt <- fpkm_dt2

stopifnot(names(pipeline_geneList) %in% rownames(fpkm_dt))
stopifnot(!duplicated(pipeline_geneList[names(pipeline_geneList) %in% rownames(fpkm_dt)]))
fpkm_dt <- fpkm_dt[ rownames(fpkm_dt) %in% names(pipeline_geneList),]
rownames(fpkm_dt) <- sapply( rownames(fpkm_dt), function(x) pipeline_geneList[names(pipeline_geneList) == x ])
stopifnot(rownames(fpkm_dt) %in% pip_g2t_dt$entrezID)

## ----prep_purity__data---------------------------------------------------
library(TCGAbiolinks)
head_sq(Tumor.purity)
purityType <- "ESTIMATE"
Tumor.purity$Sample.ID <- as.character(Tumor.purity$Sample.ID)
Tumor.purity$Sample.ID <- substr(Tumor.purity$Sample.ID, start=1, stop=15)
purity_dt <- aggregate(.~Sample.ID, FUN=mean, data=Tumor.purity[c("Sample.ID", purityType)])
mean(cond1_ID %in% purity_dt$Sample.ID) # available samples for cond1
mean(cond2_ID %in% purity_dt$Sample.ID) # available samples for cond2

## ----purityCorr----------------------------------------------------------
tads_all_purity_dt <- get_meanPurityCorr(exprTable=fpkm_dt, purityTable=purity_dt, g2tTable=pip_g2t_dt, 
                                      all_samples=c(cond1_ID, cond2_ID), purityCol="ESTIMATE" )
head(tads_all_purity_dt)

# average at the domain-level
aggTAD_purity_dt <- aggregate(purityCorr~region, FUN=mean, data=tads_all_purity_dt)
head(aggTAD_purity_dt)



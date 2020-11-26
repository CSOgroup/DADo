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

## ----prep_data-----------------------------------------------------------
# table from gene-level DE analysis:
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_fpkmDT") # this loads fpkmDT
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_pipeline_geneList") # this loads pipeline_geneList

# table with gene-to-domain assignment
gene2tad_dt <- read.delim(system.file("extdata", "ENCSR489OCU_NCI-H460_all_genes_positions.txt", package = "DADo"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)

# take only the genes used in the pipeline
pip_g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
head(pip_g2t_dt)

# table with domain positions
tad_pos_dt <- read.delim(system.file("extdata", "ENCSR489OCU_NCI-H460_all_assigned_regions.txt", package = "DADo"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("chromosome", "region", "start", "end"))
stopifnot(pip_g2t_dt$region %in% tad_pos_dt$region)
pip_tad_pos_dt <- tad_pos_dt[tad_pos_dt$region %in% pip_g2t_dt$region,]
head(pip_tad_pos_dt)

## ----do_g2t_permut-------------------------------------------------------
nExprClass <- 5
nPermut <- 10
aggFunc <- "median"

permut_dt <- get_gene2tad_multiPermut(g2TADdt=pip_g2t_dt, 
                                          RNAdt=fpkmDT, 
                                          geneIDlist=pipeline_geneList, 
                                          nClass = nExprClass, 
                                          TADonly=FALSE, # already filtred 
                                          nSimu=nPermut, 
                                          withExprClass=TRUE, 
                                          nCpu=nCpu, 
                                          aggregFun = aggFunc)
head_sq(permut_dt)
rownames(permut_dt) <- permut_dt$entrezID
permut_dt$entrezID <- NULL

stopifnot(ncol(permut_dt) == nPermut)
stopifnot(setequal(rownames(permut_dt), pipeline_geneList))
save(permut_dt, file="package_permut_dt.RData", version=2)

## ----do_acrossBD_permut--------------------------------------------------

sampAcrossBD_data <- getSampleAcrossBD(g2t_DT=pip_g2t_dt, 
                                       tadpos_DT=pip_tad_pos_dt)

str(sampAcrossBD_data[[1]])
stopifnot(length(sampAcrossBD_data) == length(unique(pip_g2t_dt$region)))
save(sampAcrossBD_data, file="package_sampAcrossBD_data.RData", version=2)


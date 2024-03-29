## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
rm(list=ls())

# devtools::install_github("CSOgroup/DADo", subdir="DADo")
# 
# if(!require(DADo))
#   devtools::install_github("CSOgroup/DADo", subdir="DADo")
  # alternatively: 
  # install.packages("DADo_0.0.0.1.tar.gz", repos = NULL, type ="source")
 # data("norm_ID")
library(DADo)

## ----plot_conserved, fig.height=8, fig.width=14--------------------------
data("conserved_region_130_genes_plot_dt") # this loads genes_plot_dt
head(genes_plot_dt)

data("conserved_region_130_tads_plot_dt") # this loads tads_plot_dt
head(tads_plot_dt)

plot_conservedRegions(genes_dt=genes_plot_dt, 
                      tads_dt=tads_plot_dt,
                      dsCat_cols = setNames(c("firebrick3", "navy", "gray50"), c("wt_vs_mut", "norm_vs_tumor", "subtypes")))

## ----plot_volcano, fig.height=6, fig.width=8-----------------------------
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanCorr_TAD") # this loads all_meanCorr_TAD
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanLogFC_TAD") # this loads all_meanLogFC_TAD
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_emp_pval_combined") # this loads emp_pval_combined

plot_volcanoTADsCorrFC(meanCorr=all_meanCorr_TAD, 
                       meanFC=all_meanLogFC_TAD, 
                       comb_pval=emp_pval_combined,
                       tads_to_annot = "chr11_TAD390")

## ----prep_obs------------------------------------------------------------
## Prepare the observed ratioDown data
# table from gene-level DE analysis:
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_DE_topTable") # this loads DE_topTable
DE_topTable$genes <- as.character(DE_topTable$genes)
# list of genes used in the pipeline
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_pipeline_geneList") # this loads pipeline_geneList
# for those genes, I have logFC data
stopifnot(names(pipeline_geneList) %in% DE_topTable$genes)
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(pipeline_geneList),]
DE_topTable$entrezID <- pipeline_geneList[DE_topTable$genes]
stopifnot(!is.na(DE_topTable$entrezID))
# table with gene-to-domain assignment
gene2tad_dt <- read.delim(system.file("extdata", "ENCSR489OCU_NCI-H460_all_genes_positions.txt", package = "DADo"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
# take only the genes used in the pipeline
pip_g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
# merge to match gene-to-TAD and logFC
merged_dt <- merge(pip_g2t_dt[,c("entrezID", "region")], DE_topTable[,c("logFC", "entrezID")], by="entrezID", all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(merged_dt))
# compute the ratioDown for each domain (cf. domain-level statistics vignette !)
ratioDown_dt <- aggregate(logFC ~ region, FUN=get_ratioDown, data=merged_dt)
colnames(ratioDown_dt)[colnames(ratioDown_dt) == "logFC"] <- "ratioDown"
obs_ratioDown <- setNames(ratioDown_dt$ratioDown, ratioDown_dt$region)
save(obs_ratioDown, file="package_obs_ratioDown.RData", version=2)

## ----ratioDown_permut_example--------------------------------------------
require(foreach)
## Prepare the ratioDown for permutation data
data("cut1000_ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_permutationsDT") # this loads permutationsDT
head_sq(permutationsDT)
stopifnot(setequal(rownames(permutationsDT), pipeline_geneList))
tad_levels <- as.character(ratioDown_dt$region)
all_permut_ratioDown_dt <- foreach(i = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
  perm_g2t <- data.frame(entrezID=as.character(rownames(permutationsDT)),
                         region = as.character(permutationsDT[,i]),
                         stringsAsFactors = FALSE)
  perm_merged_dt <- merge(perm_g2t, DE_topTable[,c("logFC", "entrezID")], by="entrezID", all.x=TRUE, all.y=FALSE)
  stopifnot(!is.na(perm_merged_dt))
  # compute the FCC for each domain
  permut_rD_dt <- aggregate(logFC ~ region, FUN=get_ratioDown, data=perm_merged_dt)
  colnames(permut_rD_dt)[colnames(permut_rD_dt) == "logFC"] <- "ratioDown"
  rownames(permut_rD_dt) <- as.character(permut_rD_dt$region)
  stopifnot(setequal(rownames(permut_rD_dt), tad_levels))
  out_dt <- permut_rD_dt[tad_levels, "ratioDown", drop=FALSE]
  stopifnot(out_dt$ratioDown >= 0 & out_dt$ratioDown <= 1) # just to check...
  out_dt
}
head_sq(all_permut_ratioDown_dt)

## ----results="hide", fig.show="hide"-------------------------------------
## Do the plots:
myplots <- plot_smileRatioDownConcord(observed_rD=obs_ratioDown,
                           permut_rD_dt=all_permut_ratioDown_dt,
                           plotTit1=paste0("NCI-H460 - TCGA normal vs. LUAD")
                           )

## ----plot_line, fig.height=6, fig.width=6--------------------------------
myplots[[1]]


## ----plot_hist, fig.height=6, fig.width=6--------------------------------
myplots[[2]]


#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.RData
# - script0: rna_geneList.RData
# - script0: pipeline_geneList.RData
# - script0: rna_madnorm_rnaseqDT.RData
################################################################################

################  OUTPUT
# - all_meanCorr_TAD.RData
################################################################################

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")

script1_name <- "1_prepGeneData"
script_name <- "4_runMeanTADCorr"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  # for toTitleCase
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu) # from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

#*************************************************************** !!! HARD-CODED SETTINGS 
# imported from setting file
txt <- paste0(toupper(script_name), "> Take the diagonal when computing the mean of lower triangle: ", as.character(withDiag), "\n")
printAndLog(txt, pipLogFile)

# imported from setting file
txt <- paste0(toupper(script_name), "> Correlation method: ",  toTitleCase(corrMethod), "\n")
printAndLog(txt, pipLogFile)

################################****************************************************************************************
################################********************************************* PREPARE INPUT
################################****************************************************************************************
#*******************************************************************************
# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
if(microarray) {
  norm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script1_name, "rna_madnorm_rnaseqDT.RData"))))
} else{
  norm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script1_name, "rna_qqnorm_rnaseqDT.RData")))) 
}
initList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "rna_geneList.RData"))))
geneList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_geneList.RData"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

norm_rnaseqDT <- norm_rnaseqDT[names(geneList),]    
stopifnot(all(rownames(norm_rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
# ADDED FOLLOWING LINE 02.07.2019, just to be sure
stopifnot(as.character(geneList) %in% gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_regionList.RData"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

all_regions <- unique(as.character(gene2tadDT$region))

################################****************************************************************************************
################################********************************************* iterate over regions to compute intraTAD correlation
################################****************************************************************************************

if(useTADonly) {
  initLen <- length(all_regions)
  all_regions <- all_regions[grep("_TAD", all_regions)]
  txt <- paste0(toupper(script_name), "> Take only the TAD regions: ", length(all_regions),"/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

cat(paste0("... start intra TAD correlation\n"))

all_meanCorr_TAD <- foreach(reg=all_regions, .combine='c') %dopar% {
  # cat(paste0("... doing region: ", reg, "/", length(all_regions), "\n"))
  reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]

# NB 02.07.2019 -> this is ok because I have reordered norm_rnaseqDT : norm_rnaseqDT <- norm_rnaseqDT[names(geneList),]    

  rowsToKeep <- which(geneList %in% reg_genes)
  subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
  stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
  # columns => the genes
  # rows => the samples
  ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
  ### UPDATE: this should not happen in the latest version !!!
  # stopifnot(ncol(subData) == length(reg_genes))
  # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
  stopifnot(ncol(subData) == length(reg_genes))
  #### CORRELATION
  corrMatrix_all <- cor(subData, method = corrMethod)
  # should be correlation of the genes
  ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
  ### UPDATE: this should not happen in the latest version !!!
  # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
  # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
  stopifnot(nrow(corrMatrix_all) == length(reg_genes))
  stopifnot(ncol(corrMatrix_all) == length(reg_genes))
  mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
}

cat(paste0("... end intra TAD correlation\n"))

names(all_meanCorr_TAD) <- all_regions
stopifnot(length(all_meanCorr_TAD) == length(all_regions))

save(all_meanCorr_TAD, file= file.path(curr_outFold, "all_meanCorr_TAD.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "all_meanCorr_TAD.RData"), "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



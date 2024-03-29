#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.RData
# - script0: pipeline_geneList.RData
# - script1: DE_topTable.RData
# - script1: DE_geneList.RData
################################################################################

################  OUTPUT
# - all_FCC_TAD.RData
################################################################################

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(".")

script1_name <- "1_prepGeneData"
script2_name <- "2_runGeneDE"
script_name <- "11_runFCC"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))

registerDoMC(nCpu) # loaded from main_settings.R

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

################################***********************************************************************************
############ LOAD INPUT DATA
################################***********************************************************************************
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(file.path(pipOutFold, script2_name, "DE_topTable.RData"))))
DE_geneList <- eval(parse(text = load(file.path(pipOutFold, script2_name, "DE_geneList.RData"))))

pipeline_geneList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_geneList.RData"))))
pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_regionList.RData"))))

if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}

stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))

entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
names(entrezList) <- DE_topTable$genes
stopifnot(length(entrezList) == length(DE_topTable$genes))

# replace the gene symbol rownames by ensemblID rownames
logFC_DT <- data.frame(entrezID =  entrezList,
                       logFC = DE_topTable$logFC, stringsAsFactors = F)

rownames(logFC_DT) <- NULL
initNrow <- nrow(logFC_DT)
logFC_DT <- logFC_DT[logFC_DT$entrezID %in% pipeline_geneList,]
txt <- paste0(toupper(script_name), "> Take only filtered genes: ", nrow(logFC_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

### take only the filtered data according to initial settings
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(pipeline_geneList),]
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

################################***********************************************************************************
################################********************************************* get observed FCC for all regions
################################***********************************************************************************

cat(paste0("... start computing FCC by TAD \n"))

head(logFC_DT)

mergedDT <- merge(logFC_DT, gene2tadDT[,c("entrezID", "region")], by="entrezID", all.x=TRUE, all.y=FALSE)

stopifnot(nrow(mergedDT) == nrow(na.omit(mergedDT)))
stopifnot(!is.na(mergedDT))
stopifnot(nrow(mergedDT) == length(pipeline_geneList))

fcc_DT <- aggregate(logFC ~ region, data=mergedDT, FUN=get_fcc)
all_FCC_TAD <- setNames(fcc_DT$logFC, fcc_DT$region)
stopifnot(length(all_FCC_TAD) == length(unique(gene2tadDT$region)))
txt <- paste0(toupper(script_name), "> Number of regions for which FCC computed: ", length(all_FCC_TAD), "\n")
printAndLog(txt, pipLogFile)

stopifnot(setequal(pipeline_regionList, names(all_FCC_TAD)))

if(useTADonly) {
    initLen <- length(all_FCC_TAD)
    all_FCC_TAD <- all_FCC_TAD[grep("_TAD", names(all_FCC_TAD))]
    txt <- paste0(toupper(script_name), "> Take only the TAD regions: ", length(all_FCC_TAD),"/", initLen, "\n")
    printAndLog(txt, pipLogFile)
}

save(all_FCC_TAD, file= file.path(curr_outFold, "all_FCC_TAD.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "all_FCC_TAD.RData"), "\n"))


################################***********************************************************************************
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))





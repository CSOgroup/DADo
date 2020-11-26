#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script1: pipeline_regionList.RData
# - script1: pipeline_geneList.RData
# - script1: rangeTADgenes.RData
# - script2: DE_topTable.RData
# - script2: DE_geneList.RData
# - script5fc: permutationsDT.RData
################################################################################

################  OUTPUT
# FCC_permDT.RData
################################################################################

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")

script1_name <- "1_prepGeneData"
script2_name <- "2_runGeneDE"
script5fc_name <- "5fc_runPermutationsMedian"
script_name <- "12_runPermutationsFCC"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu) # from main_settings.R

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=FALSE, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(file.path(pipOutFold, script2_name, "DE_topTable.RData"))))
DE_geneList <- eval(parse(text = load(file.path(pipOutFold, script2_name, "DE_geneList.RData"))))

stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))

entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
stopifnot(length(entrezList) == length(DE_topTable$genes))

# replace the gene symbol rownames by ensemblID rownames
logFC_DT <- data.frame(entrezID =  entrezList,
                       logFC = DE_topTable$logFC, stringsAsFactors = F)
rownames(logFC_DT) <- NULL

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% entrezList,]

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(file.path(pipOutFold, script5fc_name, "permutationsDT.RData"))))
if(ncol(permutationsDT) != nRandomPermut)
  stop(paste0("! NEED TO CHECK: different settings were used for running the permutations !\nncol(permutationsDT) != nRandomPermut\nncol(permut_FCC)\t=\t", ncol(permutationsDT),"\nnRandomPermut\t=\t", nRandomPermut, "\n"))

pipeline_geneList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_geneList.RData"))))
if(!setequal(pipeline_geneList, rownames(permutationsDT))) {
  txtWarningGene <- paste0(toupper(script_name), "> Not the same set of genes in permutDT and pipeline_geneList\n")
  printAndLog(txtWarningGene, pipLogFile)    
  stop(txtWarningGene)
} else{
  txtWarningGene <- ""
}

pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_regionList.RData"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
if(!setequal(pipeline_regionList, permutationsDT[,2])) {
  txtWarningRegion <- paste0(toupper(script_name), "> Not the same set of regions in permutDT and pipeline_regionList\n")
  printAndLog(txtWarningRegion, pipLogFile)    
  stop(txtWarningRegion)
} else {
  txtWarningRegion <- ""
}

### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
initNrow <- nrow(permutationsDT)
permutationsDT <- permutationsDT[which(rownames(permutationsDT) %in% logFC_DT$entrezID),]
txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(permutationsDT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

all_regions <- sort(unique(as.character(permutationsDT[,2])))

if(useTADonly) {
  if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
    stop("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
  initLen <- length(all_regions)
  all_regions <- all_regions[grep("_TAD", all_regions)]
  # if want to use only TADs, make more sens if permutations were run without the TADs
  txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)    
  if(length(all_regions) < initLen){
    stop("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
}

################################****************************************************************************************
######################################################## COMPUTE MEAN LOG FC BY TAD FOR THE PERMUTATIONS
################################****************************************************************************************

### REGIONS ARE STORED IN ROWNAMES OF PERMUTdt

cat("... start FCC permutDT \n")
fcc_permDT <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
  cat(paste0("...TAD FCC for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))  
  
  g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                           region = permutationsDT[,i_col], stringsAsFactors = F)
  g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
  g2t_permDT$region <-  as.character(g2t_permDT$region)

  unlist(sapply(unique(all_regions), function(x) {
    reg_genes <- g2t_permDT$entrezID[which(g2t_permDT$region == x)]
    # head(reg_genes)
    stopifnot(reg_genes %in% logFC_DT$entrezID)
	fc_vect <- logFC_DT$logFC[logFC_DT$entrezID %in% reg_genes]
	stopifnot(is.numeric(fc_vect))
	stopifnot(!is.na(fc_vect))
    get_fcc(fc_vect)
  }))
}

cat("... end FCC permutDT \n")

FCC_permDT <- as.data.frame(fcc_permDT)
stopifnot(ncol(FCC_permDT) == ncol(permutationsDT))
colnames(FCC_permDT) <- paste0("permutation",  c(1:ncol(permutationsDT)))
stopifnot(nrow(FCC_permDT) == length(all_regions))
rownames(FCC_permDT) <- all_regions

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

txt <- paste0(toupper(script_name), "> Number of permutations for which FCC computed: ", ncol(FCC_permDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of regions for which FCC computed: ", nrow(FCC_permDT), "\n")
printAndLog(txt, pipLogFile)

#save(FCC_permDT, file= paste0(curr_outFold, "/FCC_permDT.RData"))
# update 16.08.2019 => faster save version
my_save.pigz(FCC_permDT, pigz_exec_path = pigz_exec_path, file= file.path(curr_outFold, "FCC_permDT.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "FCC_permDT.RData"), "\n"))

cat(paste0(txtWarningGene))
cat(paste0(txtWarningRegion))

################################****************************************************************************************
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



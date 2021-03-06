#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

#### UPDATE: do not take raw counts but fpkm data !!!

set.seed(20180202) # this row was added 08.03.18, the files in OUTPUTFOLDER so far without set.seed => but not reproducible on multiple cores ?? cf. trial of the 16.08.2019

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.RData
# - script0: rna_geneList.RData
# - script0: pipeline_geneList.RData
# - script0: rna_madnorm_rnaseqDT.RData
# - script0: rna_fpkmDT.RData # UPDATE 
################################################################################

################  OUTPUT
# - permutationsDT.RData
################################################################################

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")

script1_name <- "1_prepGeneData"
script_name <- "5fc_runPermutationsMedian"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))

#source(file.path(pipScriptDir, "TAD_DE_utils_fasterPermut.R")) # UPDATE 16.08.2019 -> modified function for tad shuffling => faster permuts

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

nRandom <- nRandomPermut

if(withExprClass)
  nClass <- permutExprClass  # number of class of expression


# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> withExprClass\t=\t", withExprClass, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> nClass\t=\t", nClass, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
  
#******************************************************************************* !! HARD CODED
aggregFction <- "median"
# withExprClass <- TRUE # loaded from main_settings.R !
#*******************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
if(withExprClass) {
    # rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/rna_rnaseqDT.RData"))))
    rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script1_name, "rna_fpkmDT.RData"))))
}
initList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "rna_geneList.RData"))))
geneList <- eval(parse(text = load(file.path(pipOutFold, script1_name, "pipeline_geneList.RData"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))
#*******************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
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

if(useTADonly) {
    initLen <- length(geneList)
    TAD_gene2tadDT <- gene2tadDT[grep("_TAD",gene2tadDT$region),]
    stopifnot(all(rownames(rnaseqDT) == names(geneList)))
    rowsToKeep <- which(geneList %in% as.character(TAD_gene2tadDT$entrezID))
    geneList <- geneList[rowsToKeep]
    rnaseqDT <- rnaseqDT[rowsToKeep,]
    gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]
    txt <- paste0(toupper(script_name), "> Take only genes that are within TADs: ", length(geneList), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)
}

################################****************************************************************************************
####################################################### RUN PERMUTATIONS
################################****************************************************************************************

cat("... Start permutations\n")
if(withExprClass) {
    shuffleData <- get_multiShuffledPositions_vFunct(g2TADdt=gene2tadDT, RNAdt=rnaseqDT, 
                                              geneIDlist=geneList, nClass = nClass, TADonly=F, nSimu=nRandom, 
                                              withExprClass=withExprClass, nCpu=nCpu, aggregFun = aggregFction)
} else {
    shuffleData <- get_multiShuffledPositions_vFunct(g2TADdt=gene2tadDT, RNAdt=NULL, 
                                              geneIDlist=geneList, nClass = NULL, TADonly=F, nSimu=nRandom,
                                              withExprClass=withExprClass, nCpu=nCpu, aggregFun = aggregFction)
}

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************
permutationsDT <- shuffleData
rownames(permutationsDT) <- permutationsDT[,1]
permutationsDT <- permutationsDT[,-1]

txt <- paste0(toupper(script_name), "> Number of permutations: ", nRandom, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> With expression classes: ", as.character(withExprClass), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> WARNING: genes that map to duplicated entrezID are removed ! \n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> -> number of genes retained for permutations: ", nrow(permutationsDT), "/", length(geneList), "\n")
printAndLog(txt, pipLogFile)

#save(permutationsDT, file = paste0(curr_outFold, "/permutationsDT.RData"))
# update 16.08.2019 => faster save version
my_save.pigz(permutationsDT, pigz_exec_path = pigz_exec_path, file = file.path(curr_outFold, "permutationsDT.RData") )
cat(paste0("... written: ", file.path(curr_outFold, "permutationsDT.RData"), "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

stopifnot(ncol(permutationsDT) == (nRandom))

cat("dim(permutationsDT)\n")
cat(dim(permutationsDT))
cat("\n")

cat("... using faster save and permut\n")


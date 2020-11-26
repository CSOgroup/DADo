#!/usr/bin/Rscript

### -> UPDATE 13.08.2019
# => retain the genes for which <min_sampleRatio*100>% of the samples have at least <min_counts> counts 
# => default settings: min_sampleRatio = 0.8 and min_counts = 5
# (cf. plot ratio genes ~ ratio samples in 2_Yuanlong_Cancer_HiC_data_TAD_DA/)
# (settings in run_pipeline.sh -> written in the setting file)

cat(paste0("> START ", "1_prepGeneData",  "\n"))

startTime <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")

# do not change the path of script_name, is used for the output folder name
script_name <- paste0("1_prepGeneData")
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

# cat(paste0("setDir = ", setDir, "\n"))
cat("source main_settings \n")
source("main_settings.R") # setDir is the main_settings not in run_settings
cat("source settingF \n")
source(settingF)
source(file.path(pipScriptDir,"TAD_DE_utils.R"))
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

plotType <- "svg"

stopifnot(exists("inputDataType"))
stopifnot(inputDataType %in% c("raw", "RSEM", "FPKM", "DESeq2", "microarray"))
stopifnot(inputDataType == "RSEM") # for release july 2020

cat(paste0("> inputDataType: ", inputDataType, "\n"))

######### 0_prepGene
# ->  according to the user settings, create a geneList object 
# ... with the same genes as rownames in input data
# ... geneList["geneName"] = entrezID
# -> added 21.02.2018 -> output the FPKM table to create the expression in the permutations
#####################

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# load the expression data
if(inRdata) {
#  x <- load(paste0(setDir, rnaseqDT_file))
  x <- load(paste0(rnaseqDT_file)) # for release July 2020
  rnaseqDT <- eval(parse(text=x))
} else {
  if(geneID_loc == "rn"){
#    rnaseqDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T, row.names = 1)
    rnaseqDT <- read.delim(file.path(rnaseqDT_file), sep=my_sep, header=T, row.names = 1) # for release July 2020
  } else {
#    rnaseqDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T)
    rnaseqDT <- read.delim(file.path(rnaseqDT_file), sep=my_sep, header=T) # for release July 2020
    rownames(rnaseqDT) <- rnaseqDT[, geneID_loc]
    rnaseqDT <- rnaseqDT[,-geneID_loc]
  }
}

refGenes <- as.character(rownames(rnaseqDT))
cat(paste0("... initial number of rows (genes): ", length(refGenes), "\n"))

# ADDED 13.08.2019 to change the way counts are filtered
# set in run_pipeline.sh -> written in the setting file #min_sampleRatio <- 0.8 #min_counts <- 5
stopifnot(exists("min_sampleRatio"))
stopifnot(exists("min_counts"))
stopifnot(is.numeric(min_counts))
stopifnot(is.numeric(min_sampleRatio))
stopifnot(min_sampleRatio >= 0 & min_sampleRatio <= 1)

# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> min_sampleRatio\t=\t", min_sampleRatio, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> min_counts\t=\t", min_counts, "\n")
printAndLog(txt, pipLogFile)
### NA ARE NOT ALLOWED WHEN COMPUTING CPM, CHECK THE RNASEQDT DO NOT CONTAIN NA
txt <- paste0(toupper(script_name), "> replace NA with 0 before computing cpm: ", 
              sum(is.na(rnaseqDT)), "/", dim(rnaseqDT)[1]*dim(rnaseqDT)[2], " (",
              round((sum(is.na(rnaseqDT))/(dim(rnaseqDT)[1]*dim(rnaseqDT)[2]) * 100),2), "%)\n")
printAndLog(txt, pipLogFile)
rnaseqDT[is.na(rnaseqDT)] <- 0

stopifnot(nrow(rnaseqDT) > 0)

########################################################################################
# FILTER 1: FILTER GENES FOR WHICH POSITION IS AVAILABLE AND REMOVE THE ONES WITH DUPLICATED ENTREZ ID
# (mandatory filter)
########################################################################################

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

# could not be duplicated in refGenes because it is the rownames

cat("... filter the genes for which position information is available\n")
# gene2tadDT is used here for filtering position information, not related to the set of TADs used !!!
if(geneID_format == "entrezID"){
  rna_geneList <- get_geneList_fromEntrez(refList=refGenes, g2t=gene2tadDT, histDT_file=historyDT_file)
} else if(geneID_format == "ensemblID"){
  rna_geneList <- get_geneList_fromEnsembl(refList=refGenes, g2t=gene2tadDT, ensDT_file=ensemblDT_file)  
} else if(geneID_format == "geneSymbol"){
  rna_geneList <- get_geneList_fromSymbol(refList=refGenes, g2t=gene2tadDT, symbDT_file=symbolDT_file)
} else {
  stop("should never happen")
}

head(refGenes)
head(gene2tadDT)

stopifnot(length(rna_geneList) > 0)

# check the returned list
stopifnot(all(rna_geneList %in% gene2tadDT$entrezID)) # should be TRUE otherwise no information about their position !
stopifnot(all(names(rna_geneList) %in% refGenes))
stopifnot(!any(duplicated(rna_geneList)))
stopifnot(!any(duplicated(names(rna_geneList))))

txt <- paste0(toupper(script_name), "> Number of genes retained: ", length(rna_geneList), "/", length(refGenes), "\n")
printAndLog(txt, pipLogFile)
# discard from the rnaseqDT the rows for which I have no information
rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% names(rna_geneList),]
stopifnot(all(rownames(rnaseqDT) == names(rna_geneList)))

# in the new version, I systematically discard the ambiguous ones
if(removeDupGeneID) {
  cat("... remove ambiguous genes with duplicated entrezID\n")
  initLen <- length(rna_geneList)
  stopifnot(initLen == nrow(rnaseqDT))
  dupEntrez <- unique(as.character(rna_geneList[which(duplicated(rna_geneList))]))
  rna_geneList <- rna_geneList[!rna_geneList %in% dupEntrez]
  rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% names(rna_geneList),]
  stopifnot(all(rownames(rnaseqDT) == names(rna_geneList)))
  stopifnot(length(rna_geneList) == nrow(rnaseqDT))
  txt <- paste0(toupper(script_name), "> Remove geneID with ambiguous entrezID, final # of rows  : ", length(rna_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
  stopifnot(!any(duplicated(rna_geneList)))
}

##### TAKE ONLY THE SAMPLES OF INTEREST
#samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
#samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
samp1 <- eval(parse(text=load(file.path(sample1_file)))) # for release July 2020
samp2 <- eval(parse(text=load(file.path(sample2_file)))) # for release July 2020
rnaseqDT <- rnaseqDT[, c(samp1, samp2)]
stopifnot( ncol(rnaseqDT) == (length(samp1) + length(samp2)) )

# rna_geneList will hold the list of genes for which I have 1) genomic positions 2) expression data
save(rna_geneList, file = file.path(curr_outFold, "rna_geneList.RData"))
cat(paste0("... written: ", file.path(curr_outFold,  "rna_geneList.RData"), "\n"))

rna_rnaseqDT <- rnaseqDT
stopifnot(nrow(rna_rnaseqDT) > 0)
save(rna_rnaseqDT, file = file.path(curr_outFold, "rna_rnaseqDT.RData"))
cat(paste0("... written: ", file.path(curr_outFold,  "rna_rnaseqDT.RData"), "\n"))

#### PREPARE THE QQNORM DATA FOR THE GENES I USED FOR DE ANALYSIS
cat("... qqnorm the data for other analyses \n")
rna_qqnorm_rnaseqDT <- t(apply(rna_rnaseqDT, 1, quantNorm))
stopifnot(all(dim(rna_qqnorm_rnaseqDT) == dim(rna_rnaseqDT)))
rownames(rna_qqnorm_rnaseqDT) <- rownames(rna_rnaseqDT)
colnames(rna_qqnorm_rnaseqDT) <- colnames(rna_rnaseqDT)
save(rna_qqnorm_rnaseqDT, file = file.path(curr_outFold,  "rna_qqnorm_rnaseqDT.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "rna_qqnorm_rnaseqDT.RData"), "\n"))

# => rna_geneList: all the genes for which I have position information
# => countFilter_geneList: all the genes for which I have position information and that are >= CPM threshold

pipeline_geneList <- rna_geneList
stopifnot(length(pipeline_geneList) > 0)

########################################################################################
# FILTER 2: FILTER GENES THAT PASS MIN CPM THRESHOLD
# (only if useFilterCountData == TRUE)
########################################################################################

### !!! changed here 13.08.2019
#keep genes for which at least min_sampleRatio of the samples have >= min_counts

if(useFilterCountData) {
  stopifnot(is.numeric(rna_rnaseqDT[1,1]))
  # take only the genes for which I have position (filter 1)
  countFilter_rnaseqDT <- rna_rnaseqDT[which(rownames(rna_rnaseqDT) %in% names(pipeline_geneList)),]
  stopifnot(all(rownames(countFilter_rnaseqDT) == names(pipeline_geneList)))
  # ensure the samples are present in the column names
  stopifnot(all(samp1 %in% colnames(countFilter_rnaseqDT)))
  stopifnot(all(samp2 %in% colnames(countFilter_rnaseqDT)))
  countFilter_rnaseqDT <- countFilter_rnaseqDT[,c(samp1, samp2)]

  totSamples <- length(samp1) + length(samp2)
  stopifnot(ncol(countFilter_rnaseqDT) == totSamples)
  
  if(inputDataType == "raw" | inputDataType == "RSEM") {

    # CHANGE: at least min_counts reads in at least min_sampleRatio of samples
    # for each gene, how many samples have enough read (number)
    genes_nSamples_atLeastMinReads <- apply(countFilter_rnaseqDT, 1, function(x) sum(x >= min_counts)) 
    stopifnot(names(genes_nSamples_atLeastMinReads) == rownames(countFilter_rnaseqDT))
    # convert # of samples to ratio of samples
    genes_ratioSamples_atLeastMinReads <- genes_nSamples_atLeastMinReads/totSamples
    stopifnot(genes_ratioSamples_atLeastMinReads >= 0)
    stopifnot(genes_ratioSamples_atLeastMinReads <= 1)
    stopifnot(length(genes_ratioSamples_atLeastMinReads) == nrow(countFilter_rnaseqDT))
    # which genes have a ratio of samples with at least min counts above the min. ratio of samples
    rowsToKeep <- genes_ratioSamples_atLeastMinReads >= min_sampleRatio
    stopifnot(length(rowsToKeep) == nrow(countFilter_rnaseqDT))


    stopifnot(names(pipeline_geneList) == rownames(countFilter_rnaseqDT))
    keptRatio <- sum(rowsToKeep)/length(pipeline_geneList)
    txt <- paste0(toupper(script_name), "> useFilterCountData is TRUE -> minCount-filtered geneList; to keep: ", sum(rowsToKeep), "/", length(pipeline_geneList), "(", round(keptRatio*100,2), "%)\n")
    printAndLog(txt, pipLogFile)
  } else {
    stop("ERROR\n")
  }

  pipeline_geneList <- pipeline_geneList[rowsToKeep]
  stopifnot(length(pipeline_geneList) == sum(rowsToKeep))
}

stopifnot(length(pipeline_geneList) > 0)

########################################################################################
# FILTER 3: FILTER TO USE ONLY TAD REGIONS (AND ONLY GENES BELONGING TO) NOT INTER-TAD
# (only if useTADonly == TRUE)
########################################################################################
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT_intact <- gene2tadDT

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]

stopifnot(length(pipeline_geneList) == nrow(gene2tadDT))

initRegionsLen <- length(unique(gene2tadDT$region))

#setequal
stopifnot(all(gene2tadDT$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT$entrezID))

if(useTADonly) {
  initNrow <- nrow(gene2tadDT)
  gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> only TAD regions: ", nrow(gene2tadDT), "/", initNrow, " regions\n")
  printAndLog(txt, pipLogFile)
  initLen <- length(pipeline_geneList)
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT$entrezID]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> TAD-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

########################################################################################
# BEFORE FILTER 4: LOOK AT TAD GENE NUMBER DISTRIBUTION
# (only if useFilterSizeData == TRUE)
########################################################################################

tadGeneNbr_DT <- aggregate(entrezID ~ region, data=gene2tadDT, FUN=length)
colnames(tadGeneNbr_DT) <- c("regionID", "nbrGenes") 
tadGeneNbr_DT$region <- ifelse(regexpr("_TAD", tadGeneNbr_DT$regionID) > 0, "TAD", "BOUND")
tadGeneNbr_DT$nbrLog10 <- log10(tadGeneNbr_DT$nbrGenes)

if(useTADonly) {
  stopifnot(all(tadGeneNbr_DT$region == "TAD"))
}

p1 <- ggdensity(tadGeneNbr_DT, x = "nbrLog10",
          title ="Distribution log10(# genes) before filtering",
          # add = "mean",
          rug = TRUE,
          xlab="log10 # of genes/region",
          color = "region", fill = "region",
          palette = c("#00AFBB", "#E7B800"))
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))

if(! useFilterSizeData) {
  maxQuantGeneTAD <- 1
  minNbrGeneTAD <- 0
}
upperLimit <- as.numeric(quantile(tadGeneNbr_DT$nbrGenes, probs=maxQuantGeneTAD))
txt <- paste0(toupper(script_name), "> Current threhsold for TAD size: >= ", minNbrGeneTAD, " and <= ", upperLimit, " (", round(maxQuantGeneTAD*100, 2), "%)\n")
printAndLog(txt, pipLogFile)
tadGeneNbr_DT <- tadGeneNbr_DT[tadGeneNbr_DT$nbrGenes >= minNbrGeneTAD & tadGeneNbr_DT$nbrGenes <= upperLimit,]

# ADDED 15.08.2019 to easily retrieve the tad size limits
outFile <- file.path(curr_outFold, "gene_nbr_filter.RData")
gene_nbr_filter <- c(minNbrGeneTAD, maxQuantGeneTAD)
save(gene_nbr_filter, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

rangeTADgenes <- c(minNbrGeneTAD, upperLimit)
save(rangeTADgenes, file = file.path(curr_outFold,  "rangeTADgenes.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "rangeTADgenes.RData"), "\n"))

p2 <- ggdensity(tadGeneNbr_DT, x = "nbrLog10",
          title ="Distribution log10(# genes) after filtering",
          # add = "mean", 
          rug = TRUE,
          xlab="log10 # of genes/region",
          color = "region", fill = "region",
          palette = c("#00AFBB", "#E7B800"))
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))

p1_p2 <- ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

ggsave(filename = paste0(curr_outFold, "/", "TADsize_density.", plotType), plot = p1_p2, height=7, width=10)


########################################################################################
# FILTER 4: FILTER TAD SIZE TO USE ONLY TAD (AND ONLY GENES BELONGING TO) >= LOWER AND <= UPPER BOUND
# (only if useFilterSizeData == TRUE)
########################################################################################

#setequal
stopifnot(all(gene2tadDT$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT$entrezID))

if(useFilterSizeData){
  gene2tadDT <- gene2tadDT[gene2tadDT$region %in% as.character(tadGeneNbr_DT$regionID),]
  initLen <- length(pipeline_geneList)
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT$entrezID]
  txt <- paste0(toupper(script_name), "> useFilterSizeData is TRUE -> size-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

txt <- paste0(toupper(script_name), "> pipeline_geneList compared to available genes: ", length(pipeline_geneList), "/", length(rna_geneList), "\n")
printAndLog(txt, pipLogFile)
save(pipeline_geneList, file = file.path(curr_outFold, "pipeline_geneList.RData"))
cat(paste0("... written: ", file.path(curr_outFold,  "pipeline_geneList.RData"), "\n"))

pipeline_regionList <- unique(as.character(gene2tadDT$region))
txt <- paste0(toupper(script_name), "> pipeline_regionList compared to available regions: ", length(pipeline_regionList), "/", initRegionsLen, "\n")
printAndLog(txt, pipLogFile)
save(pipeline_regionList, file = file.path(curr_outFold,  "pipeline_regionList.RData"))
cat(paste0("... written: ", file.path(curr_outFold, "pipeline_regionList.RData"), "\n"))


########################################################################################
# added 21.02 => create a FPKM data table ! [if not already fpkm provided]
########################################################################################
rna_rnaseqDT <- eval(parse(text = load(file.path(curr_outFold,  "rna_rnaseqDT.RData"))))
rna_geneList <- eval(parse(text = load(file.path(curr_outFold,  "rna_geneList.RData"))))
stopifnot(length(rna_geneList) == nrow(rna_rnaseqDT))
stopifnot(all(names(rna_geneList) == rownames(rna_rnaseqDT)))
stopifnot(all(rna_geneList %in% entrezDT$entrezID))

# Giovanni: RSEM already corrected for gene length !

#("raw", "RSEM", "FPKM", "DESeq2", "microarray")

###**** ADDED HERE FOLLOWING DISCUSSION WITH MARCO RSEM CONDITION  --- 07.12.2018
if(inputDataType == "RSEM") {
  cat("... RSEM data: fpkm file provided in setting file\n")  
  stopifnot(exists("rna_fpkmDT_file"))
  loaded_rna_fpkmDT <- eval(parse(text = load(rna_fpkmDT_file)))
  stopifnot(rownames(rna_rnaseqDT) %in% rownames(loaded_rna_fpkmDT))
  stopifnot(colnames(rna_rnaseqDT) %in% colnames(loaded_rna_fpkmDT))
  rna_fpkmDT <- loaded_rna_fpkmDT[rownames(rna_rnaseqDT), colnames(rna_rnaseqDT)]
} else {
  stop("error\n")
}

stopifnot(dim(rna_fpkmDT) == dim(rna_rnaseqDT))

save(rna_fpkmDT, file = file.path(curr_outFold,  "rna_fpkmDT.RData"))
cat(paste0("... written: ", file.path(curr_outFold,  "rna_fpkmDT.RData"), "\n"))


################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


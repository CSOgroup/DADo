cat("load main_settings.R\n")

########################################################
### set paths to main files needed accross the workflow + general settings
########################################################

historyDT_file <- paste0("/mnt/ed4/marie/entrez2synonym/entrez/data/gene_history_reordered.txt")

# file with mapping from entrez to chromosomic positions
#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.RData") # previous also held symbols
# holds only netrezID and positions
entrezDT_file <- paste0("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")

# file with symbol synonyms to entrezID
#synoDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/all_entrez_syno_1.RData")
# mapping entrezID and all possible symbols
symbolDT_file <- paste0("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_SYMBOL/final_entrez2syno.txt")

# mapping entrezID and ensemblID
ensemblDT_file <- paste0("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_ENSEMBL/final_entrez2ensembl.txt")

nCpu <- 30

########################################################
### which kind of data to prepare and use across pipeline
########################################################

# if TRUE, from scripts 2_ to further, use only the genes that pass the DE analysis threshold for all analyses
useFilterCountData <- TRUE

# use only TADs (and corresponding genes)>= minSize et < maxQuantile and genes belonging to those TADs
useFilterSizeData <- TRUE

# if TRUE, from scripts 2_ to further, only focus on  TAD regions
useTADonly <- TRUE


########################################################
### TAD filter number of genes - for filter data (above) 
########################################################

# keep TAD with >= minNbrGeneTAD
minNbrGeneTAD <- 3

# keep TAD with <= maxQuant Genes 
maxQuantGeneTAD <- 0.99


########################################################
### for computing correlations - 4_runMeanTADCorr.R and 7_runPermutationsMeanTADCorr.R
########################################################
 
# which correlation method to use for computing correlations
corrMethod <- "pearson"

# should include correlation with itself (step 4 only)
withDiag <- FALSE



########################################################
### permutation settings - 5fc_runPermutationsMedian.R        => STEP 5fc
########################################################

withExprClass <- TRUE
# number of classes of expression used for permutate expression data
permutExprClass <- 5


########################################################
### for fast saving in step 5fc and 6
########################################################

pigz_exec_path <- file.path("/mnt/ed4/marie/pigz-2.4/pigz")


########################################################
### matching pattern to retrieve permutation correlation files step 9
########################################################
corrMatchPattern <- "meanCorr_sample_around_TADs_sameNbr.RData" 

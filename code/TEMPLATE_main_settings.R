cat("load main_settings.R\n")

########################################################
### general settings
########################################################

nCpu <- <NBR_AVAILABLE_CPU> # (we used 30)

########################################################
### set paths to main files needed accross the workflow
########################################################

historyDT_file <- <PATH_TO_ANCIENT_ENTREZ_FILE>
#entrezID        mappingID
#7003    8
#51596   42
#1261    44
#54714   45
#7291    57

# file with mapping from entrez to chromosomic positions
#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.Rdata") # previous also held symbols
# holds only netrezID and positions
entrezDT_file <- <PATH_TO_GFF_FILE>
#entrezID        chromo  start   end     assembly        strand  symbol
#100287102       chr1    11874   14409   GRCh37.p13      +       DDX11L1
#653635  chr1    14362   29370   GRCh37.p13      -       WASH7P
#100302278       chr1    30366   30503   GRCh37.p13      +       MIR1302-2

# file with symbol synonyms to entrezID
symbolDT_file <- <PATH_TO_ENTREZ_TO_SYMBOL_FILE>
#entrezID        symbol
#1       A1BG
#1       A1B
#1       ABG
#1       GAB
#1       HYST2477
#2       A2M
#2       A2MD
#2       CPAMD5
#2       FWP007

# mapping entrezID and ensemblID
ensemblDT_file <-  <PATH_TO_ENTREZ_TO_ENSEMBLE_FILE>
#entrezID        ensemblID
#1       ENSG00000121410
#2       ENSG00000175899
#3       ENSG00000256069


########################################################
### which kind of data to prepare and use across pipeline
########################################################

# if TRUE, from scripts 2_ to further, use only the genes that pass the DE analysis threshold for all analyses
useFilterCountData <- <TRUE/FALSE # (we used TRUE)

# use only TADs (and corresponding genes)>= minSize et < maxQuantile and genes belonging to those TADs
useFilterSizeData <- <TRUE/FALSE> # (we used TRUE)

# if TRUE, from scripts 2_ to further, only focus on  TAD regions
useTADonly <- <TRUE/FALSE># (we used TRUE)


########################################################
### TAD filter number of genes - for filter data (above) 
########################################################

# keep TAD with >= minNbrGeneTAD
minNbrGeneTAD <- <MIN_NBR_GENES_BY_TAD> # (we used 3)

# keep TAD with <= maxQuant Genes 
maxQuantGeneTAD <- <MAX_QUANT_KEEP> # we used 0.99>


########################################################
### for computing correlations - 4_runMeanTADCorr.R and 7_runPermutationsMeanTADCorr.R
########################################################
 
# which correlation method to use for computing correlations
corrMethod <- <"pearson"/"spearman"/"kendall"> # we used "pearson"

# should include correlation with itself (step 4 only)
withDiag <- <TRUE/FALSE> # we used FALSE


########################################################
### permutation settings - 5fc_runPermutationsMedian.R        
########################################################

# if TRUE, perform gene-to-TAD permutations inside class of expression level
withExprClass <- <TRUE/FALSE>  # (we used TRUE)
# if withExprClass, number of classes of expression used for permutate expression data
permutExprClass <- <NBR_EXPR_CLASS> #  (we used 5)

########################################################
### for fast saving in step 5fc and 6
########################################################

pigz_exec_path <- <PATH_TO_PIGZ_EXECUTABLE>

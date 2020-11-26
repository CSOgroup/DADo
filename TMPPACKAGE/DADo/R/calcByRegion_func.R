

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' RatioDown by domain
#'
#' Function that calculates the ratioDown score for each region.
#'
#' @param g2TADdt Gene-to-domain dataframe.
#' @param DEdt Differential expression table.
#' @return The ratioFC score.
#' @export

# this version uses left_joining and aggregate, return the ratio for all regions
get_downByRegion_v2 <- function(g2TADdt, DEdt) {
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- dplyr::left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  nGenesDT <- aggregate(genes ~ region, data = mergedDT, FUN = length)
  nGenes <- setNames(nGenesDT$genes, nGenesDT$region)
  ratioDownDT <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(x<0))
  ratioDown <- setNames(ratioDownDT$logFC, ratioDownDT$region)
  all_ratio <- ratioDown/nGenes[names(ratioDown)]
  stopifnot(all_ratio >= 0 & all_ratio <= 1)
  return(all_ratio)
}  


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' RatioFC by domain
#'
#' Function that calculates the ratioFC score for each region.
#'
#' @param g2TADdt Gene-to-domain dataframe.
#' @param DEdt Differential expression table.
#' @return The ratioFC score.
#' @export

# instead of getting the ratio of genes that are down-regulated, returns the ratio of FC that is down 
get_FCdownByRegion_v2 <- function(g2TADdt, DEdt) {
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- dplyr::left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  stopifnot(is.numeric(mergedDT$logFC[1]))
  # get the absolute value of tot. logFC by region
  totLogFC <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(abs(x)))
  totLogFC <- setNames(totLogFC$logFC, totLogFC$region)
  # get the absolute value of negative logFC  
  negLogFC <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(abs(x[x<0])))
  negLogFC <- setNames(negLogFC$logFC, negLogFC$region)
  stopifnot(all(names(totLogFC) == names(negLogFC)))
  all_ratio <- negLogFC/totLogFC[names(negLogFC)]
  
  # CHANGED 24.01.2018: PROBLEM IF A DOMAIN CONTAINS ALL GENES WITH 0 LOGFC (HAPPENS TOPDOM GSE71119 DEDIFF MFSM)
  # in this case ratioDown should be zero not NA
  poss_NA <- which(negLogFC==0 & totLogFC==0)
  
  if(length(poss_NA) > 0){
    stopifnot(is.na(all_ratio[poss_NA]))
    all_ratio[poss_NA] <- 0
  }
  
  stopifnot(all_ratio >= 0 & all_ratio <= 1)
  return(all_ratio)
}  


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# !!!! WARNING PIPELINE VERSION GENE NAMES ARE ROWNAMES NOT 1ST COL

#' RatioFC by domain
#'
#' Function that calculates the ratioFC score for each region.
#'
#' @param DEdt Differential expression table.
#' @param shuffData The permutation data.
#' @param stat The score to compute.
#' @param ncpu The number of CPUs to use.
#' @param TADonly If only genes from domains should be used.
#' @return The score for each permtuation.
#' @export


# in this version also I only take as reference the regions I have in the permut data  (not the regions in the gene2tadDT)
get_statFromShuffle_para  <- function(DEdt, shuffData, stat, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(stat %in% c("ratioDown", "FCdown"))
  doMC::registerDoMC(cores = ncpu)
  # ensure always have the same
  regions <- as.character(unique(shuffData[,1]))
  if(TADonly)
    regions <- regions[grep("_TAD", regions)]
  regions <- sort(regions)

  statDT <- foreach::foreach(i_perm = 1:ncol(shuffData), .combine='cbind') %dopar% {
	cat(paste0("...... ratioDown permut: ", i_perm, "/", ncol(shuffData), "\n"))
    shuff_g2TAD <- data.frame(entrezID = rownames(shuffData), region = shuffData[,i_perm], stringsAsFactors =F)
    shuff_g2TAD$entrezID <- as.character(shuff_g2TAD$entrezID)
    shuff_g2TAD$region <- as.character(shuff_g2TAD$region)
    if (stat == "ratioDown") {
      x <- get_downByRegion_v2(shuff_g2TAD, DEdt) 
    } else if(stat == "FCdown") { 
      x <- get_FCdownByRegion_v2(shuff_g2TAD, DEdt)     
    } else  {
      stop("should not happen\n")
    }
    x[regions]
  }
  stopifnot(all(rownames(statDT) == regions))
  return(statDT)
}




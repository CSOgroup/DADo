################################################################# !!!!!!!!!!!!!!!!!!! warning hard-coded
# - genome assembly used with biomart: GRCh37
##########################

# FLUX: (instead of installing flux package) if not working
#auc <- function (x, y, thresh = NULL, dens = 100, sort.x = TRUE) {
#    x <- x[!is.na(x)]
#    y <- y[!is.na(x)]
#    x <- x[!is.na(y)]
#    y <- y[!is.na(y)]
#    if (sort.x) {
#        ord <- order(x)
#        x <- x[ord]
#        y <- y[ord]
#    }
#    idx = 2:length(x)
#    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], 
#        x[2], length.out = dens)))
#    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], 
#        x[2], length.out = dens)))
#    if (!is.null(thresh)) {
#        y.0 <- y <= thresh
#        y[y.0] <- thresh
#    }
#    idx = 2:length(x)
#    integral <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + 
#        y[idx - 1]))/2
#    integral
#}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

quantNorm <- function(x) {qqnorm(x, plot.it=F)$x}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# this was used for microarray data
#madNorm <- function(x) {
#  rna_rnaseqDT_tmp <- x
#  # median center
#  all_meds <- apply(rna_rnaseqDT_tmp, 1, median)
#  rna_rnaseqDT_tmp <- sweep(rna_rnaseqDT_tmp, 1, all_meds, "-")  # substract from each row their corresponding median
#  # mad norm
#  # In order to use the MAD as a consistent estimator for the estimation of the standard deviation σ, one takes
#  # σ ^ = k ⋅ MAD where k is a constant scale factor, which depends on the distribution.[1]
#  # For normally distributed data k is taken to be: ~1.4826
#  all_mads <-1.4826 * apply(abs(rna_rnaseqDT_tmp), 1, median)
#  rna_madnorm_rnaseqDT <- sweep(rna_rnaseqDT_tmp, 1, all_mads, "/") 
#  stopifnot( all ( dim(x) == dim(rna_madnorm_rnaseqDT)))
#  rownames(rna_madnorm_rnaseqDT) <- rownames(x)
#  colnames(rna_madnorm_rnaseqDT) <- colnames(x)
#  return(rna_madnorm_rnaseqDT)
#}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  
get_geneList_fromEntrez <- function(refList, g2t, histDT_file) {
  # look if find other history mapping
  # replace if some of the input entrezID have a correspondance to the newest one stored in entrezDT 
  historyDT <- read.delim(paste0(histDT_file), header=TRUE, stringsAsFactors = F)
  historyDT$entrezID <- as.character(historyDT$entrezID)
  historyDT$mappingID <- as.character(historyDT$mappingID)
  refList_v2 <- unlist(sapply(refList, function(x) ifelse (x %in% historyDT$mappingID,
                                                                              historyDT$entrezID[historyDT$mappingID == x], x)))
  # the name is now refList, but the value not necessarily geneList
  if(!is.null(g2t))
    refList_v2 <- refList_v2[refList_v2 %in% g2t$entrezID]
  # remove if any duplicated
  dupID <- unique(refList_v2[which(duplicated(refList_v2))])
  refList_v2 <- refList_v2[!refList_v2 %in% dupID]
  stopifnot(all(names(refList_v2) %in% refList))
  stopifnot(!any(duplicated(refList_v2)))
  # return in correct order
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  returnList != names(returnList)
  return(returnList)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  

# refList is the list of symbol
get_geneList_fromSymbol <- function(refList, g2t, symbDT_file) {
  symbDT <- read.delim(paste0(symbDT_file), header=TRUE, stringsAsFactors = F)
  symbDT$entrezID <- as.character(symbDT$entrezID)
  symbDT <- symbDT[symbDT$symbol %in% refList,]
  if(!is.null(g2t))
    symbDT <- symbDT[symbDT$entrezID %in% g2t$entrezID, ]
  
  dubSymb <- unique(symbDT$symbol[which(duplicated(symbDT$symbol))])
  symbDT <- symbDT[!symbDT$symbol %in% dubSymb, ]
  dubEntrez <- unique(symbDT$entrezID[which(duplicated(symbDT$entrezID))])
  symbDT <- symbDT[!symbDT$entrezID %in% dubEntrez, ]
  refList_v2 <- setNames(symbDT$entrezID, symbDT$symbol)
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  return(returnList)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  

get_geneList_fromEnsembl <- function(refList, g2t, ensDT_file) {

  ensDT <- read.delim(paste0(ensDT_file), header=TRUE, stringsAsFactors = F)
  ensDT$entrezID <- as.character(ensDT$entrezID)
  ensDT <- ensDT[ensDT$ensemblID %in% refList,]
  
  if(!is.null(g2t))
    ensDT <- ensDT[ensDT$entrezID %in% g2t$entrezID, ]
  
  dubEns <- unique(ensDT$ensemblID[which(duplicated(ensDT$ensemblID))])
  ensDT <- ensDT[!ensDT$ensemblID %in% dubEns, ]
  
  dubEntrez <- unique(ensDT$entrezID[which(duplicated(ensDT$entrezID))])
  ensDT <- ensDT[!ensDT$entrezID %in% dubEntrez, ]
  refList_v2 <- setNames(ensDT$entrezID, ensDT$ensemblID)
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  return(returnList)
  

}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

printAndLog <- function(mytext, mylogf) {
  cat(mytext)
  cat(mytext, append=T, file = mylogf)
}

#######################################################################################################################
####################################################################################################################### added updates from TAD_DE_utils_fasterPermut
#######################################################################################################################

#get_multiShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, nSimu, nCpu, aggregFun) {
#  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
#  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
#  if(withExprClass) {
#    stopifnot(!is.null(RNAdt) & !is.null(nClass))
#  }
#  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
#  g2TADdt$region <- as.character(g2TADdt$region)
#  registerDoMC(nCpu)
#  # need to ensure that I get the same order for the genes
#  # do the first one
#  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
#                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun) 
#  colnames(allT) <- c(colnames(allT)[1], paste0(colnames(allT)[2], "1")) # region1
#  genes1 <- allT$entrezID
#  if(nSimu >1){
#    tmpDT <- foreach(i=2:nSimu, .combine='cbind') %dopar% {
#      if(withExprClass) {
#        cat(paste0("... WITH CLASS ", aggregFun, " - shuffle: ", i, "/", nSimu, "\n"))
#      } else{
#        cat(paste0("... NO CLASS - shuffle: ", i, "/", nSimu, "\n"))
#      }
#      x <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
#                                       nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun) 
#      stopifnot(all(genes1 == x[,1]))
#      x[,2]
#    }
#    colnames(tmpDT) <- paste0("region", 2:nSimu)
#    allT <- cbind(allT, tmpDT)
#  }
#  return(allT)  
#}  



# LAST UPDATE 16.08.2019 => set.seed in get_ShuffledPositions_vFunct using permut idx to make reproducible !
# add also the stopifnot to check the foreach assignment
# besides these 2 changes -> same as TAD_DE_utils_fasterPermut.R_noSeed

### same as _vJune but can pass aggregFun for aggregating the expression values

# UPDATE 16.08.2019 => geneAggregExpression is built only once in the multiShuffled function ! and passed here
# => aggreg expression, computed only once, RNAdt and aggregFun no need to be passed anymore

get_multiShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, nSimu, nCpu, aggregFun) {

  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  registerDoMC(nCpu)


# update 16.08.2019: this was done in the child function before
  warning("geneIDlist argument should correspond to rownames of RNAdt")
  warning("duplicated - ambiguous - are removed !!! ")
  duplicatedID <- geneIDlist[duplicated(geneIDlist)]
  RNAdt <- RNAdt[which(! geneIDlist %in% duplicatedID),]
  geneIDlist <- geneIDlist[which(! geneIDlist %in% duplicatedID)]
  rownames(RNAdt) <- geneIDlist
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
    stopifnot(!is.null(aggregFun))
  }

  # take only the genes for which we have their positions
  # and subset the rnaseq data for these genes only
  if(TADonly) {
    # take only the genes that are in TAD
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
    RNAdt <- RNAdt[geneListTAD,]
  } else {
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
    RNAdt <- RNAdt[geneListTAD,]
  }






if(withExprClass) {

# update 16.08.2019: this was done in the child function before
    # define expression classes based on median expression
    geneAggregExpression <- data.frame(gene = geneListTAD, expValue = apply(RNAdt, 1, aggregFun))
    rownames(geneAggregExpression) <- NULL
    # rank by expression value
    geneAggregExpression <- geneAggregExpression[order(geneAggregExpression$expValue),]
    # split into 'nClass' groups 
    nR <- nrow(geneAggregExpression)
    nGeneByClass <- rep(nR %/% nClass, nClass)
    # add the remainder (1 to each, to avoid having a lot more in one class)
    if(nR%%nClass > 0)
      nGeneByClass[1:(nR%%nClass)] <- nR %/% nClass + 1
    stopifnot(sum(nGeneByClass) == nR)
    geneClass <- rep(c(1:nClass), nGeneByClass)
    stopifnot(length(geneClass) == nR)
    # add a column to DF with their corresponding class
    geneAggregExpression$class <- geneClass
    #save(geneAggregExpression, file="geneAggregExpression.Rdata")
    stopifnot(length(unique(geneAggregExpression$class)) == nClass)
    # add a column with their initial TAD
    geneAggregExpression$initRegion <- sapply(geneAggregExpression$gene, function(x) g2TADdt$region[g2TADdt$entrezID==x])


} else {

geneAggregExpression <- NULL
}


# save(geneAggregExpression, file="geneAggregExpression.Rdata", version=2) # used to compare with slower version


  # need to ensure that I get the same order for the genes
  # do the first one
  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, geneIDlist = geneIDlist,  #  aggregFun=aggregFun, RNAdt = RNAdt, 
                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, geneAggregExpressionDT = geneAggregExpression )  
  colnames(allT) <- c(colnames(allT)[1], paste0("permut", "1")) # permut1
  genes1 <- allT$entrezID
  if(nSimu >1){
    tmpDT <- foreach(i=2:nSimu, .combine='cbind') %dopar% {
      if(withExprClass) {
        cat(paste0("... WITH CLASS ", aggregFun, " - shuffle: ", i, "/", nSimu, "\n"))
      } else{
        cat(paste0("... NO CLASS - shuffle: ", i, "/", nSimu, "\n"))
      }
    x <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt,  geneIDlist = geneIDlist,  #  aggregFun=aggregFun, RNAdt = RNAdt, 
                                         nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, geneAggregExpressionDT = geneAggregExpression, rd_idx=i ) 


      stopifnot(all(genes1 == x[,1]))
      x[,2]
    }
    colnames(tmpDT) <- paste0("permut", 2:nSimu)
    allT <- cbind(allT, tmpDT)
  }
  return(allT)  
}  
#######################################################################################################################
####################################################################################################################### added updates from TAD_DE_utils_fasterPermut
#######################################################################################################################

#### same as _vJune but can pass aggregFun for aggregating the expression values

#get_ShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, TADonly, withExprClass, aggregFun) {
#  warning("geneIDlist argument should correspond to rownames of RNAdt")
#  warning("duplicated - ambiguous - are removed !!! ")
#  duplicatedID <- geneIDlist[duplicated(geneIDlist)]
#  RNAdt <- RNAdt[which(! geneIDlist %in% duplicatedID),]
#  geneIDlist <- geneIDlist[which(! geneIDlist %in% duplicatedID)]
#  rownames(RNAdt) <- geneIDlist
#  if(withExprClass) {
#    stopifnot(!is.null(RNAdt) & !is.null(nClass))
#    stopifnot(!is.null(aggregFun))
#  }
#  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
#  g2TADdt$region <- as.character(g2TADdt$region)
#  # take only the genes for which we have their positions
#  # and subset the rnaseq data for these genes only
#  if(TADonly) {
#    # take only the genes that are in TAD
#    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
#    RNAdt <- RNAdt[geneListTAD,]
#  } else {
#    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
#    RNAdt <- RNAdt[geneListTAD,]
#  }
#  ##########
#  ##### DO IT BY SHUFFLING THE LABELS BY CLASS OF EXPRESSION
#  ##########
#  if(withExprClass) {
#    # define expression classes based on median expression
#    geneAggregExpression <- data.frame(gene = geneListTAD, expValue = apply(RNAdt, 1, aggregFun))
#    rownames(geneAggregExpression) <- NULL
#    # rank by expression value
#    geneAggregExpression <- geneAggregExpression[order(geneAggregExpression$expValue),]
#    # split into 'nClass' groups 
#    nR <- nrow(geneAggregExpression)
#    nGeneByClass <- rep(nR %/% nClass, nClass)
#    # add the remainder (1 to each, to avoid having a lot more in one class)
#    if(nR%%nClass > 0)
#      nGeneByClass[1:(nR%%nClass)] <- nR %/% nClass + 1
#    stopifnot(sum(nGeneByClass) == nR)
#    geneClass <- rep(c(1:nClass), nGeneByClass)
#    stopifnot(length(geneClass) == nR)
#    # add a column to DF with their corresponding class
#    geneAggregExpression$class <- geneClass
#    #save(geneAggregExpression, file="geneAggregExpression.Rdata")
#    stopifnot(length(unique(geneAggregExpression$class)) == nClass)
#    # add a column with their initial TAD
#    geneAggregExpression$initRegion <- sapply(geneAggregExpression$gene, function(x) g2TADdt$region[g2TADdt$entrezID==x])
#    # now, for each class, reshuffle the TAD -> new column with the reshuffled positions
#    geneAggregExpression$shuffRegion <- foreach(i_cl = 1:nClass, .combine='c') %do% {
#      subDT <- geneAggregExpression[geneAggregExpression$class == i_cl,]
#      initPos <- subDT$initRegion
#      newPos <- sample(initPos, size=length(initPos), replace=F)  
#      stopifnot(all(unique(as.character(initPos)) %in% unique(as.character(newPos))))
#      newPos
#    }
#    shuffGenePosDT <- data.frame(entrezID = geneAggregExpression$gene, 
#                                 region = as.character(geneAggregExpression$shuffRegion),
#                                 stringsAsFactors = F)
#  } else {
#    ##########
#    ##### JUST SHUFFLE THE LABELS, IRRESPECTIVE OF GENE EXPRESSION
#    ##########
#    newPos <- sample(g2TADdt$region[g2TADdt$entrezID %in% geneListTAD], replace=F)
#    shuffGenePosDT <- data.frame(entrezID = geneListTAD, 
#                                 region = as.character(newPos),
#                                 stringsAsFactors = F)
#  }
#  # to be compatible with the functions that use gene2tadDT, return a similar DF
#  return(shuffGenePosDT)
#}
# LAST UPDATE 16.08.2019 => set.seed in get_ShuffledPositions_vFunct using permut idx to make reproducible !
# add also the stopifnot to check the foreach assignment
# besides these 2 changes -> same as TAD_DE_utils_fasterPermut.R_noSeed

### same as _vJune but can pass aggregFun for aggregating the expression values

# UPDATE 16.08.2019 => geneAggregExpression is built only once in the multiShuffled function ! and passed here
# => aggreg expression, computed only once, RNAdt and aggregFun no need to be passed anymore

get_ShuffledPositions_vFunct <- function(g2TADdt, geneIDlist, nClass, TADonly, withExprClass, geneAggregExpressionDT=NULL, rd_idx=0) { # removed aggregFun and rnaDT
set.seed(16082019+rd_idx) # added for reproducibility

  warning("geneIDlist argument should correspond to rownames of RNAdt")
  warning("duplicated - ambiguous - are removed !!! ")
  duplicatedID <- geneIDlist[duplicated(geneIDlist)]

  if(withExprClass) {
    stopifnot( !is.null(nClass))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  # take only the genes for which we have their positions
  # and subset the rnaseq data for these genes only
  if(TADonly) {
    # take only the genes that are in TAD
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
  } else {
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
  }
  ##########
  ##### DO IT BY SHUFFLING THE LABELS BY CLASS OF EXPRESSION
  ##########
  if(withExprClass) {

    stopifnot(!is.null(geneAggregExpressionDT))     # 16.08.2019 now computed once and passed from parent function
    geneAggregExpression <- geneAggregExpressionDT

    # now, for each class, reshuffle the TAD -> new column with the reshuffled positions

stopifnot(unique(as.character(geneAggregExpression$class)) == paste0(1:nClass)) # => this is why the foreach assignment works !

    geneAggregExpression$shuffRegion <- foreach(i_cl = 1:nClass, .combine='c') %do% {
      subDT <- geneAggregExpression[geneAggregExpression$class == i_cl,]
      initPos <- subDT$initRegion
      newPos <- sample(initPos, size=length(initPos), replace=F)  
      stopifnot(all(unique(as.character(initPos)) %in% unique(as.character(newPos))))
      newPos
    }
    shuffGenePosDT <- data.frame(entrezID = geneAggregExpression$gene, 
                                 region = as.character(geneAggregExpression$shuffRegion),
                                 stringsAsFactors = F)
  } else {
    ##########
    ##### JUST SHUFFLE THE LABELS, IRRESPECTIVE OF GENE EXPRESSION
    ##########
    newPos <- sample(g2TADdt$region[g2TADdt$entrezID %in% geneListTAD], replace=F)
    shuffGenePosDT <- data.frame(entrezID = geneListTAD, 
                                 region = as.character(newPos),
                                 stringsAsFactors = F)
  }
  # to be compatible with the functions that use gene2tadDT, return a similar DF
  return(shuffGenePosDT)
}



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# this version uses left_joining and aggregate, return the ratio for all regions
get_downByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
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
# instead of getting the ratio of genes that are down-regulated, returns the ratio of FC that is down 
get_FCdownByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
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
  
  # CHANGED 24.01.2018: PROBLEM IF A TAD CONTAINS ALL GENES WITH 0 LOGFC (HAPPENS TOPDOM GSE71119 DEDIFF MFSM)
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

# in this version also I only take as reference the regions I have in the permut data  (not the regions in the gene2tadDT)
get_statFromShuffle_para  <- function(DEdt, shuffData, stat, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(stat %in% c("ratioDown", "FCdown"))
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  registerDoMC(cores = ncpu)
  # ensure always have the same
  regions <- as.character(unique(shuffData[,1]))
  if(TADonly)
    regions <- regions[grep("_TAD", regions)]
  regions <- sort(regions)

  statDT <- foreach(i_perm = 1:ncol(shuffData), .combine='cbind') %dopar% {
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


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

stouffer <- function(ps, two.tails=FALSE) {
	if(two.tails) ps = ps / 2
	# transform p-values into quantiles of standard normal
	qps = qnorm(1-ps, lower.tail = TRUE)
	# take the average of the quantiles
	iqps = sum(qps) / sqrt(length(qps))
	# find p-val of the average
	p = 1 - pnorm(iqps, lower.tail = TRUE)
	if(two.tails) p = 2 * p
	return(p)
}

#the individual p-values are transformed into the quantiles of a standard normal
#the p-val of the average of the quantiles is then found

# Fisher's method treates large and small p-values asymmetrically
# e.g. combining 0.999 and 0.001 gives 0.008
# is asymmetrically sensitive to small p-values

# Z-transform test: one-to-one mapping of the standard normal curve of the p-value
# Z = standard deviate (= a number drawn from from normal distribution with mean 0 and sd 1)
# the test converts p-values into standard normal deviates
# sum of the Z divided by square root of the number of tests has a std normal distribution
# no asymmetry problem

# weighted version of the Z-transform test: a weight can be assigned to each test
# if each test is given equal weight, this reduces to Z-transform test
# ideally, each study is weighted proportional to the inverse of its error variance
# (i.e. by the reciprocal of its squared standard error)
# more generally, the weights should be the inverse of the squared standard error fo the effect
# size estimate for each study

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

get_fcc <- function(fc_vect) {
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

get_ratioDown <- function(fc_vect) {
  sum(fc_vect < 0)/length(fc_vect) 
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# function code taken from fastSave package !

my_save.pigz <- function (..., list = character(), file = stop("'file' must be specified"), 
          pigz_exec_path = "",
          envir = parent.frame(), n.cores = 4, eval.promises = TRUE, 
          precheck = TRUE) 
{
  stopifnot(file.exists(pigz_exec_path))
  if (.Platform$OS.type == "unix") {
    if (system(paste0("command -v ", pigz_exec_path), wait = T, ignore.stdout = TRUE, 
               ignore.stderr = TRUE) > 0) {
      stop("The pigz command is not available on this system!")
    }
  }
  else {
    stop("Platform is not a unix system!")
  }
  if (!is.numeric(n.cores)) 
    stop("'n.cores' mut be numeric")
  names <- as.character(substitute(list(...)))[-1L]
  if (missing(list) && !length(names)) 
    warning("nothing specified to be save()d")
  list <- c(list, names)
  if (precheck) {
    ok <- vapply(list, exists, NA, envir = envir)
    if (!all(ok)) {
      n <- sum(!ok)
      stop(sprintf(ngettext(n, "object %s not found", "objects %s not found"), 
                   paste(sQuote(list[!ok]), collapse = ", ")), domain = NA)
    }
  }
  on.exit(close(con))
  n.cores.arg <- paste0(" --processes ", n.cores)
  con <- pipe(paste0(pigz_exec_path, " ", n.cores.arg, " > ", file))
  save(..., list = list, file = con, ascii = FALSE, version = 2, 
       envir = envir, eval.promises = eval.promises, precheck = precheck)
}
# x = 2
# save.pigz(x, file="x.Rdata")
# myx=5
# my_save.pigz(myx, pigz_exec_path = pigz_exec_path, file="myx.Rdata")
# x <- get(load("x.Rdata"))
# x
# myx <- get(load("myx.Rdata"))
# myx

#############################################################################################################################
#############################################################################################################################  # added from TAD_DE_utils_meanCorr
#############################################################################################################################

get_meanCorr_value <- function(exprMatrix, inside_genes, outside_genes, cormet) {
  stopifnot(inside_genes %in% rownames(exprMatrix))
  stopifnot(outside_genes %in% rownames(exprMatrix))
  stopifnot(setequal(c(inside_genes, outside_genes), rownames(exprMatrix)))
  
  nAllGenes <- length(inside_genes) + length(outside_genes)
  
  coexprMatrix <- cor(t(exprMatrix), method = cormet)
  stopifnot(dim(coexprMatrix) == nAllGenes)
  
  coexprMatrix[lower.tri(coexprMatrix, diag = TRUE)] <- NA   # because after I filter that 1 gene should be inside, and 1 gene should be outside -> can never happen to take the diag. value of coexpression
  coexprMatrix <- na.omit(melt(coexprMatrix))
  colnames(coexprMatrix)[1:2] <- c("Var1", "Var2")
  stopifnot(colnames( coexprMatrix)[3] == "value" )
  coexprMatrix$Var1 <- as.character(coexprMatrix$Var1)
  coexprMatrix$Var2 <- as.character(coexprMatrix$Var2)
  
  stopifnot(coexprMatrix$Var1 %in% outside_genes | coexprMatrix$Var1 %in% inside_genes)
  stopifnot(coexprMatrix$Var2 %in% outside_genes | coexprMatrix$Var2 %in% inside_genes)
  stopifnot(inside_genes %in% coexprMatrix$Var1 | inside_genes %in% coexprMatrix$Var2)
  stopifnot(outside_genes %in% coexprMatrix$Var1 | outside_genes %in% coexprMatrix$Var2)
  
  # take only if one of the two genes outside and the other inside
  coexprMatrix <- coexprMatrix[!  (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% outside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  coexprMatrix <- coexprMatrix[ ! (coexprMatrix$Var1 %in% inside_genes & coexprMatrix$Var2 %in% inside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  
  stopifnot(     (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% inside_genes) | 
                   (coexprMatrix$Var2 %in% outside_genes & coexprMatrix$Var1 %in% inside_genes) )
  
  

  stopifnot(nrow(coexprMatrix) == (length(inside_genes) * length(outside_genes) ))
  
  meanCorr_value <- mean(coexprMatrix$value)
  stopifnot(!is.na(meanCorr_value))
  return(meanCorr_value)
}



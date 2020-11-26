

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# LAST UPDATE 16.08.2019 => set.seed in get_ShuffledPositions_vFunct using permut idx to make reproducible !
# add also the stopifnot to check the foreach assignment
# besides these 2 changes -> same as TAD_DE_utils_fasterPermut.R_noSeed

### same as _vJune but can pass aggregFun for aggregating the expression values

# UPDATE 16.08.2019 => geneAggregExpression is built only once in the multiShuffled function ! and passed here
# => aggreg expression, computed only once, RNAdt and aggregFun no need to be passed anymore


#' Gene-to-domain permutations
#' 
#' Performs the permutation (parallelized).
#'
#' @param g2TADdt Gene-to-domain dataframe (expected columns: entrezID, region).
#' @param RNAdt Gene expression dataframe (rownames should correspond to entrezID, samples in columns).
#' @param geneIDlist List of gene IDs that should be used.
#' @param nClass The number of classes of expression in which genes are shuffled.
#' @param withExprClass If shuffling should take place within classes.
#' @param TADonly If only genes from domains should be used.
#' @param nSimu The number of permutations.
#' @param nCpu The number of CPUs to use.
#' @param aggregFun With which function gene expression should be aggregated (across samples).
#' @return The permtuation data.
#' @export


get_gene2tad_multiPermut <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, nSimu, nCpu, aggregFun) {
  	
  if(!suppressPackageStartupMessages(require("foreach"))) stop("-- foreach package required\n")  
  if(!suppressPackageStartupMessages(require("doMC"))) stop("-- doMC package required\n")  

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
    # take only the genes that are in domain
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
    # add a column with their initial domain
    geneAggregExpression$initRegion <- sapply(geneAggregExpression$gene, function(x) g2TADdt$region[g2TADdt$entrezID==x])


  } else {

    geneAggregExpression <- NULL
  }


# save(geneAggregExpression, file="geneAggregExpression.Rdata", version=2) # used to compare with slower version


  # need to ensure that I get the same order for the genes
  # do the first one
  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, geneIDlist = geneIDlist,  #  aggregFun=aggregFun, RNAdt = RNAdt, 
                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, geneAggregExpressionDT = geneAggregExpression )  
  colnames(allT) <- c(colnames(allT)[1], paste0("permut", "1")) # region1
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
#######################################################################################################################
#######################################################################################################################
#' Gene-to-domain permutations
#' 
#' Performs the permutation (for 1 permutation).
#'
#' @param g2TADdt Gene-to-domain dataframe (expected columns: entrezID, region).
#' @param geneIDlist List of gene IDs that should be used.
#' @param nClass The number of classes of expression in which genes are shuffled.
#' @param TADonly If only genes from domains should be used.
#' @param withExprClass If shuffling should take place within classes.
#' @param geneAggregExpressionDT Dataframe with aggregated gene expression (across samples).
#' @param rd_idx For setting seed.
#' @return The permtuation data.
#' @export

### same as _vJune but can pass aggregFun for aggregating the expression values

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
    # take only the genes that are in domain
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

    # now, for each class, reshuffle the domain -> new column with the reshuffled positions

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
#' Sample genes across domain boundaries
#' 
#' Sample genes across domain boundaries (as many as there are in a given domain).
#'
#' @param g2t_DT Gene-to-domain dataframe (expected columns: entrezID, region, end, start).
#' @param tadpos_DT Dataframe with domain positions (expected columns: region, end, start).
#' @return A list with permutation data (either left or right/left only/right only).
#' @export


getSampleAcrossBD <- function(g2t_DT, tadpos_DT){
  
  stopifnot(c("start", "end", "region") %in% colnames(g2t_DT))
  stopifnot(c("start", "end", "region") %in% colnames(tadpos_DT))

  all_genes <- as.character(g2t_DT$entrezID)
  
  g2t_DT$region <- as.character(g2t_DT$region)
  tadpos_DT$region <- as.character(tadpos_DT$region)
  tadpos_DT$chromo <- as.character(tadpos_DT$chromo)
  
  g2t_DT$end <- as.numeric(g2t_DT$end)
  stopifnot(!is.na(g2t_DT$end))
  g2t_DT$start <- as.numeric(g2t_DT$start)
  stopifnot(!is.na(g2t_DT$start))
  tadpos_DT$end <- as.numeric(tadpos_DT$end)
  stopifnot(!is.na(tadpos_DT$end))
  tadpos_DT$start <- as.numeric(tadpos_DT$start)
  stopifnot(!is.na(tadpos_DT$start))
  
  tadpos_DT$mid_pos <- (tadpos_DT$start+tadpos_DT$end)/2
  g2t_DT$mid_pos <- (g2t_DT$start+g2t_DT$end)/2
  
  stopifnot(g2t_DT$region %in% tadpos_DT$region)
  
  all_tads <- sort(unique(as.character(g2t_DT$region)))
  
  
  
  sample_around_TADs <- foreach(reg = all_tads) %dopar% {
    cat("...... start domain : \t", reg, "\n")
    
    curr_chromo <- as.character(tadpos_DT$chromo[tadpos_DT$region == reg])
    
    curr_start <- (tadpos_DT$start[tadpos_DT$region == reg])
    stopifnot(is.numeric(curr_start))
    
    curr_end <- (tadpos_DT$end[tadpos_DT$region == reg])
    stopifnot(is.numeric(curr_end))
    
    stopifnot(length(curr_chromo) == 1)
    stopifnot(length(curr_start) == 1)
    stopifnot(length(curr_end) == 1)
    
    curr_midPos <- (curr_start+curr_end)/2
    stopifnot(curr_midPos == tadpos_DT$mid_pos[tadpos_DT$region == reg])
    
    reg_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
    stopifnot(length(reg_genes) > 0)
    
    curr_nGenes <- length(reg_genes)
    
    # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
    # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
    
    curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,,drop=FALSE]
    stopifnot(nrow(curr_g2t) > 0)
    
    stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
    curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),,drop=FALSE]
    
    # distance to domain center
    curr_genesOutsideDT <- curr_g2t[ !curr_g2t$entrezID %in% reg_genes,,drop=FALSE]
    stopifnot(nrow(curr_genesOutsideDT) > 0)
    
    #>>> take the same number of genes, either on left or right
    curr_genesOutsideDT$distToTAD <- abs(curr_genesOutsideDT$mid_pos - curr_midPos)
    curr_genesOutsideDT <- curr_genesOutsideDT[order(curr_genesOutsideDT$distToTAD, decreasing=FALSE),,drop=FALSE]
    
    sample_around_genes <- curr_genesOutsideDT$entrezID[1:curr_nGenes]
    
    all_dist <- curr_genesOutsideDT$distToTAD[1:curr_nGenes]
    stopifnot(!is.na(all_dist))
    
    stopifnot(!is.na(curr_genesOutsideDT))
    
    stopifnot(length(all_dist) == length(sample_around_genes))
    
    #>>> take the same number of genes on the left
    curr_genesOutsideDT_left <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos < curr_midPos,]  # SMALLER
    curr_nGenes_left <- min(curr_nGenes, nrow(curr_genesOutsideDT_left)) 
    
    if(curr_nGenes_left > 0) {
      sample_around_genes_left <- curr_genesOutsideDT_left$entrezID[1:curr_nGenes_left]
      all_dist_left <- curr_genesOutsideDT_left$distToTAD[1:curr_nGenes_left]  
      stopifnot(!is.na(all_dist_left))
      ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
      stopifnot(sample_around_genes_left %in% all_genes)
      stopifnot(!sample_around_genes_left %in% reg_genes)
      
    } else {
      sample_around_genes_left <- character(0)
      all_dist_left <- c()
    }
    
    curr_genesOutsideDT_right <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos >= curr_midPos,]  # BIGGER OR EQUAL
    curr_nGenes_right <- min(curr_nGenes, nrow(curr_genesOutsideDT_right)) 
    
    if(curr_nGenes_right > 0) {
      sample_around_genes_right <- curr_genesOutsideDT_right$entrezID[1:curr_nGenes_right]
      all_dist_right <- curr_genesOutsideDT_right$distToTAD[1:curr_nGenes_right]
      stopifnot(!is.na(all_dist_right))
      
      stopifnot(sample_around_genes_right %in% all_genes)
      stopifnot(!sample_around_genes_right %in% reg_genes)
      
    } else {
      sample_around_genes_right <- character(0)
      all_dist_right <- c()
      
    }
    
    stopifnot(length(sample_around_genes) == curr_nGenes)
    stopifnot(length(sample_around_genes_left) <= curr_nGenes)
    stopifnot(length(sample_around_genes_right) <= curr_nGenes)
    
    stopifnot(length(all_dist) == length(sample_around_genes))
    stopifnot(length(all_dist_right) == length(sample_around_genes_right))
    stopifnot(length(all_dist_left) == length(sample_around_genes_left))
    
    ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
    stopifnot(sample_around_genes %in% all_genes)
	stopifnot(! sample_around_genes_left %in% sample_around_genes_right)
	stopifnot( sample_around_genes %in% sample_around_genes_left | sample_around_genes %in% sample_around_genes_right)
    
    list(
      tad_genes = reg_genes,
      
      genes = sample_around_genes,
      nGenes = length(sample_around_genes),
      minDist = min(all_dist),
      maxDist = max(all_dist),
      
      genes_left = sample_around_genes_left,
      nGenes_left = curr_nGenes_left,
      minDist_left = min(all_dist_left),
      maxDist_left = max(all_dist_left),
      
      genes_right = sample_around_genes_right,
      nGenes_right =  curr_nGenes_right,
      minDist_right = min(all_dist_right),
      maxDist_right = max(all_dist_right)
    )
  } # end foreach-iterating over domains
  names(sample_around_TADs) <- all_tads
  return(sample_around_TADs)
}





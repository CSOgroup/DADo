#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# FLUX: (instead of installing flux package)
#' AUC calculation
#' 
#' Taken from "flux" package (not able to install package after update of R). Compute the Area Under the Curve (AUC) between two vectors of points.
#'
#' @param x Vector of x positions.
#' @param y Vector of y positions.
#' @param thresh Threshold below which area is not calculated.
#' @param dens By default the data density is artificially increased by adding 100 data points between given adjacent data points.
#' @param By default the vectors in x and y are ordered along increasing x because integration makes no sense with unordered data.
#' @return The AUC.
#' @export
auc <- function (x, y, thresh = NULL, dens = 100, sort.x = TRUE) {
    x <- x[!is.na(x)]
    y <- y[!is.na(x)]
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    if (sort.x) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
    }
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    if (!is.null(thresh)) {
        y.0 <- y <= thresh
        y[y.0] <- thresh
    }
    idx = 2:length(x)
    integral <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + 
        y[idx - 1]))/2
    integral
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' Stouffer's method for p-values combination
#'
#' Function that combines p-values using Stouffer's method.
#'
#' @param ps The p-values to be combined.
#' @return The normalized dataframe.
#' @export

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

#' Fold-change concordance (FCC) score
#'
#' Function that calculates the FCC score for a vector of fold-change (FC) values.
#'
#' @param fc_vect Vector of fold-changes.
#' @return The FCC score.
#' @export

get_fcc <- function(fc_vect) {
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' RatioDown score
#'
#' Function that calculates ratio of down-regulated genes for a vector of fold-change (FC) values.
#'
#' @param fc_vect Vector of fold-changes.
#' @return The ratioDown score.
#' @export
get_ratioDown <- function(fc_vect) {
  sum(fc_vect < 0)/length(fc_vect) 
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#' RatioFC score
#'
#' Function that calculates ratio of negative fold-change (FC) for a vector of FC values.
#'
#' @param fc_vect Vector of fold-changes.
#' @return The ratioFC score.
#' @export
get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#' AUC ratio 
#'
#' Function that calculates the AUC ratio for FCC cumsum curves (observed/permutation).
#'
#' @param fcc_vect Vector of observed FCCs.
#' @param fcc_permDT Dataframe with permutation FCCs (each column is a permutation).
#' @param permQt Quantile of the permutation (default: 0.95) used for computing the AUC ratio.
#' @param doPlot If true, plot the cumsum wave curves.
#' @param plotCex Cex value for the axis labs and titles (default: 1.2).
#' @param ploygonPermutCol Color for the area of permutation value (default: grey color).
#' @param qt95PermutCol Color for the line for the permQt permutation values (default: darker grey).
#' @param lwdObs Width of the line of the observed data (default: 1.2).
#' @param colObs Color for the observed data (default: "darkred").
#' @param ... Other arguments are passed to plot() function.
#' @return (invisible) A list containing the observed AUC (observed_auc), the AUC for the permutation (permut_auc) and the AUC ratio (auc_ratio).
#' @export

get_auc_ratio <- function(fcc_vect, fcc_permDT, permQt=0.95, doPlot=FALSE,
  							plotCex=1.2, polygonPermutCol=rgb(188/255,188/255,188/255, 0.3), qt95PermutCol=rgb(140/255,140/255,140/255),
								lwdObs=1.2, colObs="darkred", ...) {
  pointObsCol <- colObs
  # sort decreasing the FCC
  obs_fcc <- sort(fcc_vect, decreasing = TRUE)
  permut_FCC_unsort <- fcc_permDT
  permut_FCC <- apply(permut_FCC_unsort, 2, sort, decreasing=TRUE)
  rownames(permut_FCC) <- NULL
  cat(paste0("... found obs. domains:\t", length(obs_fcc), "\n"))
  cat(paste0("... found permut. domains:\t", nrow(permut_FCC), "\n"))
  maxTADs <- max(c(length(obs_fcc), nrow(permut_FCC)))
  
  maxRankPlot <- ceiling(maxTADs/1000)*1000
  x_val <- c(1:length(obs_fcc))
  cumsum_permut_dt <- apply(permut_FCC, 2, function(x) cumsum(x))
  qt95Permut_cumsum <- apply(cumsum_permut_dt, 1, quantile, probs=permQt) # for each rank, take the quantile
  
  obs_cumsum <- cumsum(obs_fcc)
  
  auc_obs <- auc(x = x_val, y = obs_cumsum)
  auc_permutQt <- auc(x = x_val, y = qt95Permut_cumsum)
  auc_ratioQt <- auc_obs/auc_permutQt
  
  
  fcc_auc_ratios <- list(
    observed_auc=auc_obs,
    permut_auc=auc_permutQt,
    auc_ratio = auc_ratioQt
  )
  
  if(doPlot){
    pct_inc_qt <- round(auc_obs/auc_permutQt,2)
    my_xlab <- paste0("Domains ranked by FCC")
	my_ylab <- paste0("FCC cumulative sum")
    par(bty="l")
    plot(obs_cumsum ~ x_val,
         type = "l",
         pch = 16, cex = 0.7,
         xlab= my_xlab, 
         ylab= my_ylab,
         cex.main = plotCex,
         cex.lab = plotCex,
         cex.axis = plotCex,
         col = colObs,
         axes=FALSE,
         lwd = lwdObs,
         bty="l", ...)
    box()
    axis(2, cex.axis=plotCex, cex.lab=plotCex)
    axis(1, cex.axis=plotCex, cex.lab=plotCex, at = seq(from=0, to=maxRankPlot, by=maxRankPlot/10)) # we used 2000 and 200 hard-coded
    polygon(x = c(x_val, rev(x_val)), 
            y = c( apply(cumsum_permut_dt, 1, function(x) min(x)), rev(apply(cumsum_permut_dt, 1, function(x) max(x)))),
            border=NA,
            col = polygonPermutCol)
    lines(
      x = x_val,
      y= qt95Permut_cumsum,
      col = qt95PermutCol
    )
    legend("topleft",
           xjust=0.5, yjust=0,
           pch = c(16, 15, -1 ),
           lty=c(-1,-1,1),
           legend = c(paste0("observed (n=", length(obs_fcc), ")"), "min-max permut.", paste0(permQt, "-qt permut.")), 
           pt.cex = c(0.7, 2),
           col = c(pointObsCol, polygonPermutCol, qt95PermutCol),
           bty="n")
    legtxt <- as.expression(bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc_qt)))
    legend("right", legend=c(legtxt), bty="n")
  }
  invisible(fcc_auc_ratios)
}




#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Mean intra-domain correlation
#'
#' Function the average pairwise correlations within domains.
#'
#' @param gene2tad_dt Dataframe with gene-to-domain assignment (required columns: "entrezID" for gene IDs, "region" for domain IDs).
#' @param expr_dt Dataframe with gene expression values for which to compute correlation (row names should be gene IDs; all gene IDs from gene2tad_dt should be available).
#' @param corrMeth Correlation method (one of "pearson", "kendall"or "spearman"; default:"pearson").
#' @param minNbrGenes Compute correlation for a domain only if >= minNbrGenes belong to it.
#' @param withDiag If correlation with itself should be included (default: FALSE).
#' @param nCpu Number available CPU.
#' @return A 2-column dataframe (region/meanCorr colums).
#' @export

get_meanCorr <- function(gene2tad_dt, exprd_dt, corrMeth="pearson",  minNbrGenes=3, withDiag=FALSE, nCpu=2){
  corrMeth <- match.arg(corrMeth, choices=c("pearson", "kendall", "spearman"))
  if(!suppressPackageStartupMessages(require("foreach"))) stop("-- foreach package required\n")  
  if(!suppressPackageStartupMessages(require("doMC"))) stop("-- doMC package required\n")  
  registerDoMC(nCpu)
  stopifnot(gene2tad_dt$entrezID %in% rownames(exprd_dt))
  tad_meanCorr_dt <- foreach(tad=unique(as.character(gene2tad_dt$region)), .combine="rbind") %dopar% {
    tad_dt <- gene2tad_dt[as.character(gene2tad_dt$region) == tad, ]
	if(nrow(tad_dt) < 3){
	  warning("! found TAD with less than ", minNbrGenes, " genes !\n")
	  return(NULL)
	} 
    sub_qq <- exprd_dt[paste0(tad_dt$entrezID),]
    corrDT <- cor(t(sub_qq), method = corrMeth)
    stopifnot(dim(corrDT) == nrow(tad_dt))
    stopifnot(!is.na(corrDT))
    data.frame(
      region=tad,
      meanCorr=mean(corrDT[lower.tri(corrDT, diag = withDiag)], na.rm=TRUE),
      stringsAsFactors = FALSE
    )
  }
  return(tad_meanCorr_dt)
}





#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Correlation between gene expression and purity by domain
#'
#' Function the correlation between gene expression and purity for all domains.
#'
#' @param exprTable Dataframe with gene expression values for which to compute correlation (row names should be gene IDs; all gene IDs from g2tTable should be available).
#' @param purityTable Dataframe with sample purity values (required columns: purityCol and sampleCol).
#' @param g2tTable Dataframe with gene-to-domain assignment (required columns: "entrezID" for gene IDs, "region" for domain IDs).
#' @param all_samples Vector of samples. Should match exprTable columns and sampleCol column of purityTable. 
#' @param purityCol Name of the column of purityTable that holds purity values.
#' @param sampleCol Name of the column of purityTable that holds sample IDs.
#' @param transfExpr Transformation to apply to gene expression before computing gene expression (should be a function name, such as log10).
#' @param logOffset Numeric value added to gene expression before applying transfExpr (might be useful for log-transformation).
#' @param corrMeth Correlation method (one of "pearson", "kendall"or "spearman"; default:"pearson").
#' @param nCpu Number available CPU.
#' @return A 4-column dataframe (nSampWithPurity/region/entrezID/purityCorr colums). The returned dataframe is not averaged at domain level.
#' @export

  get_meanPurityCorr <- function(exprTable, purityTable, g2tTable, all_samples, purityCol, sampleCol="Sample.ID",
                                 transfExpr="log10", logOffset = 0.01, corrMeth="pearson", nCpu=2) {
    
    if(!suppressPackageStartupMessages(require("foreach"))) stop("-- foreach package required\n")  
    if(!suppressPackageStartupMessages(require("doMC"))) stop("-- doMC package required\n")  
    registerDoMC(nCpu)
    
    corrMeth <- match.arg(corrMeth, choices=c("pearson", "kendall", "spearman"))
    
    stopifnot(c("entrezID", "region") %in% colnames(g2tTable))
    stopifnot(all_samples %in% colnames(exprTable))
    stopifnot(c(sampleCol, purityCol) %in% colnames(purityTable))
    stopifnot(g2tTable$entrezID %in% rownames(exprTable))
    
    av_samples <- all_samples[all_samples %in% purityTable[,c(sampleCol)]]
    cat("... found purity for:\t", length(av_samples), "/", length(all_samples), "\n")
    all_samples <- av_samples
    
    all_tads_purity_dt <- foreach(tad = unique(g2tTable$region), .combine='rbind') %dopar% {
      
      tad_entrez <- g2tTable$entrezID[g2tTable$region == tad]
      stopifnot(tad_entrez %in% rownames(exprTable))
      
      tad_fpkm_dt <- exprTable[tad_entrez,]
      stopifnot(nrow(tad_fpkm_dt) == length(tad_entrez))
      
      stopifnot(all_samples %in% colnames(tad_fpkm_dt))
      
      purity_values <- setNames(purityTable[purityTable[,c(sampleCol)] %in% all_samples,paste0(purityCol)],
                                purityTable[purityTable[,c(sampleCol)] %in% all_samples,c(sampleCol)])
      
      
      stopifnot(setequal(names(purity_values), all_samples))
      
      tad_dt <- data.frame(t(tad_fpkm_dt[,all_samples]), check.names=FALSE)
      
      stopifnot(is.numeric(unlist(tad_dt)))
      
      if(!is.null(transfExpr)) {
        if(grepl("log", transfExpr)) {
          tad_dt_2 <- do.call(transfExpr, list(tad_dt+logOffset))
          stopifnot(dim(tad_dt_2)==dim(tad_dt))
          stopifnot(rownames(tad_dt_2) == rownames(tad_dt))
          stopifnot(colnames(tad_dt_2) == colnames(tad_dt))
          tad_dt <- tad_dt_2
        } else {stop("unnknown\n")}
        labTransf <- transfExpr
      } else {
        labTransf <- ""
      }
      
      stopifnot(colnames(tad_dt) %in% g2tTable$entrezID)
      
      stopifnot(setequal(colnames(tad_dt),  tad_entrez))
      tad_dt$sampID <- rownames(tad_dt)
      rownames(tad_dt) <- NULL
      
      purity_subDT <- data.frame(
        sampID = names(purity_values),
        purity = purity_values,
        stringsAsFactors = FALSE
      )
      purity_expr_dt <- merge(purity_subDT, tad_dt, by="sampID", all=TRUE)
      stopifnot(ncol(purity_expr_dt) == length(tad_entrez) + 2)
      purity_expr_dt <- purity_expr_dt[,c("sampID", "purity", tad_entrez)]
      
      # correlation each column with purity 
      
      purity_expr_dt <- na.omit(purity_expr_dt)
      stopifnot(!duplicated(purity_expr_dt$sampID))
      if(nrow(purity_expr_dt) == 0) return(NULL)
      all_purityCors <- apply(purity_expr_dt[,tad_entrez],2, function(col) cor(col, purity_expr_dt$purity, method=corrMeth))
      stopifnot(!is.na(all_purityCors))
      data.frame(
        nSampWithPurity=length(purity_expr_dt$sampID),
        region = tad,
        entrezID= names(all_purityCors),
        purityCorr = as.numeric(all_purityCors),
        stringsAsFactors = FALSE
      )
    }
    return(all_tads_purity_dt)
  }
  
  
  
  

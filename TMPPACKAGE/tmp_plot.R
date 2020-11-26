
load("../data/ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanCorr_TAD.RData")
load("../data/ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanLogFC_TAD.RData")
load("../data/ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_emp_pval_combined.RData")

meanCorr=all_meanCorr_TAD
meanFC=all_meanLogFC_TAD
comb_pval=emp_pval_combined

plot_volcanoTADsCorrFC(meanCorr=all_meanCorr_TAD, 
                       meanFC=all_meanLogFC_TAD, 
                       comb_pval=emp_pval_combined)
                                   



plot_volcanoTADsCorrFC <- function(meanCorr, meanFC, comb_pval, padjusted=FALSE,
                                   x_lab="TAD mean LogFC",
                                   y_lab="TAD mean intraCorr",
                                   plotTit=NULL,
                                   subTit=NULL,
                                   fcUpDownThresh=0,
                                   signifThresh=0.05,
                                   verySignifThresh=0.01,
                                   strongUp="#EE0011FF",
                                   lightUp="#EE001132",
                                   strongDown="#0C5BB0FF",
                                   lightDown="#0C5BB032",
                                   fontFamily="Hershey",
                                   tads_to_annot=NULL) {
  
  if(!padjusted) {
    cat("... adjust p-values with BH method\n")
    adj_comb_pval <- p.adjust(comb_pval, method="BH")
  } else {
    adj_comb_pval <- comb_pval
  }
  
  interReg <- Reduce(intersect, list(names(meanCorr) , names(meanFC), names(comb_pval)))
  keepCorr <- meanCorr[interReg]
  keepFC <- meanFC[interReg]
  keepComb <- comb_pval[interReg]
  
  cat(paste0("... kept meanCorr values:\t", length(keepCorr), "/", length(meanCorr)))
  cat(paste0("... kept meanFC values:\t", length(keepFC), "/", length(meanFC)))
  cat(paste0("... kept comb. p-val values:\t", length(keepComb), "/", length(comb_pval)))
  
  
  ex_DT <- data.frame(
    region = interReg,
    meanLogFC = keepFC[interReg],
    meanCorr = keepCorr[interReg],
    adjPvalComb = keepComb[interReg],
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(ex_DT))
  ex_DT$adjPvalComb_log10 <- -log10(ex_DT$adjPvalComb)
  
  
  nSignif <- sum(ex_DT$adjPvalComb<=signifThresh)
  nVerySignif <- sum(ex_DT$adjPvalComb<=verySignifThresh)
  
  if(is.null(plotTit)) plotTit <- paste0("TAD meanLogFC and meanCorr")
  if(is.null(subTit)) subTit <- paste0("# TADs=", nrow(ex_DT), "; # signif.<=", signifThresh, "=",nSignif,"; # signif.<=", verySignifThresh, "=",nVerySignif)
  
  ex_DT$adjPvalComb_log10_col <- ex_DT$adjPvalComb_log10
  ex_DT$adjPvalComb_log10_col[ex_DT$adjPvalComb_log10_col <= -log10(signifThresh)] <- -log10(signifThresh)
  
  ex_DT$adjPvalComb_crop <- ex_DT$adjPvalComb 
  ex_DT$adjPvalComb_crop[ex_DT$adjPvalComb_crop >= signifThresh] <- signifThresh
  
  ex_DT$adjPvalCombSignif <- ex_DT$adjPvalComb <= signifThresh
  
  
  ex_DT$dotcol_labs1 <- ifelse(ex_DT$meanLogFC > fcUpDownThresh, paste0("logFC>", fcUpDownThresh), paste0("logFC<=", fcUpDownThresh))
  ex_DT$dotcol_labs2 <- ifelse(ex_DT$adjPvalComb <= verySignifThresh, paste0("adjP<=", verySignifThresh),
                               ifelse(ex_DT$adjPvalComb <= signifThresh, paste0("adjP<=", signifThresh), NA))
  
  ex_DT$dotcol_lab <- interaction(ex_DT$dotcol_labs1, ex_DT$dotcol_labs2)
  ex_DT$dotcol_lab <- as.character(ex_DT$dotcol_lab)
  stopifnot(which(is.na(ex_DT$dotcol_labs2)) == which(ex_DT$adjPvalComb > signifThresh))
  ex_DT$dotcol_lab[is.na(ex_DT$dotcol_labs2)] <- paste0("adjP>", signifThresh)
  stopifnot(!is.na(ex_DT$dotcol_lab))
  
  mycols <- setNames(c(strongUp, lightUp, strongDown, lightDown, "darkgrey"), 
                     c(paste0("logFC>", fcUpDownThresh, ".adjP<=", verySignifThresh),
                       paste0("logFC>", fcUpDownThresh, ".adjP<=", signifThresh),
                       paste0("logFC<=", fcUpDownThresh, ".adjP<=", verySignifThresh),
                       paste0("logFC<=",fcUpDownThresh, ".adjP<=", signifThresh), 
                       paste0("adjP>", signifThresh)))
  
  mycols_labs <- setNames(c(paste0("mean LogFC > ", fcUpDownThresh, " &\nadj. p-val <=", verySignifThresh),
                            paste0("mean LogFC > ", fcUpDownThresh, " &\nadj. p-val <=", signifThresh),
                            paste0("mean LogFC <= ", fcUpDownThresh, " &\nadj. p-val <=", verySignifThresh),
                            paste0("mean LogFC <= ", fcUpDownThresh, " &\nadj. p-val <=", signifThresh), 
                            paste0("adj. p-val > ", signifThresh)),
                          c(
                            paste0("logFC>",fcUpDownThresh, ".adjP<=", verySignifThresh),
                            paste0("logFC>", fcUpDownThresh, ".adjP<=", signifThresh),
                            paste0("logFC<=", fcUpDownThresh, ".adjP<=", verySignifThresh),
                            paste0("logFC<=", fcUpDownThresh, ".adjP<=", signifThresh), 
                            paste0("adjP>", signifThresh)))
  
  fc_corr_p <- ggplot(ex_DT, aes(x=meanLogFC, y = meanCorr, fill = dotcol_lab, color=dotcol_lab )) +
    ggtitle(plotTit, subtitle = subTit)+
    scale_fill_manual(values=mycols, labels = mycols_labs)+  
    scale_color_manual(values=mycols,
                       labels=setNames(names(mycols), mycols), guide=F)+
    guides(color=F)+
    
    geom_point(data=ex_DT[ex_DT$dotcol_lab==paste0("adjP>", signifThresh),],
               shape=21,
               aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1) +
    geom_point(data=ex_DT[ex_DT$dotcol_lab!=paste0("adjP>", signifThresh),],
               shape=21,
               aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1)+
    
    labs( size= "TAD adj. p-val", fill="")+
    geom_label_repel(data = ex_DT[ex_DT$region %in% tads_to_annot,],
                     aes(x= meanLogFC, 
                         y=meanCorr, label=region),
                     inherit.aes = F,
                     min.segment.length = unit(0, 'lines'),
                     force = 10,
                     segment.size = 0.5,
                     nudge_y = 0.1)+
    scale_x_continuous(name=paste0(x_lab),breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name=paste0(y_lab),breaks = scales::pretty_breaks(n = 10))+
    scale_size(range=c(0.5,8),#expand=c(2,0),
               breaks=c(-log10(0.05), -log10(0.01), -log10(0.005), -log10(0.001)),
               labels=c(paste0(">= 0.05"),"0.01"," 0.005","0.001"),
               guide="legend") +
    guides(size = guide_legend(order=1))+
    guides(fill = guide_legend(override.aes = list(color="white", size=4),order=2))+
    theme(
      text = element_text(family=fontFamily),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      axis.line = element_line(),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
      panel.grid = element_blank(),
      legend.key = element_rect(fill = NA),
      legend.title =  element_text(size=12,face="bold"),
      panel.grid.major.y =  element_blank(),
      panel.grid.minor.y =  element_blank(),
      legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=10)
    ) +
    geom_hline(yintercept = 0, linetype=2, color="darkgrey")+
    geom_vline(xintercept = 0, linetype=2, color="darkgrey")
  fc_corr_p
  return(fc_corr_p)
}



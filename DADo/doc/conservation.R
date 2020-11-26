## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
rm(list=ls())

if(!require(DADo))
  devtools::install_github("CSOgroup/DADo", subdir="DADo")
  # alternatively: 
  # install.packages("DADo_0.0.0.1.tar.gz", repos = NULL, type ="source")
 # data("norm_ID")
library(DADo)

## ----get_conserved, fig.height=8, fig.width=14---------------------------
data("allSignif_dt") # this loads allSignif_dt
head(allSignif_dt)

data("allSignif_g2t_dt") # this loads allSignif_g2t_dt
head(allSignif_g2t_dt)


data("allSignif_tad_pos_dt") # this loads allSignif_dt
head(allSignif_tad_pos_dt)

conserv_minIntersectGenes <- 3
conserv_geneMatching <- 0.8
conserv_minOverlapBp <- 0.8

signif_conserv_data <- get_conservedRegion(
  signif_dt=allSignif_dt, all_tad_pos_dt=allSignif_tad_pos_dt, all_g2t_dt=allSignif_g2t_dt, 
  minOverlapBpRatio=conserv_minOverlapBp,
  minIntersectGenes=conserv_minIntersectGenes,
  gene_matching_fuse_threshold = conserv_geneMatching,
  verbose=FALSE
)
save(signif_conserv_data, file="package_signif_conserv_data.Rdata", version=2)

## ----plot_barplot, fig.height=6, fig.width=10----------------------------
library(ggplot2)

minConserv <- 2

conserved_dt <- data.frame(
      region = names(signif_conserv_data[["conserved_signif_tads"]]),
      nConserv = as.numeric(lengths(signif_conserv_data[["conserved_signif_tads"]])),
      stringsAsFactors=FALSE)
stopifnot(conserved_dt$nConserv >= 2)

conserved_dt <- conserved_dt[order(conserved_dt$nConserv),]
conserved_dt$region <- factor(conserved_dt$region, levels=as.character(conserved_dt$region))
conserved_dt$region_rank <- as.numeric(conserved_dt$region)
conserved_dt <- conserved_dt[conserved_dt$nConserv >= minConserv,]

nCons <- nrow(conserved_dt)

my_ylab <- "# datasets conserved"
my_xlab <- paste0("Ranked conserved regions (", nCons, " with >= ", minConserv, " DS cons.)")
myTit <- "Signif. conserved regions"
subTit <- paste0("conserv. match ratio bp >= ", conserv_minOverlapBp, " and intersect genes >= ", conserv_minIntersectGenes)

ggplot(conserved_dt, aes(x=region_rank, y=nConserv)) +
  ggtitle(myTit, subtitle = subTit)+
  geom_bar(stat="identity", color="darkblue", fill="gray80") +
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=5), expand=c(0,0))+
  labs(x=my_xlab, y=my_ylab, color="", fill="")+
   theme(
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
	    axis.line.x = element_line(),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold"),
	    legend.text = element_text(size=12),
)

## ----plot_heatmap, fig.height=14, fig.width=14---------------------------
library(reshape2)

myTit <- "Conserved regions across datasets"

# transform and plot to a binary matrix
all_datasets <- as.character(unique(allSignif_dt$dataset))

dsByReg_dt <- do.call(rbind, lapply(signif_conserv_data[["conserved_signif_tads"]], function(x) as.numeric(all_datasets %in% unique(dirname(x)))))
colnames(dsByReg_dt) <- all_datasets
stopifnot(max(rowSums(dsByReg_dt)) == max(conserved_dt$nConserv))

plot_dt <- dsByReg_dt[names(sort(rowSums(dsByReg_dt), decreasing = TRUE)),names(sort(colSums(dsByReg_dt), decreasing = TRUE))]
plot_dt_m <- melt(plot_dt)
plot_dt_m$Var1 <- factor(plot_dt_m$Var1, levels=names(sort(rowSums(dsByReg_dt), decreasing = TRUE)))
stopifnot(!is.na(plot_dt_m$Var1))
plot_dt_m$Var2 <- factor(plot_dt_m$Var2, levels=names(sort(colSums(dsByReg_dt), decreasing = TRUE)))
stopifnot(!is.na(plot_dt_m$Var2))
plot_dt_m$value <- as.character(plot_dt_m$value)

ggplot(plot_dt_m, aes(x=Var2, y=Var1, fill=value)) + 
  ggtitle(myTit)+
  geom_raster()+
  labs(x="Ranked datasets", y="Ranked conserved regions")+
  guides(fill=F)+
  scale_fill_manual(values=setNames(c("black", "white"), c("1", "0")))+
  theme(
    plot.title=element_text(size=14, face="bold", hjust=0.5),
    axis.text = element_blank()
  )

## ----plot_conserved, fig.height=8, fig.width=14--------------------------
region_to_plot <- conserved_dt$region[conserved_dt$nConserv == max(conserved_dt$nConserv)][1]
stopifnot(region_to_plot %in% names(signif_conserv_data[["conserved_signif_intersect_genes"]]))
symbols_to_plot <- signif_conserv_data[["conserved_signif_intersect_genes"]][[paste0(region_to_plot)]]
stopifnot(region_to_plot %in% names(signif_conserv_data[["conserved_signif_tads"]]))
tads_to_plot <- signif_conserv_data[["conserved_signif_tads"]][[paste0(region_to_plot)]]
stopifnot(symbols_to_plot %in% allSignif_g2t_dt$symbol)

genes_plot_dt <- allSignif_g2t_dt[allSignif_g2t_dt$symbol %in% symbols_to_plot,c("symbol", "start", "end", "chromo")]
genes_plot_dt <- unique(genes_plot_dt)
tads_plot_dt <- data.frame(dataset = dirname(tads_to_plot),
                           region = basename(tads_to_plot), regID = tads_to_plot,
                           stringsAsFactors=FALSE)
allSignif_tad_pos_dt$regID <- file.path(allSignif_tad_pos_dt$dataset, allSignif_tad_pos_dt$region)
stopifnot( tads_to_plot %in% allSignif_tad_pos_dt$regID)
allSignif_tad_pos_dt$regID <- file.path(allSignif_tad_pos_dt$dataset, allSignif_tad_pos_dt$region)
tads_plot_dt <- allSignif_tad_pos_dt[allSignif_tad_pos_dt$regID %in% tads_to_plot,]
tads_plot_dt <- unique(tads_plot_dt)
tads_plot_dt$dsCat <- 1 # draw all 
tads_plot_dt$cond1 <- gsub("TCGA.+_(.+)_.+", "\\1", basename(tads_plot_dt$dataset))
tads_plot_dt$cond2 <- gsub("TCGA.+_.+_(.+)", "\\1", basename(tads_plot_dt$dataset))
tads_plot_dt$dataset <- gsub("_40kb", "", dirname(as.character(tads_plot_dt$dataset)))
plot_conservedRegions(genes_dt=genes_plot_dt, 
                      tads_dt=tads_plot_dt)


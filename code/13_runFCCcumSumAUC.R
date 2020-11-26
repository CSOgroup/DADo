#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script11: all_FCC_TAD.RData
# - script12: FCC_permDT.RData
################################################################################

################  OUTPUT
# fcc_auc_ratios.RData (+ wave plot)
################################################################################
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(".")

script11_name <- "11_runFCC"
script12_name <- "12_runPermutationsFCC"
script_name <- "13_runFCCcumSumAUC"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)

source(file.path(pipScriptDir, "TAD_DE_utils.R"))

registerDoMC(nCpu) # loaded from main_settings.R

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


# if not set: default is to run for gene permutation and not for randomTADs
# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


qtPermutAUC <- 0.95

### SOME SETTINGS FOR THE PLOT
# set colors for plotting:
plotType <- "png"
myHeight <- 600
myWidth <- 800
plotCex <- 1.2
polygonPermutCol <- rgb(188/255,188/255,188/255, 0.3)
qt95PermutCol <- rgb(140/255,140/255,140/255)
lwdObs <- 1.2
colObs <- "darkred"
pointObsCol <- "darkred"

my_main <- paste0("Genome-wide intra-TAD fold-change concordance")
my_sub <- paste0(cond1, " vs. " , cond2)
my_xlab <- paste0("TADs ranked by FCC")
my_ylab <- paste0("FCC cumulative sum")


obs_fcc <- eval(parse(text = load(file.path(pipOutFold, script11_name, "all_FCC_TAD.RData"))))
permut_FCC <- eval(parse(text = load(file.path(pipOutFold, script12_name,  "FCC_permDT.RData"))))

if(ncol(permut_FCC) != nRandomPermut)
  stop(paste0("! NEED TO CHECK: different settings were used for running the permutations !\nncol(permut_FCC) != nRandomPermut\nnncol(permut_FCC)\t=\t", nncol(permut_FCC),"\nnRandomPermut\t=\t", nRandomPermut, "\n"))

nPermut <- ncol(permut_FCC)

my_sub <- paste0(my_sub, " - # permut=", nPermut)

stopifnot(all(names(obs_fcc) %in% rownames(permut_FCC)))
stopifnot(all(rownames(permut_FCC) %in% names(obs_fcc)))
stopifnot(setequal(rownames(permut_FCC), names(obs_fcc)))
interReg <- intersect(names(obs_fcc),rownames(permut_FCC) )

      
############################################################
############################################################ # filter the TADs and sort
############################################################
filter_obs_fcc <- sort(obs_fcc[interReg], decreasing = T)

filter_permut_FCC_unsort <- permut_FCC[interReg,]
stopifnot(length(filter_obs_fcc) == nrow(filter_permut_FCC_unsort))

filter_permut_FCC <- apply(filter_permut_FCC_unsort, 2, sort, decreasing=T)
rownames(filter_permut_FCC) <- NULL
stopifnot(length(filter_obs_fcc) == nrow(filter_permut_FCC_unsort))

maxRankPlot <- ceiling(length(interReg)/1000)*1000

x_val <- c(1:length(filter_obs_fcc))
# CHANGED HERE 28.04.2020 -> do not take abs !
# cumsum_permut_dt <- apply(filter_permut_FCC, 2, function(x) cumsum(abs(x-0)))
cumsum_permut_dt <- apply(filter_permut_FCC, 2, function(x) cumsum(x))
meanPermut_cumsum <- apply(cumsum_permut_dt, 1, mean)
qt95Permut_cumsum <- apply(cumsum_permut_dt, 1, quantile, probs=qtPermutAUC) # for each rank, take the quantile
# CHANGED HERE 28.04.2020 -> do not take abs !
# obs_cumsum <- cumsum(abs(filter_obs_fcc - 0))
obs_cumsum <- cumsum(filter_obs_fcc)

auc_obs <- auc(x = x_val, y = obs_cumsum)
auc_permutMean <- auc(x = x_val, y = meanPermut_cumsum)
auc_permutQt <- auc(x = x_val, y = qt95Permut_cumsum)
auc_ratioQt <- auc_obs/auc_permutQt


fcc_auc_ratios <- list(
  auc_obs=auc_obs,
  auc_permutMean=auc_permutMean,
  auc_permutQt=auc_permutQt
)

saveFile <- file.path(curr_outFold, paste0("fcc_auc_ratios.RData"))      
save(fcc_auc_ratios, file= saveFile)
cat(paste0("... written: ", saveFile, "\n"))
      

pct_inc_qt <- round(auc_obs/auc_permutQt,2)
      
outFile <- file.path(curr_outFold, paste0("genomeWide_FCC_cumsum_obs_permut.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
par(mar = par()$mar + c(0,1,0,0))
plot(obs_cumsum ~ x_val,
   main= my_main,
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
   bty="l")
box()
axis(2, cex.axis=plotCex, cex.lab=plotCex)
axis(1, cex.axis=plotCex, cex.lab=plotCex, at = seq(from=0, to=maxRankPlot, by=maxRankPlot/10)) # we used 2000 and 200 hard-coded
mtext(my_sub, font=3)
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
     legend = c(paste0("observed (n=", length(filter_obs_fcc), ")"), "min-max permut.", paste0(qtPermutAUC, "-qt permut.")), 
     pt.cex = c(0.7, 2),
     col = c(pointObsCol, polygonPermutCol, qt95PermutCol),
     bty="n")
legtxt <- as.expression(bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc_qt)))
legend("right", legend=c(legtxt), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))








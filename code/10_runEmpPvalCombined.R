#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

### !!! UPDATE 24.07.19: STOUFFER ONE-SIDED

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script3: all_meanLogFC_TAD.RData
# - script8: all_obs_ratioDown.RData
# - script9: emp_pval_meanLogFC.RData
# - script10sameNbr: emp_pval_meanCorr.RData
################################################################################

################  OUTPUT
# - emp_pval_combined.RData 
################################################################################

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")


script8_name <- "8_runEmpPvalMeanTADLogFC"
script9_name <- "9_runEmpPvalMeanTADCorr"
script_name <- "10_runEmpPvalCombined"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

twoTailsStouffer <- FALSE # discussion with Marco -> do it one-sided (compared two- and one- -> similar)
# stouffer(c(emp_pval_intraCorr[x], emp_pval_logFC[x]), two.tails = twoTailsStouffer)))

# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> twoTailsStouffer\t=\t", as.character(twoTailsStouffer), "\n")
printAndLog(txt, pipLogFile)

stopifnot(exists("all_permutCorr_data"))
if(dir.exists(all_permutCorr_data)){
  txt <- paste0(toupper(script_name), "> use emp p-val. from folder\t=\t", all_permutCorr_data, "\n")
  printAndLog(txt, pipLogFile)
  filePrefix <- "fromFolder_"
} else if(file.exists(all_permutCorr_data)){
  txt <- paste0(toupper(script_name), "> use emp p-val. from file\t=\t", all_permutCorr_data, "\n")
  printAndLog(txt, pipLogFile)
  filePrefix <- "fromFile_"
}



################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# load emp. p-val logFC 
fc_file <- file.path(pipOutFold,  script8_name, "emp_pval_meanLogFC.RData")
emp_pval_logFC <- eval(parse(text = load(fc_file)))
txt <- paste0(toupper(script_name), "> retrieve emp. p-val. for FC from\t=\t", fc_file, "\n")
printAndLog(txt, pipLogFile)

# load emp. p-val intraTAD corr
corr_file <- file.path(pipOutFold, script9_name, paste0(filePrefix, "emp_pval_meanCorr.RData"))
emp_pval_intraCorr <- eval(parse(text = load(corr_file)))
txt <- paste0(toupper(script_name), "> retrieve emp. p-val. for corr. from\t=\t", corr_file, "\n")
printAndLog(txt, pipLogFile)

intersectRegions <- sort(intersect(names(emp_pval_logFC), names(emp_pval_intraCorr)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)
if(length(intersectRegions) < length(emp_pval_logFC) | length(intersectRegions) < length(emp_pval_intraCorr)) {
  txt <- paste0(toupper(script_name), "> !!! WARNING: regions for empirical p-val logFC and corr. do not match !!!\n")
  printAndLog(txt, pipLogFile)    
}

### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", names(emp_pval_logFC))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val logFC data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    if( length(grep("_TAD", names(emp_pval_intraCorr))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val intraCorr data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}

initLen <- length(emp_pval_logFC)
emp_pval_logFC <- emp_pval_logFC[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val logFC: ", length(emp_pval_logFC), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initLen <- length(emp_pval_intraCorr)
emp_pval_intraCorr <- emp_pval_intraCorr[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val intraCorr: ", length(emp_pval_intraCorr), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(is.na(emp_pval_logFC)))
stopifnot(!any(is.na(emp_pval_intraCorr)))

stopifnot(names(emp_pval_logFC) == names(emp_pval_intraCorr))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL INTRA-TAD CORR & WRITE OUTPUT
################################****************************************************************************************

# COMBINE EMPIRICAL P-VALUES
emp_pval_combined <- unlist(sapply(seq_along(intersectRegions), function(x) 
                  stouffer(c(emp_pval_intraCorr[x], emp_pval_logFC[x]), two.tails = twoTailsStouffer)))
names(emp_pval_combined) <- intersectRegions

stopifnot(length(emp_pval_combined) == length(intersectRegions))

outFile <- file.path(curr_outFold, paste0(gsub("_", "", filePrefix), "PermCorr_emp_pval_combined.RData"))
save(emp_pval_combined, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

# added for the release: 7.7.2020
outFile <- file.path(curr_outFold, paste0(gsub("_", "", filePrefix), "PermCorr_adj_emp_pval_combined.RData"))
adj_emp_pval_combined <- p.adjust(emp_pval_combined, method="BH")
save(adj_emp_pval_combined, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))



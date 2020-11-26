#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script4: all_meanCorr_TAD.RData
# - script7sameNbr: meanCorr_sample_around_TADs_sameNbr.RData
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.RData
################################################################################


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- file.path(".")

script4_name <- "4_runMeanTADCorr"
script7_name <- "7_runPermutationsMeanTADCorr"
script_name <- "9_runEmpPvalMeanTADCorr"
stopifnot(file.exists(file.path(pipScriptDir, paste0(script_name, ".R"))))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

registerDoMC(nCpu) # loaded from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE


# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- file.path(pipOutFold, paste0(format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt"))
system(paste0("rm -f ", pipLogFile))

# ADDED 16.11.2018 to check using other files
txt <- paste0(toupper(script_name), "> inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


corr_type <- "meanCorr"
txt <- paste0(toupper(script_name), "> taking sample correlation for corr_type\t=\t", corr_type, "\n")
printAndLog(txt, pipLogFile)

stopifnot(exists("all_permutCorr_data"))


if(dir.exists(all_permutCorr_data)){
  ### RETRIEVE ALL THE FILES IN THE FOLDER !!!
  txt <- paste0(toupper(script_name), "> !!! take all the files recursively matching \"", corrMatchPattern, "\" in ", all_permutCorr_data, "\n")
  printAndLog(txt, pipLogFile)
  all_sampleCorr_files <- list.files(all_permutCorr_data, pattern=corrMatchPattern, full.names = TRUE, recursive = TRUE)
  if(exists("refineMatchPattern")) {
	all_sampleCorr_files <- all_sampleCorr_files[grepl(refineMatchPattern, all_sampleCorr_files)]
  }
  all_sampleCorr_files <- all_sampleCorr_files[!grepl(corrDiscardPattern, all_sampleCorr_files)]
  txt <- paste0(toupper(script_name), "> sampleCorr_files used (n=", length(all_sampleCorr_files), "):\n")
  printAndLog(txt, pipLogFile)
  txt <- paste0(paste0(all_sampleCorr_files, collapse="\n"), "\n")
  printAndLog(txt, pipLogFile)
  if(exists("nbrCorrPermutCheck")) {
	stopifnot(is.numeric(nbrCorrPermutCheck))
	stopifnot(length(all_sampleCorr_files) == nbrCorrPermutCheck)
  }
  if(length(all_sampleCorr_files) == 0) stop("could not find any file with correlations from permutation data")
  all_permut_corrValues <- foreach(corr_file = all_sampleCorr_files, .combine='c') %dopar% {
    stopifnot(file.exists(corr_file))
    corr_data <- eval(parse(text = load(corr_file)))
    all_samp_corrs <- as.numeric(sapply(corr_data, function(x) x[[paste0(corr_type)]]))
    stopifnot(!is.null(all_samp_corrs))
    all_samp_corrs <- na.omit(all_samp_corrs)  
    all_samp_corrs
  }
  filePrefix <- "fromFolder_" 
} else if(file.exists(all_permutCorr_data)) {
  txt <- paste0(toupper(script_name), "> use provided all_permutCorr_data\t=\t",all_permutCorr_data, "\n")
  printAndLog(txt, pipLogFile)
  all_permut_corrValues <- get(load(file.path(all_permutCorr_data)))
  if(exists("nbrCorrPermutCheck")) {
	stopifnot(is.numeric(nbrCorrPermutCheck))
	stopifnot(length(all_permut_corrValues) == nbrCorrPermutCheck)
  }
  all_permut_corrValues <- unlist(all_permut_corrValues)
  stopifnot(is.numeric(all_permut_corrValues))
  filePrefix <- "fromFile_" 
} else  {
  stop(paste0("error with permutation correlation data\nall_permutCorr_data should be a file or a folder\nall_permutCorr_data\t=\t",all_permutCorr_data,"\n"))
}
txt <- paste0(toupper(script_name), "> # of permutation correlation values\t=\t", length(all_permut_corrValues),"\n")
printAndLog(txt, pipLogFile)


# RETRIEVE THE OBSERVED CORR DATA
obs_corr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.RData")
stopifnot(file.exists(obs_corr_file))
all_obs_corr <- eval(parse(text = load(obs_corr_file)))
 
emp_pval_meanCorr <- sapply(all_obs_corr, function(x) {
  (sum(all_permut_corrValues >= x) + 1)/(length(all_permut_corrValues) + 1)
})
names(emp_pval_meanCorr) <- names(all_obs_corr)
stopifnot(all(emp_pval_meanCorr > 0 & emp_pval_meanCorr <= 1 ))


outFile <- file.path(curr_outFold, paste0(filePrefix, "emp_pval_meanCorr.RData"))
save(emp_pval_meanCorr, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))
       
          
        
        

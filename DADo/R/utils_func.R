#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Print text and log
#'
#' Function that "cat" a message and print it to log file.
#'
#' @param mytext The text to be printed.
#' @param mylogf Location where to print the text (e.g. log file or "").
#' @export
printAndLog <- function(mytext, mylogf) {
  cat(mytext)
  cat(mytext, append=T, file = mylogf)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Return the first lines/columns of a dataframe
#'
#' Returns the first lines and columns of a dataframe (as many columns as lines).
#'
#' @param dt A two-dimensional dataframe or matrix.
#' @param ntop The maximal number of lines/columns to display.
#' @return The first rows/columns of the dataframe.
#' @export
head_sq <- function(dt, ntop=5) {
  nRow <- min(c(ntop, nrow(dt)))	
  nCol <- min(c(ntop, ncol(dt)))
  dt[seq_len(nRow),seq_len(nCol)]
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Print text and log
#'
#' Function that "cat" or not a message and print it or not to log file.
#'
#' @param txt The text to be printed.
#' @param verbose If the text should be "cat"."
#' @param logFile File where to append the text (if not NULL)
#' @export
outTxt <- function(txt, verbose, logFile) {
	if(verbose)   cat(txt)
	if(!is.null(logFile))  cat(txt, file=logFile, append=TRUE)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Adapt pretty breaks with no zeros
#'
#' Adaptation to scales::pretty breaks to remove the zeros of the scale.
#'
#' @param n The number of breaks wanted.
#' @param ... Other parameters passed to scales::force_all.
#' @return The breaks without zeros
#' @export

noZero_breaks <- function (n = 5, ...) {
  library(scales)
  scales:::force_all(n, ...)
  function(x) {
    breaks <- pretty(x, n, ...)
    breaks <- breaks[breaks > 0]
    names(breaks) <- c(attr(breaks, "labels"))
    c(1,breaks)
  }
}


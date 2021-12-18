#' @title Calculate on-target sgRNA activity scores for enCas12a using enPAM+GB
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas12a-induced knockout using the enPAM+GB scoring method.
#'     Currently not supported on Windows machines.
#' 
#' @param sequences Character vector of 34bp sequences needed for enPAM+GB
#'     scoring, see details below.
#' 
#' @details The input sequences for enPAM+GB scoring require 4 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (4bp PAM sequence + 23bp spacer sequence) and 3 nucleootides 
#'     downstream of the protospacer sequence, for a total of 34 nucleotides.
#'     Both canonical and non-canonical PAM sequences can be provided.
#' 
#' @return \strong{getEnPAMGBScores} returns a data.frame with \code{sequence}
#'     and \code{score} columns.
#'
#' @references 
#' DeWeirdt, P.C., Sanson, K.R., Sangree, A.K. et al. Optimization of AsCas12a
#'     for combinatorial genetic screens in human cells.
#'     Nat Biotechnol 39, 94–104 (2021).
#'     \url{https://doi.org/10.1038/s41587-020-0600-6}.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' if (interactive()){
#' flank5 <- "CATG" #4bp
#' pam    <- "TTTT" #4bp
#' spacer <- "TTTGGGAACCAATCGATAATCAC" #23bp
#' flank3 <- "ATT" #3bp
#' input  <- paste0(flank5, pam, spacer, flank3) 
#' results <- getEnPAMGBScores(input)
#' }
#' @inheritParams getAzimuthScores
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getEnPAMGBScores <- function(sequences, fork=FALSE){
    if (.Platform$OS.type=="windows"){
        stop("EnPAMGB is not available for Windows at the moment.")
    }
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=34){
        stop("Provided sequences must have length 34nt",
             " ([4nt][TTTV][23mer][3nt]).")
    }
    results <- basiliskRun(env=env_enpamgb,
                           shared=FALSE,
                           fork=fork,
                           fun=.enpamgb_python, 
                           sequences=sequences)
    return(results)
}

#' @importFrom reticulate import_from_path
#' @importFrom reticulate np_array
.enpamgb_python <- function(sequences){

    dir <- system.file("python",
                       "enpamgb",
                       package="crisprScore",
                       mustWork=TRUE)
    enpamgb <- import_from_path("getEnPAMGB", path=dir)
    
    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]
    if (length(sequences.valid)>0){
        scores <- enpamgb$getEnPAMGB(np_array(sequences.valid))
        scores[scores>1]  <- 1
        scores[scores<=0] <- 0
        df$score[good] <- scores    
    }
    return(df)
}


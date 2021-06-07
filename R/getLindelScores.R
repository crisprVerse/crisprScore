#' @title Predict frameshift ratios from CRISPR/Cas9 indel prediction using
#'     Lindel
#' @description Predict frameshift ratios from CRISPR/Cas9 indel prediction
#'     using the Lindel prediction algorithm.
#' 
#' @param sequences Character vector of 65bp sequences needed for Lindel
#'     scoring, see details below. 
#' 
#' @details The input sequences for Lindel scoring require 13 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (23 nucleotides) and 29 nucleootides downstream of the protospacer
#'     sequence, for a total of 65 nucleotides. Note that only canonical PAM
#'     sequences (NGG) are accepted by Lindel.
#' 
#' @return A data.frame with predicted frameshift ratio (between 0 and 1).
#'     A higher ratio indicates a greater chance of a frameshift indel
#'     introduced by CRISPR/Cas9-induced double-strand breaks. 
#' 
#' @references
#' Wei Chen, Aaron McKenna, Jacob Schreiber, Maximilian Haeussler, Yi Yin,
#'     Vikram Agarwal, William Stafford Noble, Jay Shendure, Massively parallel
#'     profiling and predictive modeling of the outcomes of CRISPR/Cas9-mediated
#'     double-strand break repair, Nucleic Acids Research, Volume 47, Issue 15,
#'     05 September 2019, Pages 7989–8003,
#'     \url{https://doi.org/10.1093/nar/gkz487}.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' \donttest{
#' flank5 <- "ACCTTTTAATCGA" #13bp
#' spacer <- "TGCTGATGCTAGATATTAAG" #20bp
#' pam    <- "TGG" #3bp
#' flank3 <- "CTTTTAATCGATGCTGATGCTAGATATTA" #29bp
#' input <- paste0(flank5, spacer, pam, flank3)
#' results <- getLindelScores(input)
#' }
#' @inheritParams getAzimuthScores
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getLindelScores <- function(sequences, fork=FALSE){
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=65){
        stop("Sequences must have length 65nt ([33nt][NGG][29nt]).")
    }
    results <- basiliskRun(env=env_lindel,
                           shared=FALSE,
                           fork=fork,
                           fun=.lindel_python,
                           sequences=sequences)
    return(results)
}






#' @importFrom stringr str_extract
.lindel_python <- function(sequences){
    program <- system.file("python",
                           "lindel",
                           "Lindel_prediction.py",
                           package="crisprScore",
                           mustWork=TRUE)

    #weights_dir <- "/Users/fortinj2/crisprScoreData/inst/scripts/out/Lindel"
    weights_dir  <- system.file(package="crisprScoreData", "scripts/out/Lindel")
    weights_file <- file.path(weights_dir, "Model_weights.pkl")
    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]
    if (length(sequences.valid)>0){
        scores <- rep(NA_real_, length(sequences.valid))
        for (i in seq_along(sequences.valid)){
            dir  <- tempdir()
            file <- basename(tempfile())
            file.full <- file.path(dir, file)
            seq <- sequences.valid[i]
            cmd <- paste0("python ", program, " ",
                          seq, " ",
                          weights_file, " ",
                          file.full)
            system(cmd, ignore.stdout=TRUE, ignore.stderr=FALSE)
            outputs <- list.files(dir)
            outputs <- outputs[grepl(file, outputs)]
            file.remove(file.path(dir,outputs))
            score <- str_extract(outputs, "fs=[0-9]+\\.[0-9]+\\.txt")
            score <- as.numeric(gsub("fs=|.txt", "", score))
            score <- as.numeric(score)
            scores[i] <- score
        }
        df$score[good] <- scores
    }
    return(df)
}





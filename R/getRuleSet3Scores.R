#' @title Calculate on-target sgRNA activity scores for SpCas9 using Rule Set 3
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas9-induced knockout using the Rule Set 3 scoring method.
#' 
#' @param sequences Character vector of 30bp sequences needed for Rule Set 3
#'     scoring, see details below.
#' @param tracrRNA String specifying which tracrRNA is used. 
#'     Must be either "Hsu2013" (default) or "Chen2013".
#' @param mode String specifying which prediction mode is used.
#'     Must be either "sequence" (default) or "target".
#' 
#' @details The input sequences for Rule Set 3 scoring require 4 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (20bp spacer sequence + 3bp PAM sequence ) and 3 nucleootides 
#'     downstream of the protospacer sequence, for a total of 30 nucleotides.
#' 
#' @return \strong{getRuleSet3Scores} returns a data.frame with
#'     \code{sequence} and \code{score} columns. The getRuleSet3Scores score
#'     is similar to a Z-score.A higher score indicates higher
#'     knockout efficiency.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @references doi: https://doi.org/10.1101/2022.06.27.497780
#' 
#' @examples
#' 
#' if (interactive()){
#' flank5 <- "ACCG" #4bp
#' spacer <- "AATCGATGCTGATGCTAGAT" #20bp
#' pam    <- "AGG" #3bp
#' flank3 <- "AAT" #3bp
#' input  <- paste0(flank5, spacer, pam, flank3) 
#' results <- getRuleSet3Scores(input)
#' }
#' 
#' @export 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getRuleSet3Scores <- function(sequences,
                              tracrRNA=c("Hsu2013","Chen2013"),
                              mode=c("sequence", "target")
){

    tracrRNA <- match.arg(tracrRNA)
    mode <- match.arg(mode)
    if (mode=="sequence"){
        results <- .getRuleSet3Scores_sequence(sequences=sequences,
                                               tracrRNA=tracrRNA)
    } else if (mode=="target"){
        results <- .getRuleSet3Scores_target(sequences=sequences,
                                             tracrRNA=tracrRNA)
    }
    return(results)
}

.getRuleSet3Scores_sequence <- function(sequences,
                                        tracrRNA=c("Hsu2013","Chen2013")
){
    tracrRNA <- match.arg(tracrRNA)
    .dumpToFile <- function(sequences,
                            file){
        write.table(sequences,
                    quote=FALSE,
                    row.names=FALSE,
                    col.names=FALSE,
                    file=file)
    }

    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=30){
        stop("Provided sequences must have length 30nt",
             " ([4nt][20mer][NGG][3nt]).")
    }
  
    #Output data.frame
    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]

    #Saving to disk:
    dir <- tempdir()
    inputfile  <- file.path(dir, "input.txt")
    outputfile <- file.path(dir, "output.txt")
    
    # Ready to get the scores
    env <- basilisk:::.obtainEnvironmentPath(env_rs3)
    basilisk.utils::activateEnvironment(env)
    programFile <- system.file("python",
                               "rs3/getRuleSet3ScoresSequence.py",
                               package="crisprScore",
                               mustWork=TRUE)
    modelFile <- system.file("rs3",
                             "RuleSet3.pkl",
                             package="crisprScoreData")
    cmd <- paste0("python ",
                  programFile, " ",
                  inputfile, " ",
                  modelFile, " ",
                  outputfile, " ",
                  shQuote(tracrRNA))
    if (sum(good)>0){
        .dumpToFile(sequences.valid,
                    inputfile)
        system(cmd)
        #scores <- readLines(outputfile, sep="\t")[,4]
        scores <- readLines(outputfile)
        df$score[good] <- scores
    }
    if (file.exists(inputfile)){
        file.remove(inputfile)
    }
    if (file.exists(outputfile)){
        file.remove(outputfile)
    }
    return(df)
}


.getRuleSet3Scores_target <- function(sequences,
                                      tracrRNA=c("Hsu2013","Chen2013")
){
    tracrRNA <- match.arg(tracrRNA)
    stop("Target mode is not yet implemented.")
}








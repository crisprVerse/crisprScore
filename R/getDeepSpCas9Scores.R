#' @title Calculate on-target sgRNA activity scores for SpCas9 using DeepSpCas9
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas9-induced knockout using the DeepSpCas9 scoring method.
#' 
#' @param sequences Character vector of 30bp sequences needed for DeepSpCas9
#'     scoring, see details below.
#' @param fork Set to \code{TRUE} to preserve changes to the R
#'     configuration within the session.
#' 
#' @details The input sequences for DeepSpCas9 scoring require 4 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (20bp spacer sequence + 3bp PAM sequence ) and 3 nucleootides 
#'     downstream of the protospacer sequence, for a total of 30 nucleotides.
#' 
#' @return \strong{getDeepSpCas9Scores} returns a data.frame with
#'     \code{sequence} and \code{score} columns. The getDeepSpCas9Scores score
#'     takes on a value between 0 and 1. A higher score indicates higher
#'     knockout efficiency.
#' 
#' @references
#' Kim HK, Kim Y, Lee S, et al. SpCas9 activity prediction by DeepSpCas9, 
#'     a deep learning–base model with high generalization performance.
#'     Science advances. 2019 Nov 6;5(11):eaax9249.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' if (interactive()){
#' flank5 <- "ACCG" #4bp
#' spacer <- "AATCGATGCTGATGCTAGAT" #20bp
#' pam    <- "AGG" #3bp
#' flank3 <- "AAT" #3bp
#' input  <- paste0(flank5, spacer, pam, flank3) 
#' results <- getDeepSpCas9Scores(input)
#' }
#' @export 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getDeepSpCas9Scores <- function(sequences,
                                fork=FALSE){

    .dumpToFile <- function(sequences, file){
        col1 <- "Target number"
        col2 <- "30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)"
        out <- data.frame(col1=seq_along(sequences),
                          col2=sequences)
        colnames(out) <- c(col1, col2)
        write.table(out,
                    sep="\t",
                    file=file,
                    quote=FALSE,
                    col.names=TRUE,
                    row.names=FALSE)
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
    env <- basilisk::obtainEnvironmentPath(env_deepspcas9)
    envls <- basilisk.utils::activateEnvironment(env)
    on.exit(basilisk.utils::deactivateEnvironment(envls))
    programFile <- system.file("python",
                               "deepspcas9/getDeepSpCas9Scores.py",
                               package="crisprScore",
                               mustWork=TRUE)
    modelDir <- system.file("python",
                            "deepspcas9/DeepCas9_Final",
                            package="crisprScore")
    modelDir <- paste0(modelDir, "/")
    #cmd <- paste0("python ",
    #              programFile, " ",
    #              inputfile, " ",
    #              modelDir, " ",
    #              outputfile)
    if (sum(good)>0){
        .dumpToFile(sequences.valid, inputfile)
        pyBinary <- basilisk.utils:::getPythonBinary(env)

        system2(c(pyBinary,
                  programFile,
                  inputfile,
                  modelDir,
                  outputfile))
        scores <- read.table(outputfile)[,1]
        #scores <- read.table(outputfile)
        scores <- scores/100
        df$score[good] <- scores
    }
    if (file.exists(inputfile)){
        file.remove(inputfile)
    }
    if (file.exists(outputfile)){
        file.remove(outputfile)
    }
    #return(scores)
    return(df)
}







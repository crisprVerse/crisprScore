#' @title Calculate on-target sgRNA activity scores for Cas12a using DeepCpf1
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas12a-induced knockout using the DeepCpf1 scoring method.
#'     Currently not supported on Windows machines.
#' 
#' @param sequences Character vector of 34bp sequences needed for DeepCpf1
#'     scoring, see details below.
#' @param convertPAM Should non-canonical PAM sequences be converted to
#'     TTTC? TRUE by default. 
#' 
#' @details The input sequences for DeepCpf1 scoring require 4 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (4bp PAM sequence + 23bp spacer sequence) and 3 nucleootides 
#'     downstream of the protospacer sequence, for a total of 34 nucleotides.
#'     If \code{convertPAM} is set to \code{TRUE}, any non-canonical PAM
#'     sequence will be convert to TTTC for scoring purposes.
#' 
#' @return \strong{getDeepCpf1Scores} returns a data.frame with \code{sequence}
#'     and \code{score} columns. The DeepCpf1 score takes on a value between 0
#'     and 1. A higher score indicates higher knockout efficiency.
#' 
#' @references 
#' Kim, H., Min, S., Song, M. et al. Deep learning improves prediction of
#'     CRISPR–Cpf1 guide RNA activity. Nat Biotechnol 36, 239–241 (2018).
#'     \url{https://doi.org/10.1038/nbt.4061}.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' if (interactive()){
#' flank5 <- "ACC" #3bp
#' pam    <- "TTTT" #4bp
#' spacer <- "AATCGATGCTGATGCTAGATATT" #23bp
#' flank3 <- "AAGT" #4bp
#' input  <- paste0(flank5, pam, spacer, flank3) 
#' results <- getDeepCpf1Scores(input)
#' }
#' 
#' @inheritParams getAzimuthScores
#' @export 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom basilisk.utils activateEnvironment
#' @importFrom basilisk.utils deactivateEnvironment
getDeepCpf1Scores <- function(sequences,
                              convertPAM=TRUE,
                              fork=FALSE
){
    
    if (.Platform$OS.type=="windows"){
        stop("DeepCpf1 is not available for Windows at the moment.")
    }
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=34){
        stop("Provided sequences must have length 34nt",
             " ([4nt][TTTV][23mer][3nt]).")
    }
    if (convertPAM){
        pams <- substr(sequences, 5,8)
        wh <- which(!pams %in% c("TTTC", "TTTG", "TTTA"))
        if (length(wh)>0){
            sequences[wh] <- vapply(sequences[wh], function(x){
                paste0(substr(x,1,4), "TTTC",substr(x, 9,34), collapse="")
            }, FUN.VALUE="character")
        }
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
    env <- basilisk::obtainEnvironmentPath(env_deepcpf1)
    envls <- basilisk.utils::activateEnvironment(env)
    on.exit(basilisk.utils::deactivateEnvironment(envls))
    programFile <- system.file("python",
                               "deepcpf1/getDeepCpf1.py",
                               package="crisprScore",
                               mustWork=TRUE)
    if (sum(good)>0){
        .dumpToFile(sequences.valid, inputfile)
        
        pyBinary <- basilisk.utils:::getPythonBinary(env)
        system2(c(pyBinary,
                  programFile,
                  inputfile,
                  outputfile))
        scores <- read.table(outputfile)[,1]
        scores <- scores/100
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

.dumpToFile <- function(sequences, file){
    write.table(sequences,
              file=file,
              quote=FALSE,
              col.names=FALSE,
              row.names=FALSE)
}





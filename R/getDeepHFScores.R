#' @title Calculate on-target sgRNA activity scores for Cas9 using DeepHF
#' @description Calculate on-target sgRNA activity scores for 
#'     CRISPR/Cas9-induced knockout using the DeepHF scoring method. Both U6
#'     and T7 promoters are supported. Three different versions of the SpCas9
#'     nuclease are supported: wildtype (WT-SpCas9), high-fidelity Cas9
#'     (SpCas9-HF1) and enhanced Cas9 (eSpCas9).  
#'     
#' @param sequences Character vector of 23bp protospacer sequences.
#' @param enzyme Character string specifying the Cas9 variant.
#'     Wildtype Cas9 (WT) by default, see details below. 
#' @param promoter Character string speciyfing promoter used for expressing 
#'     sgRNAs for wildtype Cas9 (must be either "U6" or "T7").
#'     "U6" by default. 
#' 
#' @details Input sequences for DeepHF scoring must be 23bpprotospacer
#'     sequences (20bp spacer sequences + 3bp PAM sequences).
#'     Only canonical PAM sequences (NGG) are allowed.
#'     Users can specify for which Cas9 they wish to score sgRNAs
#'     by using the argument \code{enzyme}: "WT" for Wildtype Cas9 (WT-SpCas9),
#'     "HF" for high-fidelity Cas9 (SpCas9-HF), or "ESP" for enhanced 
#'     Cas9 (eSpCas9). For wildtype Cas9, users can also specify the promoter 
#'     used for expressing sgRNAs using the argument \code{promoter} ("U6" by
#'     default).
#' 
#' @return \strong{getDeepHFScores} returns a data.frame with \code{sequence}
#'     and \code{score} columns. The DeepHF score takes on a value
#'     between 0 and 1. A higher score indicates higher knockout efficiency.
#' 
#' @references 
#' Wang, D., Zhang, C., Wang, B. et al. Optimized CRISPR guide RNA design for
#'     two high-fidelity Cas9 variants by deep learning.
#'     Nat Commun 10, 4284 (2019).
#'     \url{https://doi.org/10.1038/s41467-019-12281-8}
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' \donttest{
#' spacer  <- "ATCGATGCTGATGCTAGATA" #20bp
#' pam     <- "AGG" #3bp 
#' input   <- paste0(spacer, pam) 
#' 
#' # Wiltype Cas9 using U6 promoter:
#' results <- getDeepHFScores(input)
#' 
#' # Wiltype Cas9 using T7 promoter:
#' results <- getDeepHFScores(input, promoter="T7")
#' 
#' #' High-fidelity Cas9:
#' results <- getDeepHFScores(input, enzyme="HF")
#' 
#' #' Enhanced Cas9:
#' results <- getDeepHFScores(input, enzyme="ESP")
#' }
#' 
#' @inheritParams getAzimuthScores
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getDeepHFScores <- function(sequences,
                            enzyme=c("WT", "ESP", "HF"),
                            promoter=c("U6", "T7"),
                            fork=FALSE){
    enzyme   <- match.arg(enzyme)
    promoter <- match.arg(promoter)
    if (enzyme=="WT" & promoter=="U6"){
        model_type <- "wt_u6"
        model_file <- "DeepWt_U6.hdf5"
    } else if (enzyme=="WT" & promoter=="T7"){
        model_type <- "wt_t7"
        model_file <- "DeepWt_T7.hdf5"
    } else if (enzyme=="ESP"){
        model_type <- "esp"
        model_file <- "esp_rnn_model.hdf5"
    } else if (enzyme=="HF"){
        model_type <- "hf"
        model_file <- "hf_rnn_model.hdf5"
    } 
    model_dir <- "/Users/fortinj2/crisprScoreData/inst/scripts/out/DeepHF"
    model_file <- file.path(model_dir, model_file)
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=23){
        stop("Provided sequences must have length 23nt (20nt-spacer + PAM)")
    }
    pams  <- substr(sequences,22,23)
    valid <- pams=="GG"
    if (sum(valid)!=length(pams)){
        stop("Positions 22 and 23 of the sequences must be G nucleotides",
             " (canonical PAM sequences).")
    }
    results <- basiliskRun(env=env_deephf,
                           shared=FALSE,
                           fork=fork,
                           fun=.deephf_python, 
                           sequences=sequences,
                           model_type=model_type,
                           model_file=model_file)
    return(results)
}


#' @importFrom reticulate import_from_path np_array
.deephf_python <- function(sequences,
                           model_type,
                           model_file){ 

    dir <- system.file("python",
                       "deephf",
                       package="crisprScore",
                       mustWork=TRUE)
    deephf <- import_from_path("getDeepHF", path=dir)

    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]
    if (length(sequences.valid)>0){
        scores <- deephf$getDeepHF(np_array(sequences.valid),
                                   model_type,
                                   model_file)
        scores <- scores[,c(1,4,5,6),drop=FALSE]
        scores$index <- scores$index+1
        scores <- scores[order(scores$index),,drop=FALSE]
        scores$index <- NULL
        colnames(scores) <- c("PAM", "sgrna.20mer", "score")
        df$score[good] <- scores$score
    }
    return(df)
}




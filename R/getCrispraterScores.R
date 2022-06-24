#' @title Calculate on-target sgRNA activity scores for Cas9 using CRISPRater
#' @description Calculate on-target sgRNA activity scores for 
#'     CRISPR/Cas9-induced knockout using the DeepHF scoring method. Both U6
#'     and T7 promoters are supported. Three different versions of the SpCas9
#'     nuclease are supported: wildtype (WT-SpCas9), high-fidelity Cas9
#'     (SpCas9-HF1) and enhanced Cas9 (eSpCas9). Currently not supported
#'     on Windows machines.
#'     
#' @param sequences Character vector of 20bp protospacer sequences.
#' 
#' @details Input sequences for CRISPRater scoring must be 20 spacer sequences.
#' 
#' @return \strong{getCrisprRaterScores} returns a data.frame with
#'     \code{sequence} and \code{score} columns. The CRISPRater score takes on
#'     a value between 0 and 1. A higher score indicates higher knockout 
#'     efficiency.
#' 
#' 
#' @references
#' Labuhn M, Adams FF, Ng M, et al. Refined sgRNA efficacy prediction improves
#' large-and small-scale CRISPR–Cas9 applications. Nucleic acids research. 2018
#' Feb 16;46(3):1375-85.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' spacer  <- "ATCGATGCTGATGCTAGATA" #20bp
#' results <- getCRISPRaterScores(spacer)
#' 
#' @export
getCRISPRaterScores <- function(sequences){
    
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=20){
        stop("Provided sequences must have length 20nt (20nt-spacer")
    }
    results <- .getModelScoreCRISPRater(sequences)
    return(results)
}



# Must be 20nt long
#spacers <- c("GGTGCTGATGCTGTGTGATG",
#             "GGTGCTGATAAAGTGTGATG")
.getModelScoreCRISPRater <- function(sequences){

    model_weights <- c(0.14177385,
                       0.06966514,
                       0.04216254,
                       0.03303432,
                       0.02355430,
                       -0.04746424,
                       -0.04878001,
                       -0.06981921,
                       -0.07087756,
                       -0.08160700)
    model_offset <- 0.6505037

    features <- .extractFeaturesForCRISPRater(sequences)
    score <- model_offset + model_weights %*% t(features)
    return(score)
}




#' @importFrom Biostrings DNAStringSet letterFrequency
.extractFeaturesForCRISPRater <- function(sequences){
    gc <- rowSums(letterFrequency(DNAStringSet(substr(sequences,4,14)),
                                  c("G", "C")))/10
    sequences <- DNAStringSet(sequences)
    mat <- as.matrix(DNAStringSet(sequences))
    features <- list()
    features[[1]] <- gc
    features[[2]] <- mat[,20,drop=FALSE]=="G"
    features[[3]] <- mat[,3,drop=FALSE]=="T" | mat[3]=="A"
    features[[4]] <- mat[,12,drop=FALSE]=="G"| mat[12]=="A"
    features[[5]] <- mat[,6,drop=FALSE]=="G"
    features[[6]] <- mat[,4,drop=FALSE]=="T" | mat[4]=="A"
    features[[7]] <- mat[,18,drop=FALSE]=="G" | mat[18]=="A"
    features[[8]] <- mat[,5,drop=FALSE]=="C" | mat[5]=="A" 
    features[[9]] <- mat[,14,drop=FALSE]=="G"
    features[[10]] <- mat[,15,drop=FALSE]=="A"
    features <- do.call(cbind,features)
    return(features)
}

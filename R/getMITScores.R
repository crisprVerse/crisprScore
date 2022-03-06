#' @title Calculate MIT off-target specificity scores for CRISPR/Cas9
#' @description Calculate MIT off-target specificity scores for CRISPR/Cas9.
#' 
#' @param spacers Character vector of 20bp spacer sequences.
#' @param protospacers Character vector of 20bp protospacer sequences 
#'     for off-targets.
#' @param pams Character vector of 3nt PAM sequences. 
#' @param includeDistance Should distance between mismatches be considered
#'     during scoring? TRUE by default.
#' 
#' @return \strong{getMITScores} returns a data.frame with \code{spacer},
#'     \code{protospacer}, and \code{score} columns. The MIT score takes on a
#'     value between 0 and 1. For a given pair (on-target, off-target), a 
#'     higher MIT score indicates a higher likelihood for the Cas9 nuclease to
#'     cut at the off-target. Non-canonical PAM sequences are taken into account
#'     by the MIT algorithm.
#'     
#' @references Hsu, P., Scott, D., Weinstein, J. et al. DNA targeting
#'     specificity of RNA-guided Cas9 nucleases.
#'     Nat Biotechnol 31, 827–832 (2013).
#'     \url{https://doi.org/10.1038/nbt.2647}.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' # Calculating MIT scores for two off-targets with respect to
#' # one spacer sequence:
#' spacer <- "AGGTGTAGTGTGTGTGATAA"
#' protospacer1 <- "CGGTGTAGTGTGTGTGATAA"
#' protospacer2 <- "CGGTGTCGTGTGTGTGATAA"
#' pams <- c("AGG", "CGG")
#' results <- getMITScores(spacers=spacer,
#'     protospacers=c(protospacer1, protospacer2),
#'     pams=pams
#' )
#' 
#' @importFrom Biostrings DNAStringSet
#' @export
getMITScores <- function(spacers,
                         protospacers,
                         pams,
                         includeDistance=TRUE
){
    spacers       <- .checkSequenceInputs(spacers)
    protospacers  <- .checkSequenceInputs(protospacers)
    if (unique(nchar(protospacers))!=20){
        stop("Protospacer sequences must have length 20nt.")
    } 
    if (unique(nchar(spacers))!=20){
        stop("Spacer sequences must have length 20nt.")
    }
    if (unique(nchar(pams))!=3){
        stop("PAM sequences must have length 3nt.")
    }
    if (length(spacers)==1){
        spacers <- rep(spacers, length(protospacers))
    } else {
        if (length(spacers) != length(protospacers)){
            stop("Input vectors 'spacers' and 'protospacers' must",
                 " have the same length.")
        }
    }
    pams <- substr(pams, 2,3)
    
    spacers.wt  <- spacers
    spacers.off <- protospacers 
    spacers.wt   <- DNAStringSet(spacers.wt)
    spacers.off  <- DNAStringSet(spacers.off)
    x <- as.matrix(spacers.wt)
    y <- as.matrix(spacers.off)
    wh <- x!=y
    mms <- lapply(seq_len(nrow(wh)), function(i){
        which(wh[i,]) 
    })
    mm.scores  <- vapply(mms,
                         .calculateSpacerSequenceScore,
                         includeDistance=includeDistance,
                         FUN.VALUE=0)
    pam.scores <- cfd.pam.weights.cas9[pams]
    mit.scores <- mm.scores*pam.scores
    mit.scores[is.na(mit.scores)] <- 1
    mit.scores <- as.numeric(mit.scores)
    data.frame(spacer=spacers, 
               protospacer=protospacers,
               score=mit.scores)
}


.calculateSpacerSequenceScore <- function(indices,
                                          includeDistance=TRUE
){
    if (length(indices)==0){
        score <- 1
    } else {
        m <- length(indices)
        if (m==1){
            d <- 19
        } else {
            d <- (max(indices)-min(indices))/(m-1)
        }
        t1 <- prod(mit.weights[as.character(indices)])
        t2 <- 1/(m^2)
        t3 <- 1/((19-d)/19*4+1)
        if (includeDistance){
            score <- t1*t2*t3
        } else {
            score <- t1*t2
        }
    }
    return(score)
}



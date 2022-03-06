#' @title Calculate CFD off-target specificity scores 
#' @description Calculate cutting frequency determination (CFD) off-target
#'     specificity scores for CRISPR/Cas9 or CRISPR/CasRX.
#' 
#' @param spacers Character vector of 20bp spacer sequences.
#'     Must be in 5' to 3' direction.
#'     For SpCas9, must be of length 20bp.
#'     For CasRx, must be at most of length 27bp.
#' @param protospacers Character vector of 20bp protospacer sequences
#'     (target sequences). Must be in 5' to 3' direction.
#' @param pams Character vector of PAM sequences.
#' @param nuclease String specifying the nuclease. Either "SpCas9" (default)
#'     or "CasRx".
#' 
#' @return \strong{getCFDScores} returns a data.frame with \code{spacer},
#'     \code{protospacer}, and \code{score} columns. The CFD score takes
#'     on a value between 0 and 1. For a given pair (on-target, off-target), 
#'     a higher CFD score indicates a higher likelihood for the nuclease
#'     to cut at the off-target. Non-canonical PAM sequences are taken into
#'     account by the CFD algorithm.
#' 
#' @references
#' Doench, J., Fusi, N., Sullender, M. et al. Optimized sgRNA design to
#'     maximize activity and minimize off-target effects of CRISPR-Cas9.
#'     Nat Biotechnol 34, 184–191 (2016).
#'     \url{https://doi.org/10.1038/nbt.3437}.
#' 
#' @author Jean-Philippe Fortin
#'
#' 
#' @examples
#' # Calculating MIT scores for two off-targets with respect to
#' # one spacer sequence:
#' spacer <- "AGGTGTAGTGTGTGTGATAA"
#' protospacer1 <- "CGGTGTAGTGTGTGTGATAA"
#' protospacer2 <- "CGGTGTCGTGTGTGTGATAA"
#' results <- getCFDScores(spacers=spacer,
#'     protospacers=c(protospacer1, protospacer2),
#'     pams=c("AGG", "CGG")
#' )
#' 
#' @importFrom Biostrings width
#' @export
getCFDScores <- function(spacers,
                         protospacers,
                         pams,
                         nuclease=c("SpCas9", "CasRx")
){
    nuclease <- match.arg(nuclease)
    spacers       <- .checkSequenceInputs(spacers)
    protospacers  <- .checkSequenceInputs(protospacers)
    if (length(spacers)==1){
        spacers <- rep(spacers, length(protospacers))
    } else {
        if (length(spacers) != length(protospacers)){
            stop("Input vectors 'spacers' and 'protospacers' must", 
                 " have the same length.")
        }
    }

    if (nuclease=="SpCas9"){
        if (unique(nchar(protospacers))!=20){
            stop("Protospacer sequences must have length 20nt.")
        } 
        if (unique(nchar(spacers))!=20){
            stop("Spacer sequences must have length 20nt.")
        }
        if (unique(nchar(pams))!=3){
            stop("PAM sequences must have length 3nt.")
        }
        pams <- substr(pams,2,3)
    } else if (nuclease=="CasRx"){
        if (unique(nchar(protospacers))!=27){
            stop("Protospacer sequences must have length 27nt.")
        } 
        if (unique(nchar(spacers))!=27){
            stop("Spacer sequences must have length 27nt.")
        }
        if (unique(nchar(pams))!=1){
            stop("PAM sequences must have length 1nt.")
        }
    }
   
    spacers.wt  <- spacers
    spacers.off <- protospacers
    spacers.wt   <- DNAStringSet(spacers.wt)
    spacers.off  <- DNAStringSet(spacers.off)
    x <- as.matrix(spacers.wt)
    y <- as.matrix(spacers.off)
    wh <- x!=y
    spacer.len <- unique(width(spacers.wt))
    MM <- matrix(rep(seq_len(spacer.len), each=nrow(x)),
                 nrow=nrow(x), ncol=ncol(x))
    MM[wh] <- paste0(x[wh], y[wh], MM[wh])
    mm <- lapply(seq_len(nrow(MM)), function(i){
        x <- MM[i,]
        xx <- x[grepl("[ACGT]+",x)]
        if (length(xx)==0){
            xx <- NA
        }
        return(xx)
    })

    pam.weights <- .getPamWeights(nuclease)
    mm.weights  <- .getMmWeights(nuclease)
    mm.scores <- vapply(mm,
                        function(x) prod(mm.weights[x]),
                        FUN.VALUE=0)
    mm.scores[is.na(mm.scores)] <- 1
    pam.scores <- pam.weights[pams]
    cfd.scores <- as.numeric(mm.scores*pam.scores)
    data.frame(spacer=spacers, 
               protospacer=protospacers,
               score=cfd.scores)
}


.getPamWeights <- function(nuclease){
    if (nuclease=="SpCas9"){
        ws <- cfd.pam.weights.cas9
    } else if (nuclease=="CasRx"){
        ws <- cfd.pam.weights.casrx
    }
    return(ws)
}

.getMmWeights <- function(nuclease){
    if (nuclease=="SpCas9"){
        ws <- cfd.mm.weights.cas9
    } else if (nuclease=="CasRx"){
        ws <- cfd.mm.weights.casrx
    }
    return(ws)
}






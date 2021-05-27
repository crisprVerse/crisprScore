#' @title Calculate CFD off-target specificity scores for CRISPR/Cas9
#' @description Calculate cutting frequency determination (CFD) off-target
#'     specificity scores for CRISPR/Cas9.
#' 
#' @param spacers Character vector of 20bp spacer sequences.
#' @param protospacers Character vector of 23bp protospacer sequences. 
#' 
#' @return \strong{getCFDScores} returns a data.frame with \code{spacer},
#'     \code{protospacer}, and \code{score} columns. The CFD score takes
#'     on a value between 0 and 1. For a given pair (on-target, off-target), 
#'     a higher CFD score indicates a higher likelihood for the Cas9 nuclease
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
#' protospacer1 <- paste0("CGGTGTAGTGTGTGTGATAA", "AGG")
#' protospacer2 <- paste0("CGGTGTCGTGTGTGTGATAA", "CGG")
#' results <- getCFDScores(spacers=spacer,
#'     protospacers=c(protospacer1, protospacer2)
#' )
#' 
#' @importFrom Biostrings width
#' @export
getCFDScores <- function(spacers, protospacers){
    spacers       <- .checkSequenceInputs(spacers)
    protospacers  <- .checkSequenceInputs(protospacers)
    if (unique(nchar(protospacers))!=23){
        stop("Protospacer sequences must have length 23nt (20nt-spacer + PAM).")
    } 
    if (unique(nchar(spacers))!=20){
        stop("Spacer sequences must have length 20nt.")
    }
    if (length(spacers)==1){
        spacers <- rep(spacers, length(protospacers))
    } else {
        if (length(spacers) != length(protospacers)){
            stop("Input vectors 'spacers' and 'protospacers' must", 
                 " have the same length.")
        }
    }
    spacers.wt  <- spacers
    spacers.off <- substr(protospacers, 1,nchar(protospacers) - 3) 
    spacers.wt   <- DNAStringSet(spacers.wt)
    spacers.off  <- DNAStringSet(spacers.off)
    pams <- substr(protospacers, nchar(protospacers) - 1, nchar(protospacers))
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
    mm.scores <- vapply(mm,
                        function(x) prod(cfd.mm.weights.cas9[x]),
                        FUN.VALUE=0)
    mm.scores[is.na(mm.scores)] <- 1
    pam.scores <- cfd.pam.weights.cas9[pams]
    cfd.scores <- as.numeric(mm.scores*pam.scores)
    data.frame(spacer=spacers, 
               protospacer=protospacers,
               score=cfd.scores)
}



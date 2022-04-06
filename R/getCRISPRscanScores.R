#' @title Calculate on-target sgRNA activity scores for Cas9 using CRISPRscan
#' 
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas9-induced knockout using the CRISPRscan scoring method.
#'     The method is also know as the Moreno-Mateos score. 
#'     The CRISPRscan  algorithm was trained using in vitro transcription
#'     of sgRNAs using a T7 promoter, and might therefore be 
#'     suboptimal to predict sgRNA activity when expressed from U6 promoter.
#' 
#' @param sequences Character vector of 35bp sequences needed for
#'     CRISPRscan scoring, see details below.
#' 
#' @details The input sequences for Rule Set 1 scoring require 6 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (23 nucleotides) and 6 nucleootides downstream of the protospacer
#'     sequence, for a total of 35 nucleotides: [6nt][20nt-spacer][NGG][6nt].
#'     Note that a canonical PAM sequence (NGG) is required for CRISPRscan.
#' 
#' @return \strong{getCRISPRscanScores} returns a data.frame with \code{sequence} 
#'     and \code{score} columns. The CRISPRscan score takes on a value between 0
#'     and 1. A higher score indicates higher knockout efficiency.
#' 
#' 
#' @references 
#' Moreno-Mateos MA, et al. CRISPRscan: designing highly efficient
#'     sgRNAs for CRISPR-Cas9 targeting in vivo. Nature methods. 2015 Oct;12(10):982-8.
#' 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' 
#' flank5 <- "ACCTGG" #6bp
#' spacer <- "ATCGATGCTGATGCTAGATA" #20bp
#' pam    <- "AGG" #3bp 
#' flank3 <- "TTGAGC" #6bp
#' input  <- paste0(flank5, spacer, pam, flank3) 
#' results <- getCRISPRscanScores(input)
#' 
#' @export
getCRISPRscanScores <- function(sequences){
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=35){
        stop("Provided sequences must have length 30nt ",
             "([6nt][20nt-spacer][PAM][6nt]).")
    }
    pams  <- substr(sequences,28,29)
    valid <- pams=="GG"
    if (sum(valid)!=length(pams)){
        stop("Positions 28 and 29 of the sequences must be G",
             " nucleotides (canonical PAM sequences required).")
    }
    results <- .crisprscan(sequences)
    return(results)
}


.crisprscan <- function(sequences){

    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]
    ns <- length(sequences.valid)
    if (ns>0){
        scores <- .crisprscan_engine(sequences.valid)
        df$score[good] <- scores
    }
    return(df)
}




.crisprscan_engine <- function(seqs){
    
    nseq <- length(seqs)

    .getWeights <- function(){
        tempWs <- read.csv(system.file("crisprscan",
                                       "crisprscan_coefficients.csv",
                                       package="crisprScore"))
        ws <- tempWs[,2]
        names(ws) <- tempWs[,1]
        return(ws)
    }

    .getFeaturesScore <- function(seqs){
        miscFeatures <- c("Intercept")
        nucFeatures <- setdiff(names(ws), miscFeatures)
        motifs <- gsub("[0-9]+", "", nucFeatures)
        start  <- as.integer(gsub("[ACGT]+", "", nucFeatures))
        end    <- start + nchar(motifs) - 1
        roster <- data.frame(start=start,
                             end=end,
                             motif=motifs,
                             w=ws[nucFeatures])
        roster <- roster[order(roster$start),]

        nfeatures <- nrow(roster)
        out <- lapply(seq_len(nseq), function(j){
            vapply(seq_len(nfeatures), function(i){
                x <- substr(seqs[j],
                            roster$start[i],
                            roster$end[i])
                as.numeric(x==roster$motif[i])
            }, FUN.VALUE=1)
        })
        out <- do.call(rbind, out)
        out <- out %*% as.numeric(roster$w)
        return(out)
    }

    ws <- .getWeights()
    score <- ws["Intercept"]
    score <- score + .getFeaturesScore(seqs)
    return(score)
}



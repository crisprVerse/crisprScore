#' @title Calculate on-target sgRNA activity scores for Cas9 using Rule Set 1
#' 
#' @description Calculate on-target sgRNA activity scores for
#'     CRISPR/Cas9-induced knockout using the Rule Set 1 scoring method.
#'     The Rule Set 1 algorithm was an early on-target efficiency method
#'     developed by the Doench lab. 
#' 
#' @param sequences Character vector of 30bp sequences needed for
#'     Rule Set 1 scoring, see details below.
#' 
#' @details The input sequences for Rule Set 1 scoring require 4 nucleotides
#'     upstream of the protospacer sequence, the protospacer sequence
#'     itself (23 nucleotides) and 3 nucleootides downstream of the protospacer
#'     sequence, for a total of 30 nucleotides: [4nt][20nt-spacer][NGG][3nt].
#'     Note that a canonical PAM sequence (NGG) is required for Rule Set 1. 
#' 
#' @return \strong{getRuleSet1Scores} returns a data.frame with \code{sequence} 
#'     and \code{score} columns. The Rule Set 1 score takes on a value between 0
#'     and 1. A higher score indicates higher knockout efficiency.
#' 
#' 
#' @references 
#' Doench, John G., et al. Rational design of highly active sgRNAs
#'     for CRISPR-Cas9–mediated gene inactivation. Nat Biotech 32, 1262-1267 (2014).
#'     \url{https://doi.org/10.1038/nbt.3026}.
#' 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' 
#' flank5 <- "ACCT" #4bp
#' spacer <- "ATCGATGCTGATGCTAGATA" #20bp
#' pam    <- "AGG" #3bp 
#' flank3 <- "TTG" #3bp
#' input  <- paste0(flank5, spacer, pam, flank3) 
#' results <- getRuleSet1Scores(input)
#' 
#' @export
getRuleSet1Scores <- function(sequences){
    sequences <- .checkSequenceInputs(sequences)
    if (unique(nchar(sequences))!=30){
        stop("Provided sequences must have length 30nt ",
             "([4nt][20nt-spacer][PAM][3nt]).")
    }
    pams  <- substr(sequences,26,27)
    valid <- pams=="GG"
    if (sum(valid)!=length(pams)){
        stop("Positions 26 and 27 of the sequences must be G",
             " nucleotides (canonical PAM sequences required).")
    }
    results <- .ruleset1(sequences)
    return(results)
}


.ruleset1 <- function(sequences){

    df <- data.frame(sequence=sequences,
                     score=NA_real_,
                     stringsAsFactors=FALSE)
    good <- !grepl("N", sequences)
    sequences.valid <- sequences[good]
    ns <- length(sequences.valid)
    if (ns>0){
        scores <- .ruleset1_engine(sequences.valid)
        df$score[good] <- scores
    }
    return(df)
}




#' @importFrom Biostrings DNAStringSet letterFrequency
.ruleset1_engine <- function(seqs){
    
    spacerLen <- 20
    offset <- 4
    nseq <- length(seqs)

    .getWeights <- function(){
        tempWs <- read.csv(system.file("ruleset1",
                                       "coefficients.csv",
                                       package="crisprScore"))
        ws <- tempWs[,2]
        names(ws) <- tempWs[,1]
        return(ws)
    }

    .getGCContentScore <- function(seqs){
        spacers <- substr(seqs,
                          offset + 1,
                          offset + spacerLen)
        spacers <- DNAStringSet(spacers)
        gc <- rowSums(letterFrequency(spacers, c("G", "C")))
        score1  <- as.numeric(gc<10) * (10 - gc) * ws["gc_low"]
        score2  <- as.numeric(gc>10) * (gc - 10) * ws["gc_high"]
        return(score1+score2)
    }

    .getFeaturesScore <- function(seqs){

        nucFeatures <- setdiff(names(ws),
                               c("Intercept", "gc_low", "gc_high"))
        motifs <- gsub("[0-9]+", "", nucFeatures)
        start  <- as.numeric(gsub("[ACGT]+", "", nucFeatures))
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
    score <- score + .getGCContentScore(seqs)
    score <- score + .getFeaturesScore(seqs)
    score <- 1/(1 + exp(-score))
    return(score)
}




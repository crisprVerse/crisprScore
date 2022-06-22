

# Must be 20nt long
#spacers <- c("GGTGCTGATGCTGTGTGATG",
#             "GGTGCTGATAAAGTGTGATG")
.getModelScore <- function(spacers){
    features <- .extractFeatures(spacers)
    score <- model_offset + model_weights%*%t(features)
    return(score)
}

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


#' @importFrom Biostrings DNAStringSet letterFrequency
.extractFeatures <- function(spacers){
    gc <- rowSums(letterFrequency(DNAStringSet(substr(spacers,4,14)),
                                  c("G", "C")))/10
    spacers <- DNAStringSet(spacers)
    mat <- as.matrix(DNAStringSet(spacers))
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

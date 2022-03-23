#' @title Calculate on-target sgRNA activity scores for CasRx using CasRx-RF
#' 
#' @description Calculate on-target sgRNA activity scores 
#'     for CasRx (RfxCas13d) using the CasRx-RF algorithm.
#' 
#' @param mrnaSequence A \code{DNAStringSet} represeting the mRNA sequence
#'     for which to extract spacer sequences and calculate scores.
#' @param directRepeat String specifying the direct repeat used in the 
#'     CasRx construct.
#' @param binaries Named list of paths for binaries needed for 
#'     CasRx-RF. Names of the list must be "RNAfold", "RNAhybrid",
#'     and "RNAplfold". Each list element is a string specifying
#'     the path of the binary. If NULL (default), binaries must be
#'     available on the PATH.
#' @param sort Should spacers be sorted by score? FALSE by default.
#' @param verbose Should messages be printed to console? TRUE by default.
#' 
#' @details The function first extracts all 23mer spacer sequences targeting
#'    the mRNA sequence, and scores them for on-target activity. 
#'    
#' @return A data.frame with the following columns:
#' \itemize{
#' \item \code{ID} Character vector specifying spacer ID.
#' \item \code{spacer} 23-mer spacer sequence.
#' \item \code{pfs_site} coordinate of the protospacer flanking sequence (PFS).
#' \item \code{protospacer} 23-mer protospacer sequence (reverse complement of the spacer sequence).
#' \item \code{PFS} PFS nucleotide.
#' \item \code{score} Raw score (not standardized).
#' \item \code{standardizedScore} Score standardized between 0 and 1.
#' \item \code{quartile} Quartile score (1 to 4, with 4 being the best quartile.
#' }
#' 
#' A scores closer to 1 indicates higher predicted on-target activity.
#' 
#' @references 
#' Wessels HH, Méndez-Mancilla A, Guo X, et al.
#' Massively parallel Cas13 screens reveal principles 
#' for guide RNA design. Nat biotechnol. 2020 Jun;38(6):722-7.
#' 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' 
#' \donttest{
#' fasta <- file.path(system.file(package="crisprScore"),
#' "casrxrf/test.fa"
#' mrnaSequence <- Biostrings::readDNAStringSet(filepath=fasta,
#'     format="fasta",
#'     use.names=TRUE)
#' results <- getCasRxRFScores(mrnaSequence)
#' }
#' 
#' @export
getCasRxRFScores <- function(mrnaSequence,
                             directRepeat="aacccctaccaactggtcggggtttgaaac",
                             binaries=NULL,
                             sort=FALSE,
                             verbose=TRUE
){
    if (!is(mrnaSequence, "DNAStringSet")){
        stop("mrnaSequence must be a DNAStringSet")
    }
    if (length(mrnaSequence)>1){
        stop("mrnaSequence must be of length 1.")
    }
    spacers <- .prepareCasRxSpacers(mrnaSequence)
  
    if (is.null(binaries)){
         binaries <- list(RNAfold="RNAfold",
                          RNAplfold="RNAplfold",
                          RNAhybrid="RNAhybrid")
    } else {
        requiredNames <- c("RNAfold", "RNAhybrid", "RNAplfold")
        if (!all.equal(sort(names(binaries)), requiredNames)){
            stop("binaries must be either NULL or a named list ",
                 "with the following names: RNAfold, RNAhybrid, and RNAplfold.")
        }
    }
    spacers <- .addSpacerFeatures(spacers=spacers,
                                  directRepeat=directRepeat,
                                  mrnaSequence=mrnaSequence,
                                  binaries=binaries,
                                  verbose=verbose)
    if (verbose){
        cat("[getCasRxRFScores] Calculating scores \n")
    }
    spacers <- .addRawCasRxScores(df=spacers, 
                                  sort=sort)
    spacers <- .addStandardizedCasRxScores(df=spacers)
    spacers <- .cleanColumns(spacers)
    return(spacers)
}
#fasta <- "/Users/fortinj2/Cas13design/inst/fasta/test.fa" # Argument 1
#mrnaSequence <- Biostrings::readDNAStringSet(filepath=fasta,
                                 #format="fasta",
                                 #use.names=TRUE)
#spacers <- getCasRxRFScores(mrnaSequence)



#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges Views
#' @importFrom stringr str_pad
#' @importFrom BiocGenerics start end 
.prepareCasRxSpacers <- function(mrnaSequence,
                                 spacerLength=23,
                                 GC.min=0.1,
                                 GC.max=0.9
){    
    spacers <- .extractSpacersFromDNAString(mrnaSequence[[1]],
                                            spacerLength=spacerLength)
    spacers <- .filterSpacers(spacers,
                              GC.min=GC.min,
                              GC.max=GC.max,
                              HomopolymerLengthT=4,
                              HomopolymerLengthnonT=5)
    spacers <- reverseComplement(spacers)
    spacerSet <- DNAStringSet(spacers) 

    # Assign names
    spacerIndex <- str_pad(seq(1,length(spacers),1),
                           nchar(length(spacers)),
                           pad = "0")
    names(spacerSet) <- paste0("crRNA", spacerIndex,":",
                               width(mrnaSequence)-end(spacers)+1, "-" ,
                               width(mrnaSequence)-start(spacers)+1) 
    # Transform to data frame
    spacers <- as.data.frame(spacerSet)
    colnames(spacers) <- "spacer"
    spacers$pos <- vapply(rownames(spacers), .extractPosition, FUN.VALUE=0)

     # Adding pfs_site column:
    spacers$pfs_site <- spacers$pos+1

    # Adding protospacer column:
    protospacer <- reverseComplement(DNAStringSet(spacers$spacer))
    spacers$protospacer <- as.character(protospacer)

    # Adding PFS column
    dna <- as.character(as.vector(mrnaSequence[[1]]))
    spacers$PFS <- dna[spacers$pfs_site]

    # Filtering out spacers out of range:
    good <- spacers$pfs_site<=width(mrnaSequence)
    spacers <- spacers[good,,drop=FALSE]

    return(spacers)
}



.filterSpacers <- function(spacers,
                           GC.min=0.1,
                           GC.max=0.9,
                           HomopolymerLengthT=4,
                           HomopolymerLengthnonT=5
){
    GC <- vapply(spacers, .getGC, FUN.VALUE=0)
    consecutiveBases <- lapply(spacers, .getConsecutiveBases)
    consecutiveBases <- do.call(rbind, consecutiveBases)
    GC.index <- which(GC > GC.min & GC < GC.max)
    consecutiveBases.index <- which(apply(consecutiveBases[,c("T","C","G")],
                                          MARGIN=1,max) <= HomopolymerLengthnonT & consecutiveBases[,c("A")] <= HomopolymerLengthT)
  

    if (length(GC.index) > 0 & length(consecutiveBases.index) > 0){
        idx <- intersect(GC.index, consecutiveBases.index)
        if (length(idx) > 0){
            spacers <- return(spacers[idx])
        } else{
            stop("no intersect: GC and Homopolymers out of range")
        }
    } else{
        if ( length(GC.index) == 0 & length(consecutiveBases.index) > 0){
            stop("GC out of range")
        } else if (  length(GC.index) > 0 & length(consecutiveBases.index) == 0 ){
            stop("Homopolymers out of range")
        } else{
            stop("GC and Homopolymers out of range")
        }
    }  
    return(spacers)
}


.extractSpacersFromDNAString <- function(dnaString,
                                         spacerLength=23
){
    start <- seq(1, length(dnaString)-spacerLength+1, 1)
    end   <- seq(spacerLength, length(dnaString), 1)
    subStrings <- Views(dnaString,
                        start=start,
                        end=end)
    return(subStrings)
}


.cleanColumns <- function(spacers){

    spacers <- spacers[, -which( colnames(spacers) %in% c(.FEATURES,
                                                         c('NTdens_min_GC','pT','pG|pC','pA|pT','pTT')) )]
    spacers <- cbind.data.frame(rownames(spacers) ,spacers)
    colnames(spacers)[1] <- 'ID'
    rownames(spacers) <- NULL
    spacers$pos <- NULL
    spacers$rank <- NULL
    return(spacers)
}






#' @importFrom stats predict
.predictRawCasRxScores  <- function(x,
                                    fittedModelList
){
    Model           <- fittedModelList[[1]]
    ModelInputMeans <- fittedModelList[[2]]
    ModelInputSDs   <- fittedModelList[[3]]
  
    # Removing incomplete entries
    FIELDS <- .FEATURES[-1]
    x <- x[complete.cases(x), FIELDS]
  
    # Scaling features:
    cond1 <- sapply(x, is.numeric)
    cond2 <- colnames(x) != "normCS"
    cond3 <- colnames(x) != "Gquad"
    cond4 <- colnames(x) != "DR" 
    cond5 <- grepl("[A,C,G,T]_",colnames(x)) == FALSE 
    numeric <-  cond1 & cond2 & cond3 & cond4 & cond5
    x[,numeric] <- scale(x[,numeric],
                         center=ModelInputMeans,
                         scale=ModelInputSDs)
  
    x[,numeric][x[,numeric] > 1] <- 1
    x[,numeric][x[,numeric] < 0] <- 0
    return(predict(Model, newdata=x))
}

.addRawCasRxScores <- function(df,
                               fittedModelList=NULL,
                               sort=FALSE
){
    if (is.null(fittedModelList)){
        fittedModelList <- .getPretrainedModel()
    }
    scores <- .predictRawCasRxScores(df, fittedModelList)
    wh <- match(names(scores) , rownames(df))
    
    if (!all(names(scores)==rownames(df)[wh])){
        stop("Exiting! Guide names do not correspond.")
    }
    
    df$score <- NA
    df$rank  <- NA
    df$score[wh] <- scores
    
    # Assign Rank
    # For guides that reside to close to the target's 5' 
    # end it may be that not all features are assigned. 
    # Thus, all guides with NA features will not be ranked.
    # The rank ranges from 0 to 1, with 1 being the highest rank
    good <- !is.na(df$score)
    df$rank[good] <- signif(1- (rank(-df$score[good],
                                     na.last = TRUE,
                                     ties.method = "first")/length(df$score[wh])),4) 
    if (sort){
        df <- df[order(df$rank, decreasing=TRUE),,drop=FALSE]
    }
    return(df)
}

.addStandardizedCasRxScores <- function(df){
    trainingData <- .getTrainingData()
    # Standardize original Screen data to bin predicted
    # guide scores according to screen efficacy quartiles
    min <- quantile(trainingData$normCS, na.rm=TRUE , probs = 0.05)
    max <- quantile(trainingData$normCS, na.rm=TRUE , probs = 0.95)
    trainingData$standardizedScore <- vapply(trainingData$normCS,
                                             FUN=.standardizeScores,
                                             min=min,
                                             max=max,
                                             FUN.VALUE=0)

    # Standardize GuideScores relative to ModelInput
    df$standardizedScore <- vapply(df$score,
                                   FUN=.standardizeScores,
                                   min=min,
                                   max=max,
                                   FUN.VALUE=0)
    df$standardizedScore[df$standardizedScore > 1] <- 1
    df$standardizedScore[df$standardizedScore < 0] <- 0
    df[["quartile"]] <- .assignQuartiles(df$standardizedScore,
                                         q=quantile(trainingData$standardizedScore,
                                                    na.rm=TRUE))
    return(df)
}


.standardizeScores <-  function(x, min, max){
    if(is.na(x)){
      out <- NA
    } else{
      out <- (x - min) / (max - min)
    }
    return(out)
}





.assignQuartiles = function(x,q = quantile(x)){
    quartiles = vector()
    for ( i in 1:length(x)){
      if (is.na(x[i]) == T){
          quartiles[i] = NA
      } else if (x[i] <= q["25%"]){
          quartiles[i] = 1
      } else if (x[i] > q["25%"] & x[i] <= q["50%"]){
          quartiles[i] = 2
      } else if (x[i] > q["50%"] & x[i] <= q["75%"]){
          quartiles[i] = 3
      } else if (x[i] > q["75%"]){
          quartiles[i] = 4
      } else{
          stop("Exiting. Value outside quartiles")
      }
    }
    return(quartiles) 
}




########## Spacer Features
.addSpacerFeatures <- function(spacers,
                               directRepeat,
                               mrnaSequence,
                               binaries,
                               verbose=TRUE
){
    bin.RNAfold   <- binaries[["RNAfold"]]
    bin.RNAplfold <- binaries[["RNAplfold"]]
    bin.RNAhybrid <- binaries[["RNAhybrid"]]

    if (verbose){
        cat("[addCasRxScores] Calculating MFE features \n")
    }
    spacers <- .addMFEFeatures(spacers,
                               directRepeat=directRepeat,
                               bin.RNAfold=bin.RNAfold)
    
    if (verbose){
        cat("[addCasRxScores] Calculating accessibility features \n")
    }
    spacers <- .addTargetSiteAccessibilityFeatures(dat=spacers,
                                                   mrnaSequence=mrnaSequence,
                                                   bin.RNAplfold=bin.RNAplfold)

    if (verbose){
        cat("[addCasRxScores] Calculating hybridization features \n") 
    }
    spacers <- .addRNAhybridFeatures(spacers,
                                     bin.RNAhybrid=bin.RNAhybrid)
    
    if (verbose){
        cat("[addCasRxScores] Calculating nucleotide density features \n")
    }
    spacers <- .addNTpointdensities(DAT=spacers,
                                    mrnaSequence=mrnaSequence)
    spacers <- .addLetterProbs(x=spacers)
    return(spacers)
}


# fold nt including the direct repeat sequence
.addMFEFeatures <- function(spacers,
                            directRepeat=directRepeat,
                            bin.RNAfold
){
    stats <- lapply(spacers$spacer,
                    FUN=.getMFE,
                    DR=directRepeat,
                    bin.RNAfold=bin.RNAfold)
    stats <- do.call(rbind, stats)
    spacers[["MFE"]]   <- stats[,1]
    spacers[["DR"]]    <- stats[,2]
    spacers[["Gquad"]] <- stats[,3]
    return(spacers)
}

.getMFE = function(spacer,
                   DR,
                   bin.RNAfold
){
    crRNA <- paste0(DR, spacer)
    cmd <- paste0( "echo ",crRNA, " | ", bin.RNAfold , " --gquad --noPS")
    output <- system(cmd , intern = TRUE)
    mfe <- as.numeric(gsub("\\)","",gsub("\\(","",strsplit(output[2], split = " ")[[1]][2])))
    gq <- ifelse( grepl("\\+",output[2]) == TRUE, 1, 0)
    dr <- .evalFOLD(output[2])
    return(c(mfe,gq,dr))
}




.evalFOLD <- function(x){
    if(substr(x, 1, 24) == "((((((.(((....))).))))))"){
        return(1)
    } else{
        return(0)
    }
}


.addRNAhybridFeatures <- function(spacers,
                                  bin.RNAhybrid
){
    spacers[["hybMFE_3.12"]] <- .getRNAhybMFE_bulk(dat=spacers,
                                                   POS=3,
                                                   WIDTH=12,
                                                   bin.RNAhybrid=bin.RNAhybrid)
    spacers[["hybMFE_15.9"]] <- .getRNAhybMFE_bulk(dat=spacers,
                                                   POS=15,
                                                   WIDTH=9,
                                                   bin.RNAhybrid=bin.RNAhybrid)
    return(spacers)
}





#' @importFrom Biostrings reverseComplement
#' @importFrom utils write.table
.getRNAhybMFE_bulk <- function(dat,
                               POS=3,
                               WIDTH=12,
                               bin.RNAhybrid
){
  Sys.setenv("PATH" = paste(Sys.getenv('PATH'),
                            bin.RNAhybrid, sep = ':'))  
  
  bashScript <- system.file(package="Cas13design")
  bashScript <- file.path(bashScript, "scripts/RNAhyb.sh")
  
  # transform to DNAstringSet
  GuideSeq = DNAStringSet(dat[["spacer"]])
  # extract guide
  g <- as.character(subseq(x=GuideSeq,
                           start=POS,
                           width=WIDTH))
  # extract target
  t <- as.character(reverseComplement(subseq(x=GuideSeq,
                                             start=POS,
                                             width=WIDTH)))
  # write tmp file to hard disk 
  tmp <- cbind(g, t)
  RanStr <- paste0(sample( letters , size = 6 , replace = F), collapse = "")
  write.table(tmp,
              file=paste0(RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt'),
              sep=",",
              quote=FALSE,
              col.names=FALSE,
              row.names=FALSE)
  # calculate RNA hyb MFE in bulk 
  cmd <- paste0("bash ",
                bashScript, " ",
                paste0(RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt'))
  output <- system(cmd, intern=TRUE)
  # extract the MFE
  hybMFE <- as.numeric(sapply( output , FUN=function(j){as.numeric(strsplit(j, split = ":")[[1]][5])})) # Gets MFE)
  tmp <- file.remove(paste0(RanStr,'_hybMFE_',POS,'.',WIDTH,'.txt'))
  return(hybMFE)
}







#' @importFrom Biostrings writeXStringSet
.addTargetSiteAccessibilityFeatures <- function(dat,
                                                mrnaSequence,
                                                bin.RNAplfold
){
  
  #generate random string for tmp file that will be written
  #to the har drive to avoid colisions
  RanStr = paste0(sample( letters , size = 6 , replace = FALSE), collapse = "")
  names(mrnaSequence) = RanStr
  
  # writing the fa file back to hard drive as a tmp file. This is done to be independent of any naming issues
  writeXStringSet(mrnaSequence,
                  filepath=paste0('./',RanStr,'.fa'),
                  append=FALSE,
                  format="fasta")
  
  cmd <- paste0("cat ",
                paste0('./',RanStr,'.fa')  , " | ",
                bin.RNAplfold , " -L 40 -W 80 -u 50 ")
  output <- system(cmd, intern=TRUE)
  
  UnpairedProbabilities <- .readUnpairedProbabilities(x=paste0('./',RanStr,'_lunp'))
  UnpairedProbabilities.tranformed <- .transform_RNAplfold_predictions(UnpairedProbabilities)
  
 
  # As there was no clear pattern, this will
  # only record the unpaired probability covering the entire guide match
  .getUnpairedProbSimple = function(x, MA){
      dens <- MA[23,][x$pos - 11]
      return(dens)
  }
  Log10_Unpaired <- .getUnpairedProbSimple(x=dat,
                                           MA=log10(UnpairedProbabilities.tranformed))  # in use
  
  # clean up
  tmp = file.remove( paste0('./',RanStr,'.fa') ,paste0('./',RanStr,'_lunp'),  paste0('./',RanStr,'_dp.ps'))
  dat[["Log10_Unpaired"]] <- Log10_Unpaired
  return(dat)
}


#' @importFrom utils read.delim
.readUnpairedProbabilities = function(x){
    out <- read.delim(x,
                      sep="\t",
                      skip=1,
                      row.names="X.i.")[,1:50]
    colnames(out) <- seq(1,50,1)
    out <- t(out)
    return(out)
}

.transform_RNAplfold_predictions <- function(x){
  
  
  ma <- matrix(NA,ncol = ncol(x), nrow = nrow(x) )  
  colnames(ma) <- colnames(x)
  rownames(ma) <- rownames(x)
  for (i in 1:nrow(x)){
    
    if((i %% 2) == 0) {
      
      d = (i/2)-1
      
      ma[i, 1:(ncol(x)-d) ]  <- x[i, (1+d):ncol(x) ]
      
    } 
    else {
      
      d = (i-1)/2
      
      ma[i, 1:(ncol(x)-d) ]  <- x[i, (1+d):ncol(x) ]
      
    }
  }
  return(ma)
}




.FEATURES <-  c('normCS', 'MFE','DR','Gquad', 'Log10_Unpaired',
                'hybMFE_3.12','hybMFE_15.9',
                'NTdens_max_A','NTdens_max_C','NTdens_max_G',
                'NTdens_max_T','NTdens_max_AT','NTdens_max_GC',
                'NTdens_min_A','NTdens_min_C','NTdens_min_G',
                'NTdens_min_T','NTdens_min_AT',
                'pA','pC','pG',
                'pAA','pAC','pAG',
                'pAT','pCA','pCC',
                'pCG','pCT','pGA',
                'pGC','pGG','pGT',
                'pTA','pTC','pTG')





.extractPosition <- function(y,
                             end="end"
){
    coord <- strsplit(strsplit(y, split = ":")[[1]][2] ,split = "_")[[1]][1]
    if (end == "end" ){
        out <- as.integer(strsplit(coord, split = "-")[[1]][2])
    } else{
        out <- as.integer(strsplit(coord, split = "-")[[1]][1])
    }
    return(out)
}


#Gets the GC content of a guide given an XStringset
#' @importFrom Biostrings alphabetFrequency
.getGC <- function(y){
    return(sum(alphabetFrequency(y,
                                 baseOnly=TRUE,
                                 as.prob=TRUE)[c("C", "G")]))
}

#' @importFrom Biostrings longestConsecutive
.getConsecutiveBases <- function(y){
    consecutiveBases <- vector()
    consecutiveBases[1] <- longestConsecutive(as.character(y), "A")
    consecutiveBases[2] <- longestConsecutive(as.character(y), "C")
    consecutiveBases[3] <- longestConsecutive(as.character(y), "G")
    consecutiveBases[4] <- longestConsecutive(as.character(y), "T")
    names(consecutiveBases) <- c("A","C","G","T")
    return(consecutiveBases)
}





#' @importFrom Biostrings letterFrequency
#' @importFrom XVector subseq
#' @importFrom BiocGenerics width
.get_NT_density_Vector = function(mrnaSequence,
                                  NT="G",
                                  WINDOW=30
){
  
  D = WINDOW  
  ma <- rep(NA,ncol = width(mrnaSequence))  
  
  
  if((D %% 2) == 0) {
    
    d = D/2
    for (p in 1:width(mrnaSequence)){
      
      if ( (p-(d-1)) < 1 | (p+d) > width(mrnaSequence) ){
        ma[p] <- NA
      } else{
        ma[p] <-  letterFrequency( subseq(mrnaSequence, start=p-(d-1), end=p+d)  , letters = NT ,as.prob = T) 
      }
    }
  } else {
    
    d = (D-1)/2
    for (p in 1:width(mrnaSequence)){
      if ( (p-d) < 1 | (p+d) > width(mrnaSequence) ){
        ma[p] <- NA
      } else{
        ma[p] <-  letterFrequency( subseq(mrnaSequence, start=p-d, end=p+d),
                                  letters=NT,
                                  as.prob=TRUE) 
      }
    }
  }
  return(ma)
}





.addNTpointdensities <- function(DAT,
                                 Vec.list,
                                 mrnaSequence
){

  file <- system.file(package="crisprScore")
  file <- file.path(file,
                    "casrxrf/LocalNTdensityCorrelations.txt")
  cors <- read.delim(file,
                     sep='\t',
                     header=TRUE,
                     stringsAsFactors=FALSE)
  comb <- cors[grep("combined",cors$Screen),]
  max <- comb[grep("max",comb$COR),]
  min <- comb[grep("min",comb$COR),]
  Vec.list <- .getNTdensitities(mrnaSequence=mrnaSequence,
                                min=min,
                                max=max)
  
  .getVal = function(j,vec, POINT = -11){
      # if j is NA, return vector of NAs
      if (is.na(j) == T){
          out <- NA
      } else {
          if ((j + POINT) < 1){
              out <- NA
          } else {
              out <- vec[j + POINT]
          }
      }
      return(out)
  }

  ma.A = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_A"]] , POINT = max[which(max$NT  == "A"),"P"])
  ma.C = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_C"]] , POINT = max[which(max$NT  == "C"),"P"])
  ma.G = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_G"]] , POINT = max[which(max$NT  == "G"),"P"])
  ma.T = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_T"]] , POINT = max[which(max$NT  == "T"),"P"])
  ma.AT = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_AT"]] , POINT = max[which(max$NT  == "AT"),"P"])
  ma.GC = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["max_GC"]] , POINT = max[which(max$NT  == "GC"),"P"])
  
  mi.A = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_A"]] , POINT = min[which(min$NT  == "A"),"P"])
  mi.C = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_C"]] , POINT = min[which(min$NT  == "C"),"P"])
  mi.G = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_G"]] , POINT = min[which(min$NT  == "G"),"P"])
  mi.T = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_T"]] , POINT = min[which(min$NT  == "T"),"P"])
  mi.AT = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_AT"]] , POINT = min[which(min$NT  == "AT"),"P"])
  mi.GC = sapply( as.list(DAT$pos) , FUN = .getVal ,vec=Vec.list[["min_GC"]] , POINT = min[which(min$NT  == "GC"),"P"])
  L = list(ma.A,ma.C,ma.G,ma.T,ma.AT,ma.GC,
           mi.A,mi.C,mi.G,mi.T,mi.AT,mi.GC)
  names(L) = paste0( "NTdens_" , names(Vec.list))
  dens = do.call(cbind , L)
  out = cbind.data.frame(DAT , dens)
  return(out)
}

.getNTdensitities = function(mrnaSequence, min, max){
  NTs  <- c("A","C","G","T","AT","GC")
  ma.A <- .get_NT_density_Vector(mrnaSequence=mrnaSequence,
                                 NT = NTs[1] , WINDOW = max[which(max$NT  == NTs[1]),"W"])
  ma.C <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[2] , WINDOW = max[which(max$NT  == NTs[2]),"W"])
  ma.G <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[3] , WINDOW = max[which(max$NT  == NTs[3]),"W"])
  ma.T <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[4] , WINDOW = max[which(max$NT  == NTs[4]),"W"])
  ma.AT <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[5] , WINDOW = max[which(max$NT == NTs[5]),"W"])
  ma.GC <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[6] , WINDOW = max[which(max$NT == NTs[6]),"W"])
  
  mi.A <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[1] , WINDOW = min[which(min$NT  == NTs[1]),"W"])
  mi.C <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[2] , WINDOW = min[which(min$NT  == NTs[2]),"W"])
  mi.G <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[3] , WINDOW = min[which(min$NT  == NTs[3]),"W"])
  mi.T <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[4] , WINDOW = min[which(min$NT  == NTs[4]),"W"])
  mi.AT <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[5] , WINDOW = min[which(min$NT == NTs[5]),"W"])
  mi.GC <- .get_NT_density_Vector(mrnaSequence=mrnaSequence, NT = NTs[6] , WINDOW = min[which(min$NT == NTs[6]),"W"])
  
  L = list(ma.A,ma.C,ma.G,ma.T,ma.AT,ma.GC,
           mi.A,mi.C,mi.G,mi.T,mi.AT,mi.GC)
  names(L) = c(paste0('max_',NTs),
               paste0('min_',NTs)) 
  return(L)
}




#' @importFrom Biostrings letterFrequency dinucleotideFrequency
#' @importFrom Biostrings DNAStringSet
.addLetterProbs <-  function(x){
    G <- DNAStringSet(x[["spacer"]])
    
    letterSets <- list("A", "C", "G", "T", "AT", "GC")
    nucProbs <- lapply(letterSets, function(lets){
        letterFrequency(G, letters=lets, as.prob=TRUE)
    })
    diNucProbs <- dinucleotideFrequency(G, as.prob=TRUE)
    out <- cbind.data.frame(nucProbs, diNucProbs)
    rownames(out) <- names(G)
    colnames(out) <- c("pA","pC","pG","pT",
                       "pG|pC", "pA|pT", 
                       "pAA","pAC","pAG","pAT",
                       "pCA","pCC","pCG","pCT",
                       "pGA","pGC","pGG","pGT",
                       "pTA","pTC","pTG","pTT")
    x <- cbind.data.frame(x, out)
    return(x)
}



####### Training utilities ########
.getPretrainedModel <- function(){
    file <- crisprScoreData::RFcombined.rds()
    readRDS(file)
}

#' @importFrom utils read.delim
.getTrainingData <- function(){
    file <- system.file(package="crisprScore")
    file <- file.path(file,
                      "casrxrf/Cas13designGuidePredictorInput.csv.gz")
    trainingData <- read.delim(file,
                               header=TRUE,
                               sep=",",
                               stringsAsFactors=FALSE) 
    return(trainingData)
}


.getTrainingModel <- function(){
    trainingData <- .getTrainingData()
    fittedModelList <- trainModel(x=trainingData,
                                  FIELDS=.FEATURES)
    return(fittedModelList)
}



# The effect size value is normCS
#' @importFrom stats quantile
#' @importFrom randomForest randomForest
#' @importFrom stats complete.cases
trainModel <- function(x,
                       FIELDS=.FEATURES,
                       seed=1234
){
    set.seed(seed)

    # remove incomplete entries
    x <-  x[complete.cases(x), FIELDS, drop=FALSE]

    # scale numeric values for training
    numeric <-  sapply(x, is.numeric) & colnames(x) != "normCS"  & colnames(x) != "Gquad"  & colnames(x) != "DR" & grepl("[A,C,G,T]_",colnames(x)) == FALSE #DO NOT SCALE response value to [0,1] interval

    tmpmean <-  apply(x[,numeric], 2, function(x){
        quantile(x, 0.05)
    })
    tmpsd <-  apply(x[,numeric], 2, function(x){
        quantile(x, 0.95)-quantile(x, 0.05)
    }) 
    x[,numeric] = scale(x[,numeric],
                        center=tmpmean,
                        scale=tmpsd) 
    x[,numeric][x[,numeric] > 1] <- 1
    x[,numeric][x[,numeric] < 0] <-  0


    model = randomForest(x[,-1], x[,1],
                         importance=FALSE,
                         ntree=2000)

    out = list(model=model,
               tmpmean=tmpmean,
               tmpsd=tmpsd)
    return(out)
}



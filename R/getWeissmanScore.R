# Function to call python module to get Weissman sgRNA score
# 
# converts data to correct format before calling python module

# crisprai package

##

library(reticulate)
library(basilisk)

# use_condaenv(condaenv="crisprai2", required=TRUE)
# use_python_version(version="2.7.15", required = TRUE)
# use_python(python="~/miniconda3/envs/crisprai2/bin/python2.7", required=TRUE)

basilisk::configureBasiliskEnv(src = "/gstore/data/omni/crispr/piru/crisprScore/R/basilisk_min.R")

##

# input_file_tss  <- "example_input_for_R/tssTable.txt"
# input_file_grna <- "example_input_for_R/sgrnaInfoTable.txt"

input_file_tss  <- "/gstore/data/omni/crispr/piru/crisprScore/inst/python/crisprai/example_input_for_R/tssTable.txt"
input_file_grna <- "/gstore/data/omni/crispr/piru/crisprScore/inst/python/crisprai/example_input_for_R/sgrnaInfoTable.txt"



tssTable <- read.table(input_file_tss,
                       head=TRUE,
                       stringsAsFactors=FALSE)
grnaTable <- read.table(input_file_grna,
                        head=TRUE,
                        stringsAsFactors=FALSE)
tss_df <- tssTable
sgrnaInfo_df <- sgrnaInfoTable <- grnaTable


#' @export 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getWeissmanScore <- function(tss_df,
                             sgrnaInfo_df,
                             verbose=FALSE,
                             modality=c("CRISPRa", "CRISPRi")
){
  
    modality  <- match.arg(modality)
    inputList <- .prepareInputData(tss_df,
                                   sgrnaInfo_df,
                                   verbose=verbose)

    results <- basiliskRun(env=env_crisprai,
                           shared=FALSE,
                           fork=fork,
                           fun=.pyPredictWeissmanScore,
                           modality=modality,
                           tssTable=inputList[["tssTable"]],
                           p1p2Table=inputList[["p1p2Table"]],
                           sgrnaTable=inputList[["sgrnaTable"]],
                           libraryTable=inputList[["libraryTable"]],
                           verbose=verbose)
    
    # results <- .pyPredictWeissmanScore(modality=modality,
    #                                    tssTable=inputList[["tssTable"]],
    #                                    p1p2Table=inputList[["p1p2Table"]],
    #                                    sgrnaTable=inputList[["sgrnaTable"]],
    #                                    libraryTable=inputList[["libraryTable"]],
    #                                    verbose=verbose)

    
    return(results)
}


#' @importFrom reticulate import_from_path
#' @importFrom reticulate py_suppress_warnings
#' @importFrom reticulate r_to_py
.pyPredictWeissmanScore <- function(modality,
                                    tssTable,
                                    p1p2Table,
                                    sgrnaTable,
                                    libraryTable,
                                    verbose
                                    
){
  
    if (.Platform$OS.type=="windows"){
      stop("Weissman score is not available for Windows at the moment.")
    }
    
    
    # isdf <- testthat::expect_is(tss,class = "data.frame")
    tssTable <- r_to_py(tssTable)
    p1p2Table <- r_to_py(p1p2Table)
    sgrnaTable <- r_to_py(sgrnaTable)
    libraryTable <- r_to_py(libraryTable)
    # newdf <- testthat::expect_is(converted,class = "pandas.core.frame.DataFrame")
    
    dir <- system.file("python",
                       "crisprai",
                       package="crisprScore",
                       mustWork=TRUE)
    
    # dir = "/gstore/data/omni/crispr/piru/crisprScore/inst/python/crisprai"
  
    # path = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/WeissmanHorlbeck/CRISPRai/python/"
    # path = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/crisprScore/inst/python/crisprai/"
    
    pyWeissmanScore <- import_from_path("predictWeissmanScore", path=dir)
    
    scores <- py_suppress_warnings(pyWeissmanScore$predictWeissmanScore(tssTable=tssTable,
                                                                        p1p2Table=p1p2Table,
                                                                        sgrnaTable=sgrnaTable,
                                                                        libraryTable=libraryTable,
                                                                        modality=modality,
                                                                        verbose=verbose))
    
    ## test pass df to pandas
    # scores <- py_suppress_warnings(pyWeissmanScore$testfuncpd(df=tssTable))
    
    return(scores)
}



.prepareInputData <- function(tss_df,
                              sgrnaInfo_df,
                              verbose=FALSE
){
    tssTable <- .getTssTable(tss_df)
    if (verbose){
        message("Done creating TSS table.")
    }
    p1p2Table <- .getP1P2Table(tss_df)
    if (verbose){
        message("Done creating p1p2 table.")
    }
    sgrnaTable <- .getSgrnaTable(tss_df, sgrnaInfo_df)
    if (verbose){
        message("Done creating sgRNA table.")
    }
    libraryTable <- .getLibraryTable(tss_df, sgrnaInfo_df)
    if (verbose){
        message("Done creating library table.")
    }

    inputList <- list(tssTable=tssTable,
                      p1p2Table=p1p2Table,
                      sgrnaTable=sgrnaTable,
                      libraryTable=libraryTable)

    inputList <- .removeInconsistentGenes(inputList,
                                          value="strand")
    if (verbose){
        message("Done removing strand mismatches.")    
    }

    inputList <- .removeInconsistentGenes(inputList,
                                          value="chromosome")
    if (verbose){
        message("Done removing chr mismatches")    
    }
    
    inputList <- .removeMissingGenes(inputList)
    if (verbose){
        message("Done removing missing genes")    
    }
    
    return(inputList)
}



.getTssTable <- function(tssTable){

    .validateColumns <- function(tssTable){
        cols <- c("gene_symbol",
                  "promoter",
                  "position",
                  "strand",
                  "chr")
        if (!all(cols %in% colnames(tssTable))){
            missing <- setdiff(cols, colnames(tssTable))
            stop("Some of the mandatory columns are not found in ",
                 "the input TSS table: ",
                 paste0(missing, collapse=", "))
        }
        tssTable <- tssTable[,cols,drop=FALSE]
        return(tssTable)
    }

    .renameColumns <- function(tssTable){
        cols <- c("gene",
                  "transcripts",
                  "position",
                  "strand",
                  "chromosome") 
        colnames(tssTable) <- cols 
        return(tssTable)
    }

    tssTable <- .validateColumns(tssTable)
    tssTable <- .renameColumns(tssTable)
    tssTable$position <- floor(as.numeric(tssTable$position))
    tssTable[["cage peak ranges"]] <- paste0("[(",
                                             tssTable$position,
                                             ", ", 
                                             tssTable$position + 1,
                                             ")]")
    return(tssTable)
}



.getP1P2Table <- function(tssTable){

    .validateColumns <- function(tssTable){
        cols <- c("gene_symbol",
                  "promoter",
                  "chr",
                  "strand",
                  "position")
        if (!all(cols %in% colnames(tssTable))){
            missing <- setdiff(cols, colnames(tssTable))
            stop("Some of the mandatory columns are not found in ",
                 "the input TSS table: ",
                 paste0(missing, collapse=", "))
        }
        tssTable <- tssTable[,cols,drop=FALSE]
        return(tssTable)
    }

    .renameColumns <- function(tssTable){
        cols <- c("gene",
                  "transcript",
                  "chromosome",
                  "strand",
                  "position") 
        colnames(tssTable) <- cols 
        return(tssTable)
    }

    p1p2Table <- .validateColumns(tssTable)
    p1p2Table <- .renameColumns(p1p2Table)
    p1p2Table[["TSS source"]] <- "CAGE, matched peaks"
    p1p2Table$position <- floor(as.numeric(p1p2Table$position))
    p1p2Table[["primary TSS"]] <- paste0("(",
                                         p1p2Table$position,
                                         ", ", 
                                         p1p2Table$position + 1,
                                         ")")
    p1p2Table[["secondary TSS"]] <- p1p2Table[["primary TSS"]]
    p1p2Table$position <- NULL
    return(p1p2Table)
}



.getSgrnaTable <- function(tssTable,
                           sgrnaInfoTable
){

    .validateSpacerSequence <- function(sgrnaInfoTable){
        if (!"spacer_19mer" %in% colnames(sgrnaInfoTable)){
            stop("spacer_19mer must be in sgrnaInfoTable")
        }
        return(sgrnaInfoTable)
    }

    .validateGrnaColumns <- function(sgrnaInfoTable){
        cols <- c("grna_id",
                  "tss_id",
                  "pam_site",
                  "strand")
        if (!all(cols %in% colnames(sgrnaInfoTable))){
            missing <- setdiff(cols, colnames(sgrnaInfoTable))
            stop("Some of the mandatory columns are not found in ",
                 "the input sgRNA info table: ",
                 paste0(missing, collapse=", "))
        }
        sgrnaInfoTable <- sgrnaInfoTable[,cols,drop=FALSE]
        return(sgrnaInfoTable)
    }

    .renameGrnaColumns <- function(sgrnaInfoTable){
        cols <- c("sgId",
                  "tss_id",
                  "position",
                  "strand")
        colnames(sgrnaInfoTable) <- cols 
        return(sgrnaInfoTable)
    }

    # convert Sonata coordinates (*N*GG)
    # to Weissman coordinates (NG*G*)
    .pamGenentechToWeissman <- function(x){
        return(x + 2)
    }


    sgrnaInfoTable <- .validateSpacerSequence(sgrnaInfoTable)
    spacerLength <- nchar(sgrnaInfoTable$spacer_19mer)[1]
    sgrnaInfoTable <- .validateGrnaColumns(sgrnaInfoTable)
    sgrnaInfoTable <- .renameGrnaColumns(sgrnaInfoTable)
    sgrnaInfoTable[["gene_name"]] <- gsub("_.*", "",
                                          sgrnaInfoTable[["tss_id"]])
    sgrnaInfoTable <- sgrnaInfoTable[!is.na(sgrnaInfoTable$strand),,drop=FALSE]
    sgrnaInfoTable <- sgrnaInfoTable[!is.na(sgrnaInfoTable$position),,drop=FALSE]

    sgrnaInfoTable$Sublibrary <- "customLibrary"
    sgrnaInfoTable$length <- spacerLength
    sgrnaInfoTable$pass_score <- "e39m1"

    pam <- .pamGenentechToWeissman(sgrnaInfoTable[["position"]])
    sgrnaInfoTable[["pam coordinate"]] <- pam

    matching_rows <- match(sgrnaInfoTable$tss_id,
                           tssTable$tss_id)
    txCol <- paste0("['", tssTable[["promoter"]], "']")
    sgrnaInfoTable[["transcript_list"]] <- txCol[matching_rows]

    sgrnaInfoTable <- sgrnaInfoTable[,c("sgId", 
                                      "Sublibrary", 
                                      "gene_name", 
                                      "length", 
                                      "pam coordinate", 
                                      "pass_score", 
                                      "position", 
                                      "strand", 
                                      "transcript_list")]

    return(sgrnaInfoTable)
}



.getLibraryTable <- function(tssTable, sgrnaInfoTable){

    .validateGrnaColumns <- function(sgrnaInfoTable){
        cols <- c("grna_id",
                  "tss_id",
                  "spacer_19mer")
        if (!all(cols %in% colnames(sgrnaInfoTable))){
            missing <- setdiff(cols, colnames(sgrnaInfoTable))
            stop("Some of the mandatory columns are not found in ",
                 "the input sgRNA info table: ",
                 paste0(missing, collapse=", "))
        }
        sgrnaInfoTable <- sgrnaInfoTable[,cols,drop=FALSE]
        return(sgrnaInfoTable)
    }

    .renameGrnaColumns <- function(sgrnaInfoTable){
        cols <- c("sgId",
                  "tss_id",
                  "sequence")
        colnames(sgrnaInfoTable) <- cols 
        return(sgrnaInfoTable)
    }

    libraryTable <- .validateGrnaColumns(sgrnaInfoTable)
    libraryTable <- .renameGrnaColumns(libraryTable)
    libraryTable[["gene"]] <- gsub("_.*", "",
                                   libraryTable[["tss_id"]])
    libraryTable[["sublibrary"]] <- "customLibrary"

    matching_rows <- match(libraryTable$tss_id,
                           tssTable$tss_id)
    libraryTable[["transcripts"]] <- tssTable[matching_rows, "promoter"]

    # drop unnecessary columns and reorder
    libraryTable <- libraryTable[,c("sgId", 
                                  "sublibrary", 
                                  "gene", 
                                  "transcripts", 
                                  "sequence")]
    return(libraryTable)
}



.removeInconsistentGenes <- function(inputList,
                                     value=c("strand", "chromosome")
){
    value <- match.arg(value)
    tssTable <- inputList[["tssTable"]]
    dfs <- split(tssTable[[value]],
                 f=tssTable[["gene"]])
    ns <- vapply(dfs, function(x){
        length(unique(x))
    }, FUN.VALUE=1L)
    mismatch_genes <- names(ns)[ns>1]

    if (length(mismatch_genes)>0){
        # inputList <- lapply(inputList, function(df){
        #     col <- colnames(df)[grep("gene", colnames(df))]
        #     df <- df[!df[[col]] %in% mismatch_genes,,drop=FALSE]
        #     return(df)
        # })   
        inputList <- .removeGenes(inputList, mismatch_genes)
    }
    # inputList <- lapply(inputList, as.data.frame)
    return(inputList)
}



.removeMissingGenes <- function(inputList){
    
    tssTable <- inputList[["tssTable"]]
    sgrnaTable <- inputList[["sgrnaTable"]]
    libraryTable <- inputList[["libraryTable"]]
    
    # check if all genes in sgrnaTable and libraryTable are in TSS
    missing_genes <- append(setdiff(sgrnaTable$gene_name, tssTable$gene),
                            setdiff(libraryTable$gene, tssTable$gene))
    missing_genes <- unique(missing_genes)
    
    if (length(missing_genes)>0){
        # inputList <- lapply(inputList, function(df){
        #     col <- colnames(df)[grep("gene", colnames(df))]
        #     df <- df[!df[[col]] %in% missing_genes,,drop=FALSE]
        #     return(df)
        # })
        inputList <- .removeGenes(inputList, missing_genes)
    }
    
    # inputList <- lapply(inputList, as.data.frame)
    return(inputList)
}

.removeGenes <- function(inputList, genes){
    inputList <- lapply(inputList, function(df){
        col <- colnames(df)[grep("gene", colnames(df))]
        df <- df[!df[[col]] %in% genes,,drop=FALSE]
        return(df)
    })
    inputList <- lapply(inputList, as.data.frame)
    return(inputList)
}







## testing ##


scores <- getWeissmanScore(tss_df,
                            sgrnaInfo_df,
                            verbose=TRUE,
                            modality=c("CRISPRa"))

head(scores)


# head(finalList[["tssTable"]])
# head(finalList[["p1p2Table"]])
# head(finalList[["sgrnaTable"]])
# head(finalList[["libraryTable"]])
# 
# tss <- finalList[["tssTable"]]
# p1p2 <- finalList[["p1p2Table"]]
# sgrnainfo <- finalList[["sgrnaTable"]]
# libraryTable <- finalList[["libraryTable"]]



# write.table(tss, file = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/WeissmanHorlbeck/CRISPRai/input_files/sonata_hg38/min_input/new_trx_sonata_tssTable.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# write.table(p1p2, file = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/WeissmanHorlbeck/CRISPRai/input_files/sonata_hg38/min_input/new_trx_sonata_p1p2Table.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# write.table(sgrnainfo, file = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/WeissmanHorlbeck/CRISPRai/input_files/sonata_hg38/min_input/new_trx_sonata_sgrnaInfoTable.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# write.table(libraryTable, file = "/Volumes/GoogleDrive/My Drive/Projects/CRISPR screen libraries/WeissmanHorlbeck/CRISPRai/input_files/sonata_hg38/min_input/new_trx_sonata_libraryTable.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE)

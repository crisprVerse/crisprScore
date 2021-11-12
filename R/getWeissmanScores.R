#' @title Calculate on-target sgRNA activity scores using crisprai
#' @description Use the Weissman scoring method to calculate on-target sgRNA 
#'     activity scores for Cas9-based CRISPR activation (CRISPRa) and CRISPR 
#'     inactivation (CRISPRi) gene perturbation studies. The Weissman algorithm 
#'     incorporates chromatin features, transcription start site, and sequence 
#'     to predict single-guide RNA (sgRNA) activity scores as these can all 
#'     influence activity in CRISPRa and CRISPRi perturbation studies. This 
#'     method currently only works for sgRNAs designed for use with Cas9 and 
#'     hg38 genome assembly.
#'
#' @param tss_df A \code{data.frame} containing transcription start site (TSS) 
#'     data for promoters. Must have these columns: \code{gene_symbol}, 
#'     \code{promoter}, \code{transcripts}, \code{position}, \code{strand}, and
#'     \code{chr}. See below for more details about \code{tss_df}.    
#' @param sgrna_df \code{A data.frame} containing PAM start site and sequence
#'     for each sgRNA. Must have these columns: \code{grna_id}, \code{tss_id}, 
#'     \code{pam_site}, \code{strand}, and \code{spacer_19mer}. See below for 
#'     more details about \code{sgrna_df}. 
#' @param verbose Should messages be printed to the console? Default value is 
#'     \code{TRUE}.
#' @param modality Which mode of perturbation is being used? Must be a 
#'     \code{string} specifying either \code{CRISPRa} or \code{CRISPRi}.
#' @param fork Set to \code{TRUE} to preserve changes to the R
#'     configuration within the session.
#'     
#' @details \code{tss_df} details:
#'     This must be a \code{data.frame} that contains the following columns:
#'     * gene_symbol: HGNC/HUGO gene identifier.
#'     * promoter: promoter ID.
#'     * transcripts: Ensembl transcript identifier.
#'     * position: start position of transcription start site (TSS).
#'     * strand: strand location. Either _+_ or _-_.
#'     * chr: chromosome location for gene. _e.g. chr19_.
#'     
#' @details \code{sgrna_df} details:
#'     This must be a \code{data.frame} that contains the following columns:
#'     * grna_id: unique sgRNA identifier
#'     * tss_id: name of the transcription start site.
#'     * pam_site: position of the __N__ in the _NGG_ PAM sequence.
#'     * strand: strand location. Either _+_ or _-_.
#'     * spacer_19mer: sgRNA spacer sequence.
#' 
#' @return \strong{getWeissmanScore} returns a \code{data.frame} with 
#'     \code{sequence} and \code{score} columns. The Weissman score takes on a 
#'     value between 0 and 1. A higher score indicates higher sgRNA efficiency.
#' 
#' @references 
#' Horlbeck et al. Compact and highly active next-generation libraries for 
#'     CRISPR-mediated gene repression and activation
#'     eLife 2016;5:e19760.
#'     \url{https://doi.org/10.7554/eLife.19760}.
#' 
#' @author Pirunthan Perampalam
#' 
#' @examples 
#' \dontrun{
#' results <- getWeissmanScores(tss_df = tssExampleCrispra,
#'                              sgrna_df = grnaExampleCrispra,
#'                              modality = "CRISPRa")
#' 
#' results <- getWeissmanScores(tss_df = tssExampleCrispri,
#'                              sgrna_df = grnaExampleCrispri,
#'                              modality = "CRISPRi")
#' }
#' 
#' @md
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getWeissmanScores <- function(tss_df,
                              sgrna_df,
                              verbose=FALSE,
                              modality=c("CRISPRa", "CRISPRi"),
                              fork=FALSE
){

    modality  <- match.arg(modality)
    inputList <- .prepareInputData(tss_df,
                                   sgrna_df,
                                   verbose=verbose)

    dir <- system.file("crisprai",
                       "temp_data",
                       package="crisprScore",
                       mustWork=TRUE)
    dir <- "/Users/fortinj2/crisprScore/inst/crisprai/temp_data"

    pickleFile <- paste0(dir, "/", modality, "_model.pkl")
    
    # TO DO: Change these methods to pull from Experiment Hub 
    # specify fasta, chromatin data files, and pickle files
    dnasef = "/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2_lifted_hg38.bigWig"
    fairef = "/wgEncodeOpenChromFaireK562Sig_lifted_hg38.bigWig"
    mnasef = "/wgEncodeSydhNsomeK562Sig_lifted_hg38.bigWig"
    
    fastaFile <- paste0(dir, "/hg38.fa")
    chromatinFiles=c(dnase=paste0(dir, dnasef),
                     faire=paste0(dir, fairef),
                     mnase=paste0(dir, mnasef))

    results <- basiliskRun(env=env_crisprai,
                           shared=FALSE,
                           fork=fork,
                           fun=.pyPredictWeissmanScore,
                           modality=modality,
                           tssTable=inputList[["tssTable"]],
                           p1p2Table=inputList[["p1p2Table"]],
                           sgrnaTable=inputList[["sgrnaTable"]],
                           libraryTable=inputList[["libraryTable"]],
                           pickleFile=pickleFile,
                           fastaFile=fastaFile,
                           chromatinFiles=chromatinFiles,
                           verbose=verbose)

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
                                    pickleFile,
                                    fastaFile,
                                    chromatinFiles,
                                    verbose

){

    if (.Platform$OS.type=="windows"){
      stop("CRISPRai is not available for Windows at the moment.")
    }

    tssTable <- r_to_py(tssTable)
    p1p2Table <- r_to_py(p1p2Table)
    sgrnaTable <- r_to_py(sgrnaTable)
    libraryTable <- r_to_py(libraryTable)
    chromatinFiles <- r_to_py(chromatinFiles)

    dir <- system.file("python",
                       "crisprai",
                       package="crisprScore",
                       mustWork=TRUE)

    pyWeissmanScore <- reticulate::import_from_path("predictWeissmanScores",
                                                    path=dir)
    
    # TO DO: add chromatin file vectors and pickle file path vectors
    
    scores <- py_suppress_warnings(
      pyWeissmanScore$predictWeissmanScore(tssTable=tssTable,
                                           p1p2Table=p1p2Table,
                                           sgrnaTable=sgrnaTable,
                                           libraryTable=libraryTable,
                                           modality=modality,
                                           pickleFile=pickleFile,
                                           fastaFile=fastaFile,
                                           chromatinFiles=chromatinFiles,
                                           verbose=verbose))
    
    scores <- as.data.frame(scores)
    colnames(scores) <- c("score")
    return(scores)
}



.prepareInputData <- function(tss_df,
                              sgrna_df,
                              verbose=FALSE
){
    tssTable <- .getTssTable(tss_df)
    if (verbose){
        message("Done creating TSS table.")
    }
    p1p2Table <- .getP1P2Table(tss_df)
    if (verbose){
        message("Done creating P1P2 table.")
    }
    sgrnaTable <- .getSgrnaTable(tss_df, sgrna_df)
    if (verbose){
        message("Done creating sgRNA table.")
    }
    libraryTable <- .getLibraryTable(tss_df, sgrna_df)
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
                  "transcripts",
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
    sgrnaInfoTable[["pam_coordinate"]] <- pam

    tssTable$ID <- paste(tssTable$gene_symbol, tssTable$promoter, sep = "_")
    matching_rows <- match(sgrnaInfoTable$tss_id,
                           tssTable$ID)
    txCol <- paste0("['", tssTable[["transcripts"]], "']")
    txCol[which(txCol=="['NA']")] <- "['all']"
    sgrnaInfoTable[["transcript_list"]] <- txCol[matching_rows]
    
    sgrnaInfoTable <- sgrnaInfoTable[,c("sgId",
                                      "Sublibrary",
                                      "gene_name",
                                      "length",
                                      "pam_coordinate",
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

    tssTable$ID <- paste(tssTable$gene_symbol, tssTable$promoter, sep = "_")
    matching_rows <- match(libraryTable$tss_id,
                           tssTable$ID)
    libraryTable[["transcripts"]] <- tssTable[matching_rows, "transcripts"]
    libraryTable[which(is.na(libraryTable$transcripts)), "transcripts"] <- "all"

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
        inputList <- .removeGenes(inputList, mismatch_genes)
    }

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
        inputList <- .removeGenes(inputList, missing_genes)
    }

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


# library(crisprDesignGne)
# fastaFile <- getGenomeFasta()
# chromatinFiles <- getChromatinFiles()
# results <- getCrispraiScores(tss_df=tssExampleCrispra,
#                              sgrna_df=sgrnaExampleCrispra,
#                              chromatinFiles=chromatinFiles,
#                              fastaFile=fastaFile,
#                              modality="CRISPRa")


#' @title Calculate on-target sgRNA activity scores for CRISPRa and CRISPRi
#' @description Use the Weissman lab scoring method (library design v2) to
#'     calculate on-target sgRNA  activity scores for Cas9-based CRISPR
#'     activation (CRISPRa) and CRISPR inactivation (CRISPRi) gene perturbation
#'     studies. The algorithm incorporates chromatin features, transcription
#'     start site, and sequence to predict gRNA activity scores.
#'     Only sgRNAs designed for the human genome (hg38 build) using Cas9 are
#'     supported at the moment, and only spacers of length 19 are
#'     supported at the moment. 
#' 
#' @param tss_df A \code{data.frame} specifying coordinates of transcription
#'     start site (TSS) of the targeted promoter regions.
#'     Must have the following columns: \code{gene_symbol}, 
#'     \code{promoter}, \code{transcripts}, \code{position}, \code{strand}, and
#'     \code{chr}. See details section below for more information.
#' @param sgrna_df A \code{A data.frame} specifying coordinates and spacer 
#'     sequences of the sgRNAs to score.  Must have the following columns:
#'     \code{grna_id}, \code{tss_id}, \code{pam_site}, \code{strand}, and
#'     \code{spacer_19mer}. See details section below for more information.
#' @param verbose Should messages be printed to the console?
#'     TRUE by default.
#' @param modality Which mode of perturbation is being used? Must be a 
#'     \code{string} specifying either \code{CRISPRa} or \code{CRISPRi}.
#' @param fastaFile String specifying fasta file of the hg38 genome.
#' @param chromatinFiles Named character vector of length 3 specifying
#'     BigWig files containing chromatin accessibility data. 
#'     
#' @details \code{tss_df} details:
#'     This must be a \code{data.frame} that contains the following columns:
#'     * tss_id: string specifying name of the TSS.
#'     * gene_symbol: string specifying sHGNC/HUGO gene identifier.
#'     * promoter: string specifying promoter ID (e.g. "P1" or "P2").
#'     * transcripts: Ensembl transcript identifier.
#'     * position: start position of TSS in hg38 coordinates.
#'     * strand: strand of the gene/TSS. Must be either _+_ or _-_.
#'     * chr: string specifying chromosome (e.g. "chr1").
#'     
#' @details \code{sgrna_df} details:
#'     This must be a \code{data.frame} that contains the following columns:
#'     * grna_id: string specifying a unique sgRNA identifier.
#'     * tss_id: string specifying name of the TSS.
#'     * pam_site: genomic ccoordinate of the __N__ in the _NGG_ PAM sequence.
#'     * strand: strand fo the sgRNA. Must be either _+_ or _-_.
#'     * spacer_19mer: string specifying sgRNA 19mer spacer sequence.
#' 
#' @return \strong{getCrispraiScores} returns a \code{data.frame} with 
#'     \code{grna_id} and \code{score} columns. The Weissman score takes on a 
#'     value between 0 and 1. A higher score indicates higher sgRNA efficiency.
#' 
#' @references 
#' Horlbeck et al. Compact and highly active next-generation libraries for 
#'     CRISPR-mediated gene repression and activation
#'     eLife 2016;5:e19760.
#'     \url{https://doi.org/10.7554/eLife.19760}.
#' 
#' @author Pirunthan Perampalam, Jean-Philippe Fortin
#' 
#' @examples 
#' \dontrun{
#' results <- getCrispraiScores(tss_df=tssExampleCrispra,
#'                              sgrna_df=sgrnaExampleCrispra,
#'                              modality="CRISPRa")
#' 
#' results <- getCrispraiScores(tss_df=tssExampleCrispri,
#'                              sgrna_df=sgrnaExampleCrispri,
#'                              modality="CRISPRi")
#' }
#' 
#' @md
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
getCrispraiScores <- function(tss_df,
                              sgrna_df,
                              verbose=FALSE,
                              modality=c("CRISPRa", "CRISPRi"),
                              fastaFile=NULL,
                              chromatinFiles=NULL
){
    .checkChromatinFiles(chromatinFiles)
    .checkFastaFile(fastaFile)
    .checkTssFrame(tss_df)
    .checkGrnaFrame(sgrna_df)

    modality  <- match.arg(modality)
    inputList <- .prepareInputData(tss_df,
                                   sgrna_df,
                                   verbose=verbose)

    if (modality=="CRISPRa"){
        pickleFile <- crisprScoreData::CRISPRa_model.pkl()
    } else if (modality=="CRISPRi"){
        pickleFile <- crisprScoreData::CRISPRi_model.pkl()
    }
    
    results <- .pyPredictWeissmanScore(modality=modality,
                                       tssTable=inputList[["tssTable"]],
                                       p1p2Table=inputList[["p1p2Table"]],
                                       sgrnaTable=inputList[["sgrnaTable"]],
                                       libraryTable=inputList[["libraryTable"]],
                                       pickleFile=pickleFile,
                                       fastaFile=fastaFile,
                                       chromatinFiles=chromatinFiles,
                                       verbose=verbose)
    rownames(results) <- sgrna_df[["grna_id"]]
    return(results)
}




#' @importFrom reticulate import_from_path
#' @importFrom reticulate py_suppress_warnings
#' @importFrom reticulate r_to_py
#' @importFrom basilisk.utils activateEnvironment
#' @importFrom basilisk.utils deactivateEnvironment
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
    env <- basilisk::obtainEnvironmentPath(env_crisprai)
    envls <- basilisk.utils::activateEnvironment(env)
    on.exit(basilisk.utils::deactivateEnvironment(envls))

    if (.Platform$OS.type=="windows"){
        stop("CRISPRai is not available for Windows at the moment.")
    }

    # We are going to save the tables to a temporary directory
    .dumpToFile <- function(table, file){
        write.table(table,
                    sep="\t",
                    file=file,
                    quote=FALSE,
                    col.names=TRUE,
                    row.names=FALSE)
    }
    inputDir <- tempdir()
    .dumpToFile(tssTable, file=file.path(inputDir, "tssTable.tsv"))
    .dumpToFile(p1p2Table, file=file.path(inputDir, "p1p2Table.tsv"))
    .dumpToFile(sgrnaTable, file=file.path(inputDir, "sgrnaTable.tsv"))
    .dumpToFile(libraryTable, file=file.path(inputDir, "libraryTable.tsv"))
    

    # Creating a roster file of all files:
    roster <- data.frame(path=chromatinFiles)
    roster$object <- names(chromatinFiles)
    rownames(roster) <- NULL
    roster <- roster[, c("object", "path")]
    roster <- rbind(roster, c("fasta",fastaFile))
    roster <- rbind(roster, c("pickleFile",pickleFile))
    roster <- rbind(roster, c("tssTable", file.path(inputDir, "tssTable.tsv")))
    roster <- rbind(roster, c("p1p2Table", file.path(inputDir, "p1p2Table.tsv")))
    roster <- rbind(roster, c("sgrnaTable", file.path(inputDir, "sgrnaTable.tsv")))
    roster <- rbind(roster, c("libraryTable", file.path(inputDir, "libraryTable.tsv")))
    rosterFile <- file.path(inputDir, "roster.tsv")
    .dumpToFile(roster, file=rosterFile)
    outputFile <- file.path(inputDir, "scores.txt")

    # Ready to call the python to generate the scores:
    program <- system.file("python",
                           "crisprai",
                           "getWeissmanScores.py",
                           package="crisprScore",
                           mustWork=TRUE)

    pyBinary <- basilisk.utils:::getPythonBinary(env)
    system2(c(pyBinary,
              program,
              rosterFile,
              modality,
              outputFile,
              verbose))

    scores <- read.table(outputFile)
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




.checkFastaFile <- function(fasta){
    if (is.null(fasta)){
        stop("Argument fasta cannot be NULL.")
    }
    if (length(fasta)!=1){
        stop("fasta must be a character vector of length 1.")
    }
     if (!file.exists(fasta)){
        stop("Fasta file cannot be found.")
    }
    invisible(TRUE)
}

.checkChromatinFiles <- function(chromatinFiles){
    if (is.null(chromatinFiles)){
        stop("Argument chromatinFiles cannot be NULL.")
    }
    if (length(chromatinFiles)!=3){
        stop("chromatinFiles must be a character vector of length 3.")
    }
    choices <- c("mnase", "faire", "dnase")
    if (!all(sort(names(chromatinFiles))==sort(names(choices)))){
        stop("chromatinFiles must be a character vector with the",
             " following names: mnase, faire, dnase.")
    }
    if (!sum(file.exists(chromatinFiles))==3){
        stop("Some of the chromatin files cannot be found.")
    }
    invisible(TRUE)
}






.checkTssFrame <- function(tssFrame){
    cols <- c("tss_id",
              "gene_symbol",
              "promoter",
              "transcripts",
              "position",
              "strand",
              "chr")
    if (!all(cols %in% colnames(tssFrame))){
        choices <- setdiff(cols, colnames(tssFrame))
        stop("The following columns are missing in the tssFrame: \n \t",
             paste0(choices, collapse=", "),".")
    }
    
    if (sum(is.na(tssFrame$tss_id))>0){
        stop("tss_id has some missing values.")
    }
    if (sum(is.na(tssFrame$gene_symbol))>0){
        stop("gene_symbol has some missing values.")
    }
    if (sum(is.na(tssFrame$promoter))>0){
        stop("promoter has some missing values.")
    }
    #if (sum(is.na(tssFrame$transcripts))>0){
    #    stop("transcripts has some missing values.")
    #}
    if (sum(is.na(tssFrame$position))>0){
        stop("position has some missing values.")
    }
    if (sum(is.na(tssFrame$strand))>0){
        stop("strand has some missing values.")
    }
    if (sum(is.na(tssFrame$chr))>0){
        stop("chr has some missing values.")
    }

    # Check promoters:
    dfs <- split(tssFrame, f=tssFrame$gene_symbol)
    lens <- vapply(dfs, function(df){
        length(unique(df$strand))
    }, FUN.VALUE=1)
    if (any(lens>1)){
        stop("Some genes have promoters with different strand directions.")
    }
    invisible(TRUE)
}

.checkGrnaFrame <- function(grnaFrame){
    cols <- c("grna_id",
              "tss_id",
              "pam_site",
              "strand", 
              "spacer_19mer")
    if (!all(cols %in% colnames(grnaFrame))){
        choices <- setdiff(cols, colnames(grnaFrame))
        stop("The following columns are missing in the grnaFrame: \n \t",
             paste0(choices, collapse=", "),".")
    }
    if (sum(is.na(grnaFrame$grna_id))>0){
        stop("grna_id has some missing values.")
    }
    if (sum(is.na(grnaFrame$tss_id))>0){
        stop("tss_id has some missing values.")
    }
    if (sum(is.na(grnaFrame$pam_site))>0){
        stop("pam_site has some missing values.")
    }
    if (sum(is.na(grnaFrame$strand))>0){
        stop("strand has some missing values.")
    }
    if (sum(is.na(grnaFrame$spacer_19mer))>0){
        stop("spacer_19mer has some missing values.")
    }
    lens <- unique(nchar(grnaFrame$spacer_19mer))
    if (length(lens)!=1){
        stop("Sequences specified in spacer_19mer must",
             " all be of length 19.")
    } else {
        if (lens!=19){
            stop("Sequences specified in spacer_19mer must",
                 " all be of length 19.")
        }
    }
    invisible(TRUE)
}



.checkSequenceInputs <- function(sequences,
                                 return.output=TRUE
){
    stopifnot(is.character(sequences))
    len <- unique(nchar(sequences))
    if (length(len)!=1){
        stop("Provided sequences for must all have identical length.")
    }
    if (return.output){
        return(sequences)
    } else {
        return(invisible(TRUE))
    }
}

#' @import utils
NULL


getPythonBinary <- function(envpath){
    suffix <- if (.Platform$OS.type == "windows") 
        "Scripts/python.exe"
    else "bin/python"
    file.path(envpath, suffix)
}

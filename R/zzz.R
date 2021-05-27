#' @importFrom basilisk setBasiliskFork
.onLoad <- function(libname, pkgname){
    options(reticulate.useImportHook=FALSE)
}

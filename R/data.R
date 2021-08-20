#' Mismatch position weights for the CFD off-target algorithm (Cas9)
#'
#' CFD off-target scoring weights for each nucleotide position and each
#'     pair of DNA:RNA mismatches, for 20-mer spacer sequences.
#'     Position 20 refers to the spacer position adjacent to
#'     the PAM sequence. A higher weight indicates greater mismatch tolerance 
#'     by the Cas9 nuclease.
#' 
#' @format Numeric vector 
#' @source \url{https://dx.doi.org/10.1038\%2Fnbt.3437}
"cfd.mm.weights.cas9"


#' PAM sequence weights for the CFD off-target algorithm (Cas9)
#'
#' CFD off-target PAM weights for each possible 3-mer Cas9 PAM sequences.
#'     A higher weight indicates greater mismatch tolerance 
#'     by the Cas9 nuclease.
#' 
#' @format Numeric vector
#' @source \url{https://dx.doi.org/10.1038\%2Fnbt.3437}
"cfd.pam.weights.cas9"



#' Mismatch position weights for the MIT off-target algorithm (Cas9)
#'
#' MIT off-target scoring weights for each nucleotide position, for 20-mer
#'     spacer sequences. Position 20 refers to the spacer position adjacent to
#'     the PAM sequence. A higher weight indicates greater mismatch tolerance 
#'     by the Cas9 nuclease.
#' 
#' @format Numeric vector 
#' @source \url{https://dx.doi.org/10.1038\%2Fnbt.2647}
"mit.weights"



#' data.frame detailing available scoring methods
#'
#' data.frame detailing available scoring methods with information needed
#'     to extract nucleotide sequences needed by each scoring algorithm.
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{method}{name of the scoring method}
#'   \item{nuclease}{nuclease compatible with the scoring method}
#'   \item{left}{upstream offset (relative to PAM site) to extract nucleotide
#'       sequence needed for scoring}
#'   \item{right}{downstream offset (relative to PAM site) to extract nucleotide
#'        sequence needed for scoring}
#'   \item{type}{type of the scoring algorithm (on-target or off-target)}
#'   \item{len}{length of the nucleotide sequence needed for scoring}
#' }
#' @usage data(scoringMethodsInfo)
"scoringMethodsInfo"




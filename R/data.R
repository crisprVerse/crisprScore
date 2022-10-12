# #' Mismatch position weights for the CFD off-target algorithm (Cas9)
# #'
# #' CFD off-target scoring weights for each nucleotide position and each
# #'     pair of DNA:RNA mismatches, for 20-mer spacer sequences.
# #'     Position 20 refers to the spacer position adjacent to
# #'     the PAM sequence. A higher weight indicates greater mismatch tolerance 
# #'     by the Cas9 nuclease.
# #' 
# #' @format Numeric vector 
# #' @source \url{https://dx.doi.org/10.1038\%2Fnbt.3437}
# "cfd.mm.weights.cas9"


# #' PAM sequence weights for the CFD off-target algorithm (Cas9)
# #'
# #' CFD off-target PAM weights for each possible 3-mer Cas9 PAM sequences.
# #'     A higher weight indicates greater mismatch tolerance 
# #'     by the Cas9 nuclease.
# #' 
# #' @format Numeric vector
# #' @source \url{https://dx.doi.org/10.1038\%2Fnbt.3437}
# "cfd.pam.weights.cas9"


# #' Mismatch position weights for the CFD off-target algorithm (CasRx)
# #'
# #' CFD off-target scoring weights for each nucleotide position and each
# #'     pair of DNA:RNA mismatches, for 27-mer spacer sequences.
# #'     Position 1 refers to the spacer position adjacent to
# #'     the PAM sequence. A higher weight indicates greater mismatch tolerance 
# #'     by the CasRx nuclease.
# #' 
# #' @format Numeric vector 
# "cfd.mm.weights.casrx"


# #' PAM sequence weights for the CFD off-target algorithm (CasRx)
# #'
# #' CFD off-target PAM weights for each possible 1-mer CasRx PAM sequences.
# #'     A higher weight indicates greater mismatch tolerance 
# #'     by the CasRx nuclease.
# #' 
# #' @format Numeric vector
# "cfd.pam.weights.casrx"




# #' Mismatch position weights for the MIT off-target algorithm (Cas9)
# #'
# #' MIT off-target scoring weights for each nucleotide position, for 20-mer
# #'     spacer sequences. Position 20 refers to the spacer position adjacent to
# #'     the PAM sequence. A higher weight indicates greater mismatch tolerance 
# #'     by the Cas9 nuclease.
# #' 
# #' @format Numeric vector 
# #' @source \url{https://dx.doi.org/10.1038\%2Fnbt.2647}
# "mit.weights"



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
#'   \item{label}{proper case-sensitive method name for labeling}
#'   \item{len}{length of the nucleotide sequence needed for scoring}
#' }
#' @usage data(scoringMethodsInfo)
"scoringMethodsInfo"




#' Example CRISPRa gRNAs data.frame for the getCrispraiScores function
#'
#' Example CRISPRa gRNAs data.frame for the getCrispraiScores function. 
#'     The targeted TSSs are described in the object tssExampleCrispra.
#' 
#' @format A data.frame with 5 columns:
#' \describe{
#'   \item{grna_id}{String specifying gRNA unique identifier.}
#'   \item{tss_id}{String specifying the targeted TSS id.}
#'   \item{pam_site}{Genomic coordinate specifying the first nucleotide of 
#'        PAM sequence.}
#'   \item{strand}{Strand of the gRNA strand. Must be "+" or "-".}
#'   \item{spacer_19mer}{String specifying the nucleotide sequence of the 19mer
#'    spacer sequence.}
#' }
#' @usage data(sgrnaExampleCrispra)
"sgrnaExampleCrispra"


#' Example CRISPRi gRNAs data.frame for the getCrispraiScores function
#'
#' Example CRISPRi gRNAs data.frame for the getCrispraiScores function.
#'     The targeted TSSs are described in the object tssExampleCrispri.
#' 
#' @format A data.frame with 5 columns:
#' \describe{
#'   \item{grna_id}{String specifying gRNA unique identifier.}
#'   \item{tss_id}{String specifying the targeted TSS id.}
#'   \item{pam_site}{Genomic coordinate specifying the first nucleotide of 
#'        PAM sequence.}
#'   \item{strand}{Strand of the gRNA strand. Must be "+" or "-".}
#'   \item{spacer_19mer}{String specifying the nucleotide sequence of the 19mer
#'    spacer sequence.}
#' }
#' @usage data(sgrnaExampleCrispri)
"sgrnaExampleCrispri"



#' Example TSS data.frame for the getCrispraiScores function
#'
#' Example TSS data for gRNAs stored in sgrnaExampleCrispra.
#' 
#' @format A data.frame with 7 columns:
#' \describe{
#'   \item{tss_id}{String specifying the targeted TSS id.}
#'   \item{gene_symbol}{String specifying gene symbol.}
#'   \item{promoter}{String specifying the promoter suffix to add to the
#'   gene symbol columns to obtain the unique TSS id.}
#'   \item{transcripts}{Ensembl IDs of the targeted transcript.}
#'   \item{position}{Integer specifying genomic coordinate of the TSS.}
#'   \item{strand}{Strand of TSS. Must be "+" or "-".}
#'   \item{chr}{String specifying the chromosome name.}
#' }
#' @usage data(tssExampleCrispra)
"tssExampleCrispra"


#' Example TSS data.frame for the getCrispraiScores function
#'
#' Example TSS data for gRNAs stored in sgrnaExampleCrispri.
#' 
#' @format A data.frame with 7 columns:
#' \describe{
#'   \item{tss_id}{String specifying the targeted TSS id.}
#'   \item{gene_symbol}{String specifying gene symbol.}
#'   \item{promoter}{String specifying the promoter suffix to add to the
#'   gene symbol columns to obtain the unique TSS id.}
#'   \item{transcripts}{Ensembl IDs of the targeted transcript.}
#'   \item{position}{Integer specifying genomic coordinate of the TSS.}
#'   \item{strand}{Strand of TSS. Must be "+" or "-".}
#'   \item{chr}{String specifying the chromosome name.}
#' }
#' @usage data(tssExampleCrispri)
"tssExampleCrispri"



library(crisprDesignGne)
library(crisprScore)
chromatinFiles <- getChromatinFiles()
fastaFile <- getGenomeFasta()

results <- getCrispraiScores(tss_df=tssExampleCrispra,
                             sgrna_df=sgrnaExampleCrispra,
                             modality="CRISPRa",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)

results <- getCrispraiScores(tss_df=tssExampleCrispri,
                             sgrna_df=sgrnaExampleCrispri,
                             modality="CRISPRi",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)


tssExampleCrispri[3,5] <- tssExampleCrispri[3,5]-500
results <- getCrispraiScores(tss_df=tssExampleCrispri,
                             sgrna_df=sgrnaExampleCrispri,
                             modality="CRISPRi",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)




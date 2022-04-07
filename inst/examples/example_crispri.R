library(crisprDesignGne)
library(crisprScore)
chromatinFiles <- getChromatinFiles()
fastaFile <- getGenomeFasta()

results <- getWeissmanScores(tss_df=tssExampleCrispra,
                             sgrna_df=sgrnaExampleCrispra,
                             modality="CRISPRa",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)

results <- getWeissmanScores(tss_df=tssExampleCrispri,
                             sgrna_df=sgrnaExampleCrispri,
                             modality="CRISPRi",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)


tssExampleCrispri[3,5] <- tssExampleCrispri[3,5]-500
results <- getWeissmanScores(tss_df=tssExampleCrispri,
                             sgrna_df=sgrnaExampleCrispri,
                             modality="CRISPRi",
                             chromatinFiles=chromatinFiles,
                             fastaFile=fastaFile)




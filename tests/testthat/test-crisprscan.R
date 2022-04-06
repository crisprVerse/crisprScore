context('CRISPRscan score')

sequences  <- c("CCCCCCCCCCCTGCTGATGCTAGATAAGGTTGTGA",
                "GGACCTATCGATGCTGATGCTAGATATGGTTGTGA",
                "GGAAAAAAAAAAACTGATGCTAGATACGGTTGTGA",
                "GGACCTATCGATGCTCGTGCTGGGTACGGTTGTGA",
                "GCCCCCCTCGATGCTGATGCTAGATAGGGCACACA")
scores <- c(0.531, 0.531, 0.45, 0.712, 0.618)


test_that('CRISPRscan score', {
    expect_equal(round(getCRISPRscanScores(sequences)$score,3), scores)
})


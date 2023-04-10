context('Off-target scoring methods')

spacer_cas9 <- "ATCGATGCTGATGCTAGATA"
protospacers_cas9 <- c("ATCGATGCTGATGCTAGATA",
                       "ATCGATGCTGATGCTAGATA",
                       "ATCGATGCTGATGCTAGATA",
                       "ATCGATGCTGATGCTAGATA",
                       "TTCGATGCTGATGCTAGATA",
                       "TTCGATGCTGATGCTAGATG",
                       "TTCGATGCTAATCCTAGATG")
pam_cas9 <- c("AGG","AAG", "AGA", "AGT", "AGG", "AGG", "AGG")
#dput(round(getMITScores(spacer, protospacers)$score,3))
#dput(round(getCFDScores(spacer, protospacers)$score,3))
mit_scores      <- c(1, 0.259, 0.069, 0.016, 1, 0.104, 0.003)
cfd_scores_cas9 <- c(1, 0.259, 0.069, 0.016, 1, 0.765, 0.301)


test_that('MIT scores', {
    expect_equal(round(getMITScores(spacer_cas9,
                                    protospacers_cas9,
                                    pam_cas9)$score,3), mit_scores)
})

test_that('CFD scores for Cas9', {
    expect_equal(round(getCFDScores(spacer_cas9,
                                    protospacers_cas9,
                                    pam_cas9)$score,3), cfd_scores_cas9)
})



spacer_casrx <- "ATCGATGCTGATGCTAGATATGTCGTA"
protospacers_casrx <- c("CTCGATGCTGATGCTAGATATGTCGTA",
                        "CTCGATGCTGATGCTAGATATGTCGTA",
                        "ATCGATGCTGATGCTAGCTATGTCGTA",
                        "ATCGATGCTCATGCTAGAGATGTCGTA",
                        "ATCGATGCTGATGCTAGATATGTCGCA",
                        "ATCGATGCTGATGCTAGATATGTCGCA")
pams_casrx <- c("A", "C", "A", "C", "C","G")
cfd_scores_casrx <- c(0.901, 0.901, 0.021, 0.089, 0.918, 0.918)

test_that('CFD scores for CasRx', {
    expect_equal(round(getCFDScores(spacer_casrx,
                                    protospacers_casrx,
                                    pams_casrx,
                                    nuclease="CasRx")$score,3),
    cfd_scores_casrx)
})






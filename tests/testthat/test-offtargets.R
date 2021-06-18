context('Off-target scoring methods')

spacer <- "ATCGATGCTGATGCTAGATA"
protospacers <- c("ATCGATGCTGATGCTAGATAAGG",
                  "ATCGATGCTGATGCTAGATAAAG",
                  "ATCGATGCTGATGCTAGATAAGA",
                  "ATCGATGCTGATGCTAGATAAGT",
                  "TTCGATGCTGATGCTAGATAAGG",
                  "TTCGATGCTGATGCTAGATGAGG",
                  "TTCGATGCTAATCCTAGATGAGG")
#dput(round(getMITScores(spacer, protospacers)$score,3))
#dput(round(getCFDScores(spacer, protospacers)$score,3))
mit_scores <- c(1, 0.259, 0.069, 0.016, 0.583, 0, 0)
cfd_scores <- c(1, 0.259, 0.069, 0.016, 1, 0.765, 0.301)


test_that('MIT scores', {
    expect_equal(round(getMITScores(spacer, protospacers)$score,3), mit_scores)
})

test_that('CFD scores', {
    expect_equal(round(getCFDScores(spacer, protospacers)$score,3), cfd_scores)
})
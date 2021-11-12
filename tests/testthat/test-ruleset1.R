context('On-target Rule Set 1 score')

ruleset1_sequence  <- "ACCTATCGATGCTGATGCTAGATAAGGTTG"
ruleset1_score     <- c(0.144)


test_that('Rule Set 1 score', {
    expect_equal(round(getRuleSet1Scores(ruleset1_sequence)$score,3), ruleset1_score)
})

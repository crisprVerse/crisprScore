context('On-target scoring methods')

azimuth_sequences  <- "ACCTATCGATGCTGATGCTAGATAAGGTTG"
deepcpf1_sequences <- "ACCTTTTAATCGATGCTGATGCTAGATATTAAGT"
enpamgb_sequences  <- "CATGTTTTTTTGGGAACCAATCGATAATCACATT"
lindel_sequences   <- "ACCTTTTAATCGATGCTGATGCTAGATATTAAGTGGCTTTTAATCGATGCTGATGCTAGATATTA"
deephf_sequences   <- "ATCGATGCTGATGCTAGATAAGG"

azimuth_scores     <- c(0.518)
deepcpf1_scores    <- c(0.689)
enpamgb_scores     <- c(0.713)
lindel_scores      <- c(0.798)
deephf_scores_u6   <- c(0.522)
deephf_scores_t7   <- c(0.309)
deephf_scores_esp  <- c(0.170)
deephf_scores_hf   <- c(0.313)

test_that('Azimuth scores', {
  expect_equal(round(getAzimuthScores(azimuth_sequences)$score,3), azimuth_scores)
})

test_that('DeepCpf1 scores', {
  expect_equal(round(getDeepCpf1Scores(deepcpf1_sequences)$score,3), deepcpf1_scores)
})

test_that('enPAM+GB scores', {
  skip_on_os("windows")
  expect_equal(round(getEnPAMGBScores(enpamgb_sequences)$score,3), enpamgb_scores)
})

test_that('Lindel scores', {
  expect_equal(round(getLindelScores(lindel_sequences)$score,3), lindel_scores)
})

test_that('DeepHF scores', {
  skip_on_os("windows")
  expect_equal(round(getDeepHFScores(deephf_sequences, enzyme="WT")$score,3), deephf_scores_u6)
  expect_equal(round(getDeepHFScores(deephf_sequences, enzyme="ESP")$score,3), deephf_scores_esp)
  expect_equal(round(getDeepHFScores(deephf_sequences, enzyme="HF")$score,3), deephf_scores_hf)
  expect_equal(round(getDeepHFScores(deephf_sequences, enzyme="WT", promoter="T7")$score,3), deephf_scores_t7)
})



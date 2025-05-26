
fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
seq <- loadSequence(seq_path = fa_path)
G4_detected <- G4HunterDetect(seq)

mut <- GRanges(
  seqnames = c("seq1", "seq5"),
  ranges = IRanges(start = c(81, 11), end = c(88, 11)),
  strand = "*",
  ref = c("GGGTAGGG", "A"),
  alt = c("G", "AGGGGGGGGGGGGGGGG")
)

mut_G4 <- G4VarImpact(G4_detected, mut, ref_col = "ref", alt_col = "alt")

test_that("Basic input validation", {

  expect_error(filterVarImpact(GRanges()),
               "Input GRanges object 'mut_G4' is empty.")

  expect_error(filterVarImpact(data.frame()),
               "The input object 'mut_G4' must be a GRanges object.")

  expect_error(filterVarImpact(mut_G4),
               "All thresholds.*are NULL")
})

test_that("Threshold parameter validation", {

  expect_error(filterVarImpact(mut_G4, raw_score_threshold = -1),
               "must be a numeric value between 0 and 4")
  expect_error(filterVarImpact(mut_G4, raw_score_threshold = 5),
               "must be a numeric value between 0 and 4")

  expect_error(filterVarImpact(mut_G4, mut_score_threshold = -1),
               "must be a numeric value between 0 and 4")

  expect_error(filterVarImpact(mut_G4, score_diff_threshold = 1),
               "must be a numeric value between -4 and 0")
  expect_error(filterVarImpact(mut_G4, score_diff_threshold = -5),
               "must be a numeric value between -4 and 0")
})

test_that("Metadata storage", {

  filtered <- filterVarImpact(mut_G4,
                              raw_score_threshold = 2,
                              mut_score_threshold = 1.9,
                              score_diff_threshold = -0.3)

  expect_equal(metadata(filtered)$raw_score_threshold, 2)
  expect_equal(metadata(filtered)$mut_score_threshold, 1.9)
  expect_equal(metadata(filtered)$score_diff_threshold, -0.3)

  filtered_partial <- filterVarImpact(mut_G4, score_diff_threshold = -2)
  expect_null(metadata(filtered_partial)$raw_score_threshold)
  expect_equal(metadata(filtered_partial)$score_diff_threshold, -2)

})


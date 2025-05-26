
hg19 <- BSgenome.Hsapiens.UCSC.hg19
chr21_seq <- DNAStringSet(hg19$chr21)
names(chr21_seq) <- "chr21"
G4 <- G4HunterDetect(chr21_seq)
data(snp_gr)

res_snv <- SNVImpactG4(G4, snp_gr,
                       alt_col = "alt")

test_that("test filterSNVImpact with input type error ", {
  expect_error(filterSNVImpact("test"),
               "The input object 'gr' must be a GRanges object.")
})

test_that("test filterSNVImpact with empty input", {
  empty_gr <- GRanges()
  expect_error(filterSNVImpact(empty_gr),
               "Input GRanges object 'gr' is empty.")
})

test_that("test filterSNVImpact with all thresholds are NULL", {
  expect_error(filterSNVImpact(res_snv),
               paste0("All thresholds \\('raw_score_threshold', ",
                      "'mut_score_threshold', and 'score_diff_threshold'\\) ",
                      "are NULL. Please provide at least one."))
})

test_that("test filterSNVImpact with a negative raw_score_threshold", {
  expect_error(filterSNVImpact(res_snv, raw_score_threshold = -1),
               paste0("'raw_score_threshold' must ",
                      "be a numeric value between 0 and 4."))
})

test_that("test filterSNVImpact with a raw_score_threshold greater than 4", {
  expect_error(filterSNVImpact(res_snv, raw_score_threshold = 4.2),
               paste0("'raw_score_threshold' must ",
                      "be a numeric value between 0 and 4."))
})

test_that("test filterSNVImpact with a negative mut_score_threshold", {
  expect_error(filterSNVImpact(res_snv, mut_score_threshold = -1),
               paste0("'mut_score_threshold' must be a ",
                      "numeric value between 0 and 4."))
})

test_that("test filterSNVImpact with a mut_score_threshold greater than 4", {
  expect_error(filterSNVImpact(res_snv, mut_score_threshold = 4.8),
               paste0("'mut_score_threshold' must be a ",
                      "numeric value between 0 and 4."))
})

test_that("test filterSNVImpact with a positive score_diff_threshold", {
  expect_error(filterSNVImpact(res_snv, score_diff_threshold = 0.2),
               paste0("'score_diff_threshold' must be a ",
                      "numeric value between -4 and 0."))
})

test_that("test filterSNVImpact with a score_diff_threshold less than -4", {
  expect_error(filterSNVImpact(res_snv, score_diff_threshold = -4.2),
               paste0("'score_diff_threshold' must be a ",
                      "numeric value between -4 and 0."))
})

test_that("test filterSNVImpact with a defined raw_score_threshold", {
  filtered_gr <- filterSNVImpact(res_snv, raw_score_threshold = 2.52)
  expect_equal(length(filtered_gr), 1)
  expect_equal(filtered_gr$G4.info.score, 2.59)
})

test_that("test filterSNVImpact with a defined mut_score_threshold", {
  filtered_gr <- filterSNVImpact(res_snv, mut_score_threshold = 1.2)
  expect_equal(length(filtered_gr), 52)
  expect_equal(sum(filtered_gr$mut.score), -4.153)
})

test_that("test filterSNVImpact with a defined score_diff_threshold", {
  filtered_gr <- filterSNVImpact(res_snv, score_diff_threshold = -0.35)
  expect_equal(length(filtered_gr), 36)
})

test_that("test filterSNVImpact with combined thresholds", {
  filtered_gr <- filterSNVImpact(res_snv, raw_score_threshold = 1.5,
                                 mut_score_threshold = 1.1,
                                 score_diff_threshold = -0.3)
  expect_equal(length(filtered_gr), 2)
  expect_equal(filtered_gr$G4.info.score, c(1.54, 1.52))
})

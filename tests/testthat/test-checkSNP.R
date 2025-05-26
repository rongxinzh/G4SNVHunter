
test_that("test checkSNV function", {

  gr1 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = 100, width = 1),
    ref = "A",
    alt = "T"
  )
  expect_true(checkSNV(gr1, mode = "wra", ref_col = "ref", alt_col = "alt"))

  gr2 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = c(100, 200, 300), width = 1),
    ref = c("A", "C", "G"),
    alt = c("T", "T", "A")
  )
  expect_true(checkSNV(gr2, mode = "wra", ref_col = "ref", alt_col = "alt"))

  gr3 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = 100, width = 10)
  )
  expect_false(checkSNV(gr3, mode = "w"))

  gr4 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = 100, width = 1),
    ref = "AG",
    alt = "T"
  )
  expect_false(checkSNV(gr4, mode = "r", ref_col = "ref"))
  expect_error(checkSNV(gr4, mode = "r", ref_col = "reference"))

  gr5 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = 100, width = 1),
    ref = "A",
    alt = "AG"
  )
  expect_false(checkSNV(gr5, mode = "a", alt_col = "alt"))

  gr6 <- GRanges(
    seqnames = Rle("chr1"),
    ranges = IRanges(start = 100, width = 1)
  )
  expect_error(checkSNV(gr6, mode = "r", ref_col = "ref"))
  expect_error(checkSNV(gr6, mode = "a", alt_col = "alt"))

  expect_error(checkSNV(gr6, mode = "xyz", ref_col = "ref"))
  expect_error(checkSNV(gr6, mode = "xyz", ref_col = "reference"))

})


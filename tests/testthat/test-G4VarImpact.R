

test_seq <- paste0("ATCATCATGGGAGGGAGGGAGGGATGGGTATCAGTCGATGCTATAGCATCGACTA",
                  "GGGGCAGGGCATGCAGGGATGGGCATGCTAGCATCAGCTACGATCGGGATCAGCTA",
                  "GCACATCAGATCCCTACCCGCCCTCCCCTGATCGCACTATCTAC")
test_seq <- DNAStringSet(test_seq)
names(test_seq) <- "CHR"
G4 <- G4HunterDetect(test_seq)

variants<- GRanges(
  seqnames = c("CHR", "CHR", "CHR", "CHR", "CHR", "CHR"),
  ranges = IRanges(start = c(5, 33, 12, 11, 124, 389),
                   end = c(10, 33, 15, 11,124, 390)),
  strand = "*",
  rsid = c("ID1", "ID2", "ID3", "ID4", "ID5", "ID6"),
  sid = "ts1",
  ref = c("TCATGG", "A", "AGGG", "G", "C", "AG"),
  alt = c("T", "ACCCGACCATCGCACCCCCC, G", "A,T,GGG", "C", "G", "A"),
  other_info = c("test1", "test2", "test3", "test4", "test5", "test6"),
  other_info2 = c("msg1", "msg2", "msg3", "msg4", "msg5", "msg6")
)

test_that("Input validation works correctly", {

  expect_error(G4VarImpact(NULL, variants),
               "'G4' must be a non-null, non-empty GRanges object.")
  expect_error(G4VarImpact(G4, NULL),
               "'variants' must be a non-null, non-empty GRanges object.")

  expect_error(G4VarImpact(GRanges(), variants),
               "'G4' must be a non-null, non-empty GRanges object.")

  test_g4 <- G4
  test_g4$sequence <- NULL
  expect_error(G4VarImpact(test_g4, variants),
               "'G4' must contain a 'sequence' column")
  expect_error(G4VarImpact(G4, variants, ref_col = NULL),
               "'ref_col' cannot be NULL or empty")
  expect_error(G4VarImpact(G4, variants, ref_col = "ref", alt_col = NULL),
               "'alt_col' cannot be NULL or empty")

  expect_error(G4VarImpact(G4, variants,
                           ref_col = "ref",
                           alt_col = "alt",
                           mode = "x"),
               "Invalid mode. Choose 's' for single-variant mode")

  expect_error(
    G4VarImpact(G4, variants,
                ref_col = "ref",
                alt_col = "alt",
                mode = "m"),
    "When mode is 'm', 'sampleid_col' must be specified"
  )
})

test_that("single-variant mode (s) works correctly", {

  result <- G4VarImpact(G4, variants,
                        ref_col = "ref", alt_col = "alt",
                        mode = "s")

  expect_s4_class(result, "GRanges")
  expect_true(nrow(mcols(result)) > 0)
  expect_true("mutated.max_score" %in% colnames(mcols(result)))
  expect_true("score.diff" %in% colnames(mcols(result)))

  no_overlap_vars <- GRanges(
    seqnames = "chr2",
    ranges = IRanges(start = 1000, width = 1),
    ref = "G",
    alt = "A"
  )
  empty_result <- G4VarImpact(
    G4,
    no_overlap_vars,
    ref_col = "ref",
    alt_col = "alt",
    mode = "s"
  )
  expect_equal(length(empty_result), 0)
})

test_that("multi-variant mode (m) works correctly", {

  result <- G4VarImpact(G4, variants,
                        ref_col = "ref", alt_col = "alt",
                        mode = "m",
                        sampleid_col = "sid")

  expect_s4_class(result, "GRanges")
  expect_true(nrow(mcols(result)) > 0)
})

test_that("Special cases are handled properly", {

  ins_var <- GRanges(
    seqnames = "CHR",
    ranges = IRanges(start = 20, width = 1),
    ref = "-",
    alt = "GGGGGGGGGGGG",
    sampleid = "sample1",
    variantid = "var1"
  )
  expect_s4_class(
    G4VarImpact(
      G4,
      ins_var,
      ref_col = "ref",
      alt_col = "alt",
      mode = "s"
    ),
    "GRanges"
  )

  del_var <- GRanges(
    seqnames = "CHR",
    ranges = IRanges(start = 17, width = 3),
    ref = "GGG",
    alt = "-",
    sampleid = "sample1",
    variantid = "var1"
  )
  expect_s4_class(
    G4VarImpact(
      G4,
      del_var,
      ref_col = "ref",
      alt_col = "alt",
      mode = "m",
      sampleid_col = "sampleid"),
    "GRanges")
})

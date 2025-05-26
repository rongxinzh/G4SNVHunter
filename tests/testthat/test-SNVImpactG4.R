
seq <- DNAStringSet(
  c("TGTACGGGGTACGGGACTGGGCATGGGCATGCACCCTGCCCTACCCCTCCCATTACGTAGCTAGCTACGAT",
    "AGTGATCGACTAGCTAGCTAGCGACTAGCAGCTAGCGCCCTACGCCCTGACCCCTACGCCCCATACTAG",
    "GATCGATCGATCAGCGTAGCTAGCATCGATCGATCAGCATGCAG"))
names(seq) <- c("seq1", "seq2", "seq3")
G4_gr <- G4HunterDetect(seq)

snp_gr <- GRanges(seqnames = c("seq1", "seq2", "seq1"),
                  ranges = IRanges(start = c(20, 39, 7), width = 1),
                  strand = "*")
mcols(snp_gr)$REF <- c("G", "C", "G")
mcols(snp_gr)$ALT <- c("C", "G", "A")
mcols(snp_gr)$snp_id <- c("snp1", "snp2", "snp3")
mcols(snp_gr)$sampleid <- c("sample1", "sample2", "sample1")

test_that("test SNVImpactG4 function", {

  result <- SNVImpactG4(G4_gr, snp_gr, alt_col = "ALT", mode = "s")

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 3)
  expect_true("mut.G4.seq" %in% colnames(mcols(result)))
  expect_true("score.diff" %in% colnames(mcols(result)))

  result <- SNVImpactG4(G4_gr, snp_gr, alt_col = "ALT", mode = "m",
                        sampleid_col = "sampleid", snvid_col = "snp_id")

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 2)
  expect_true("snv.ids" %in% colnames(mcols(result)))
  expect_true("sample.ids" %in% colnames(mcols(result)))

})

test_that("test SNVImpactG4 function with errors", {
  expect_error(SNVImpactG4(NULL, NULL),
               "'G4' must be a non-null, non-empty GRanges object.")
  expect_error(SNVImpactG4(NULL, GRanges()),
               "'G4' must be a non-null, non-empty GRanges object.")
  expect_error(SNVImpactG4(NULL, snp_gr),
               "'G4' must be a non-null, non-empty GRanges object.")
  expect_error(SNVImpactG4(GRanges(), NULL),
               "'G4' must be a non-null, non-empty GRanges object.")
  expect_error(SNVImpactG4(G4_gr, NULL),
               "'snvs' must be a non-null, non-empty GRanges object.")
  expect_error(SNVImpactG4(G4_gr, GRanges()),
               "'snvs' must be a non-null, non-empty GRanges object.")
  expect_error(
    SNVImpactG4(G4_gr, snp_gr),
    "snvs must contain a column named 'alt' in mcols.")
  expect_error(
    SNVImpactG4(G4_gr, snp_gr, mode = "x"),
    paste0("Invalid mode. Choose 's' for variant-centric mode ",
           "or 'm' for sample-centric mode."))
  expect_error(
    SNVImpactG4(G4_gr, snp_gr, alt_col = "ALT", sampleid_col = "test"),
    "When mode is 's', 'sampleid_col' should not be specified.")
  expect_error(
    SNVImpactG4(G4_gr, snp_gr, alt_col = "ALT", snvid_col = "test"),
    "When mode is 's', 'snvid_col' should not be specified.")
  expect_error(SNVImpactG4(G4_gr, snp_gr, alt_col = "hello-world"),
               "snvs must contain a column named 'hello-world' in mcols.")
  expect_error(SNVImpactG4(G4_gr, snp_gr, mode = "m",
                           sampleid_col = "test"),
               "snvs must contain a column named 'alt' in mcols.")
  expect_error(
    SNVImpactG4(G4_gr, snp_gr, mode = "m",
                alt_col = "ALT", sampleid_col = "test"),
    "When mode is 'm', 'sampleid_col' must be specified and present in 'snvs'.")
  expect_error(
    SNVImpactG4(G4_gr, snp_gr, mode = "m",
                alt_col = "ALT", sampleid_col = "sampleid"),
    "When mode is 'm', 'snvid_col' must be specified and present in 'snvs'.")
  expect_error(
    SNVImpactG4(G4_gr,
                snp_gr, mode = "m", alt_col = "ALT",
                sampleid_col = "sampleid", snvid_col = "test"),
    "When mode is 'm', 'snvid_col' must be specified and present in 'snvs'.")

})


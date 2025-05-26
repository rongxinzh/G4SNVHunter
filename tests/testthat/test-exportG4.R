
fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
seq <- loadSequence(seq_path = fa_path)
G4_detected <- G4HunterDetect(seq)

test_that("exportG4 throws error with non-GRanges input", {
  expect_error(exportG4(G4 = "notGRanges", filename = "test.txt"),
               "'G4' must be a GRanges object")
})

test_that("exportG4 throws error with empty GRanges", {
  empty_gr <- GRanges()
  expect_error(exportG4(G4 = empty_gr, filename = "test.txt"),
               "'G4' cannot be NULL")
})

test_that("exportG4 throws error with invalid filename extension", {
  expect_error(exportG4(G4_detected, filename = "res.txt2"),
               "File extension must be one of: txt, csv, xlsx.", fixed = TRUE)
})

test_that("exportG4 throws error with invalid include_metadata", {
  expect_error(exportG4(G4_detected, filename = "res.txt",
                        include_metadata = "hello"),
               "'include_metadata' must be a single logical value.", fixed = TRUE)
})

test_that("exportG4 throws error with invalid revcomp_antisense", {
  expect_error(exportG4(G4_detected, filename = "res.txt",
                        revcomp_antisense = "hello-world"),
               "'revcomp_antisense' must be a single logical value.", fixed = TRUE)
})

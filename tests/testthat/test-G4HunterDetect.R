
test_that("test G4HunterDetect function (input)", {

  sequences <- DNAStringSet(c(
    "AGTGAATGGGATGGGAGGAGGGACGGGGTAGTACAGCATGGGTGCGGGATGGGCGTGGGATCGGAGCATG",
    "TAGGTAGCTACGACACCCTGCCCTACCCCTACCCCTATCTCCTAGGCCCTAGCCCTACCCTAGCCCA"
  ))
  names(sequences) <- c("seq1", "seq2")

  result <- G4HunterDetect(sequences)
  expect_s4_class(result, "GRanges")

  expect_true(all(c("score", "max_score", "sequence") %in%
                    names(mcols(result))))

  expect_equal(length(result), 3)

  result <- G4HunterDetect(sequences, include_sequences = FALSE)

  expect_equal(c("score", "max_score"), names(mcols(result)))
  expect_equal(length(result), 3)

  result <- G4HunterDetect(sequences, strand = "p", include_sequences = FALSE)

  expect_equal(c("score", "max_score"), names(mcols(result)))
  expect_equal(length(result), 2)
  expect_equal(as.character(unique(strand(result))), "+")


  result <- G4HunterDetect(sequences, strand = "p", include_sequences = TRUE)

  expect_equal(c("score", "max_score", "sequence"), names(mcols(result)))
  expect_equal(length(result), 2)
  expect_equal(as.character(unique(strand(result))), "+")


  sequences_short <- DNAStringSet(c(
    "AGTGACTGACAATGGGAATGCGGGAGGAGGGACGGGGGTAGTACAGCTAGTCGAC",
    "TAGGTAGCTACGACAC",
    "TAGCTAGCTACGTACGATCGTACGATCGATCGATCAGCTAGCTACGTCATTCAG",
    "TAGCTAGCTACGTACGATCGTACGA",
    "GGGTGGGTGGGTGGTTGGGTGGGTG",
    "GGGTGGGTGGGTGGTTGGGTGGGTGG",
    "GGGTGGGTGGGTGGTTGGGTGGGTGA",
    "GGGTGGGTGGGTGGTTGGGTGGGT",
    "ctcccctacccctaccctaccccct"
  ))
  names(sequences_short) <- paste0("seq", 1:9)

  result <- G4HunterDetect(sequences_short, window_size = 25)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 5)

  sequences_short <- DNAStringSet(c(
    "ATGCTAGCTAGT",
    "TAGGTAGCTACGACAC",
    "AGCTAGCTA"
  ))
  names(sequences_short) <- paste0("seq", 1:3)

  result <- G4HunterDetect(sequences_short, window_size = 25)

})

test_that("test G4HunterDetect function with errors", {

  expect_error(G4HunterDetect(sequences = "TGTAGCGAGGGGGGGACTAGGGGGGCATGCGGG"),
               "'sequences' must be a DNAStringSet object.")

  sequences <- DNAStringSet(c(
    "TAGTCGAGGGCATGCGATCAGCTACGATCGATCGATCAGTCAGCTAGCATCGATGGGTGGGTGGGTGGGG"
  ))
  names(sequences) <- c("seq1")

  expect_error(G4HunterDetect(sequences, threshold = 0.5),
               "'threshold' is not reasonable. It should be between 1 and 4.")

  expect_error(G4HunterDetect(sequences, window_size = 9),
               "'window_size' must be an integer between 10 and 50.")

  expect_error(G4HunterDetect(sequences, strands = "all"),
               "'strands' should be 'b' or 'p', but got 'all' instead.")
})

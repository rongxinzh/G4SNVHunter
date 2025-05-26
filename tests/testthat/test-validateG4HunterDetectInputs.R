
sequence <- DNAStringSet(c(
  "TAGGCATCGAGGGGCTGAGGGCTGGGCTGAGGGCTGGGGCATGG"
))
names(sequence) <- c("seq1")

test_that("test validateG4HunterDetectInputs function", {

  expect_null(validateG4HunterDetectInputs(sequence, 1.5, 25, TRUE, "b"))
  expect_null(validateG4HunterDetectInputs(sequence, 2, 15, FALSE, "p"))

  expect_error(validateG4HunterDetectInputs("ATGCTAG", 0.99, 25, FALSE, "b"),
               "'sequences' must be a DNAStringSet object.")
  expect_error(validateG4HunterDetectInputs(1.2, 0.99, 25, FALSE, "b"),
               "'sequences' must be a DNAStringSet object.")
  expect_error(validateG4HunterDetectInputs(sequence, 0.99, 25, FALSE, "b"),
               "'threshold' is not reasonable. It should be between 1 and 4.")
  expect_error(validateG4HunterDetectInputs(sequence, 6, 25, TRUE, "b"),
               "'threshold' is not reasonable. It should be between 1 and 4.")
  expect_error(validateG4HunterDetectInputs(sequence, 1.5, 20.2, TRUE, "b"),
               "'window_size' must be an integer between 10 and 50.")
  expect_error(validateG4HunterDetectInputs(sequence, 1.5, 25, "yes", "b"),
               paste0("'include_sequences' must be a logical value ",
               "\\(TRUE or FALSE\\)."))
  expect_error(validateG4HunterDetectInputs(sequence, 1.5, 25, TRUE, "x"),
               "'strands' should be 'b' or 'p', but got 'x' instead.")
  expect_error(validateG4HunterDetectInputs(sequence, 1.5, 125, TRUE, "b"),
               "'window_size' must be an integer between 10 and 50.")
})


test_that("test G4HunterScore function", {

  seq1 <- "GGGATCGGGTACGGGGACTGGGACTGGGCA"
  seq2 <- "ccctgatcccgatcgccctagcgcgagtcgactg"
  seq3 <- "gycuqgwy12ygd2yudy3wqdqwdqw"

  expect_equal(G4HunterScore(seq1), 1.57)
  expect_equal(G4HunterScore(seq2), -0.676)

  expect_error(G4HunterScore(seq3),
               "'seq' must only contain A, T, C, G, U, N characters.")

  expect_error(G4HunterScore(12345),
               "'seq' must be a single string.")

  expect_equal(G4HunterScore("gggtaagggatgggtcggg"),
               G4HunterScore("GGGTAAGGGATGGGTCGGG"))

})

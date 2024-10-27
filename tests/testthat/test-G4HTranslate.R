
test_that("test G4HTranslate function", {

  seq1 <- "GGGTGAATGGATGCCCGCCCCCCATGGTG"
  expected_result1 <-c(3, 3, 3, 0, 1, 0, 0, 0, 2, 2, 0, 0, 1, -3, -3, -3,
                       1, -4, -4, -4, -4, -4, -4, 0, 0, 2, 2, 0, 1)

  expect_equal(G4HTranslate(seq1), expected_result1)

  seq2 <- "TTTGGGGG"
  expected_result2 <-c(0, 0, 0, 4, 4, 4, 4, 4)

  expect_equal(G4HTranslate(seq2), expected_result2)

  seq3 <- "T"
  expected_result3 <-c(0)

  expect_equal(G4HTranslate(seq3), expected_result3)

})

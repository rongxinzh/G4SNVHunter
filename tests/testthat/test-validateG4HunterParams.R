
test_that("test validateG4HunterParams function", {

  expect_true(validateG4HunterParams(1.5, 25, TRUE, "b"))
  expect_true(validateG4HunterParams(2, 15, FALSE, "p"))

  expect_error(validateG4HunterParams(0.99, 25, FALSE, "b"),
               "'threshold' is not reasonable. It should be between 1 and 4.")
  expect_error(validateG4HunterParams(6, 25, TRUE, "b"),
               "'threshold' is not reasonable. It should be between 1 and 4.")
  expect_error(validateG4HunterParams(1.5, 20.2, TRUE, "b"),
               "'window_size' must be an integer between 10 and 50.")
  expect_error(validateG4HunterParams(1.5, 25, "yes", "b"),
               paste0("'include_sequences' must be a logical value ",
               "\\(TRUE or FALSE\\)."))
  expect_error(validateG4HunterParams(1.5, 25, TRUE, "x"),
               "'strands' should be 'b' or 'p', but got 'x' instead.")
  expect_error(validateG4HunterParams(1.5, 125, TRUE, "b"),
               "'window_size' must be an integer between 10 and 50.")
})

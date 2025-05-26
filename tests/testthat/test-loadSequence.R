
test_that("test loadSequence function (load from data.frame)", {

  seq_df <- data.frame(
    chr = c("seq1", "seq2"),
    sequence = c(
      paste0(rep("G", 100), collapse = ""),
      paste0(rep("A", 200), collapse = "")
    )
  )

  result <- loadSequence(genome_seq = seq_df)
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), 2)
  expect_equal(names(result), c("seq1", "seq2"))
  expect_equal(width(result), c(100, 200))
})

test_that("test loadSequence function (load from fasta)", {

  fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")

  result <- loadSequence(seq_path = fa_path)
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), 5)
  expect_equal(names(result), paste0("seq", 1:5))
  expect_equal(width(result), c(625, 83, 87, 50, 44))

})

test_that("test loadSequence function (load from txt)", {

  txt_path <- system.file("extdata", "seq.txt", package = "G4SNVHunter")

  result <- loadSequence(seq_path = txt_path)
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), 5)
  expect_equal(names(result), paste0("seq", 1:5))
  expect_equal(width(result), c(625, 83, 87, 50, 44))

})

test_that("test loadSequence function with errors", {

  seq_df <- data.frame(
    chr = c("seq1", "seq2"),
    sequence = c(
      paste0(rep("G", 100), collapse = ""),
      paste0(rep("A", 200), collapse = "")
    )
  )

  fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")

  expect_error(
    loadSequence(genome_seq = seq_df,
                 seq_path = fa_path),
    paste0("Both 'genome_seq' and 'seq_path' cannot be provided ",
           "simultaneously. Please specify only one.")
  )

  expect_error(
    loadSequence(genome_seq = list(chr = "seq1", sequence = "G")),
    "'genome_seq' must be a data.frame."
  )

  bad_df <- data.frame(chr = "seq1", sequence = "GGGG", extra = "hello-world")
  expect_error(
    loadSequence(genome_seq = bad_df),
    "The data frame must have two columns: 'chr' and 'sequence'."
  )

  empty_df <- data.frame(chr = character(), sequence = character())
  expect_error(
    loadSequence(genome_seq = empty_df),
    "The data frame must have at least one record."
  )

  expect_error(
    loadSequence(seq_path = "./non_existent_file_quwdgy2d23gdyt23d23hd34.txt"),
    "The file specified in 'seq_path' does not exist."
  )

})

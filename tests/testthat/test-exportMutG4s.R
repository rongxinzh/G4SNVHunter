
fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
seq <- loadSequence(seq_path = fa_path)
G4_detected <- G4HunterDetect(seq)

mut <- GRanges(
  seqnames = c("seq1", "seq5"),
  ranges = IRanges(start = c(81, 11), end = c(88, 11)),
  strand = "*",
  ref = c("GGGTAGGG", "A"),
  alt = c("G", "AGGGGGGGGGGGGGGGG")
  )

mut_G4 <- G4VarImpact(G4_detected, mut, ref_col = "ref", alt_col = "alt")
filtered_mut_G4 <- filterVarImpact(mut_G4, score_diff_threshold = -0.3)

test_that("exportMutG4 errors on invalid mut_G4 input", {
  expect_error(exportMutG4(NULL, "test.txt"),
               "'mut_G4' must be a GRanges object.")
  expect_error(exportMutG4(GRanges(), "test.txt"),
               "'mut_G4' cannot be NULL. Please check your input.")
})

test_that("exportMutG4 errors on invalid filename input", {
  expect_error(exportMutG4(filtered_mut_G4, NULL),
               "'filename' should be a single character string")
  expect_error(exportMutG4(filtered_mut_G4, 123),
               "'filename' should be a single character string")
  expect_error(exportMutG4(filtered_mut_G4, c("a.txt", "b.txt")),
               "'filename' should be a single character string")
  expect_error(exportMutG4(filtered_mut_G4, "test.pdf"),
               "File extension must be one of: txt, csv, xlsx")
})

test_that("exportMutG4 errors when required columns are missing", {
  missing_gr <- filtered_mut_G4
  mcols(missing_gr)$G4.info.sequence <- NULL
  expect_error(
    exportMutG4(missing_gr, "test.txt", revcomp_G4_seq = TRUE),
    "The 'mut_G4' object must contain a 'G4.info.sequence' column"
  )

  missing_gr <- filtered_mut_G4
  mcols(missing_gr)$mutated.G4.seq <- NULL
  expect_error(
    exportMutG4(missing_gr, "test.txt", revcomp_mutG4_seq = TRUE),
    "The 'mut_G4' object must contain a 'mutated.G4.seq' column"
  )
})

test_that("exportMutG4 correctly exports TXT files", {
  txt_file <- tempfile(fileext = ".txt")
  expect_s3_class(exportMutG4(filtered_mut_G4, txt_file), "data.frame")
  expect_true(file.exists(txt_file))

  lines <- readLines(txt_file)
  expect_true(grepl("^#threshold:1.5", lines[1]))

  unlink(txt_file)
})

test_that("exportMutG4 correctly exports CSV files", {
  csv_file <- tempfile(fileext = ".csv")
  expect_s3_class(exportMutG4(filtered_mut_G4, csv_file,
                              include_metadata = FALSE), "data.frame")
  expect_true(file.exists(csv_file))

  df <- read.csv(csv_file)
  expect_true("G_rich_sequence" %in% colnames(df))

  unlink(csv_file)
})

test_that("exportMutG4 correctly exports XLSX files", {
  xlsx_file <- tempfile(fileext = ".xlsx")
  expect_s3_class(exportMutG4(filtered_mut_G4, xlsx_file,
                              revcomp_G4_seq = FALSE), "data.frame")
  expect_true(file.exists(xlsx_file))

  wb <- openxlsx::loadWorkbook(xlsx_file)
  expect_true("G4SNVHunter_Output" %in% names(wb))

  unlink(xlsx_file)
})

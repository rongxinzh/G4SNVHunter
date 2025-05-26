
vcf_file <- system.file("extdata",
                        "example_variants_chr16.vcf",
                        package = "G4SNVHunter")

maf_file <- system.file("extdata",
                        "example_variants_chr16.maf",
                        package = "G4SNVHunter")

test_that("loadVariant correctly loads VCF files", {

  result <- loadVariant(vcf_file, "vcf")
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 5000)
  expect_true(all(c("ID", "REF", "ALT") %in% names(mcols(result))))

  result_no_id <- loadVariant(vcf_file, "vcf", keep_vcf_id = FALSE)
  expect_false("ID" %in% names(mcols(result_no_id)))
})

test_that("loadVariant correctly loads MAF files", {

  result <- loadVariant(maf_file, "maf")
  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 13)
  expect_true(all(c("Reference_Allele", "Tumor_Seq_Allele2") %in%
                    names(mcols(result))))
})

test_that("loadVariant validates inputs correctly", {

  expect_error(loadVariant("test.txt", "invalid_type"),
               "Invalid file type: must be 'vcf' or 'maf'")

  expect_error(loadVariant(vcf_file, "maf"),
               "The file looks like a VCF")
  expect_error(loadVariant(maf_file, "vcf"),
               "The file looks like a MAF")

  expect_error(loadVariant("nonexistent.txt", "vcf"),
               "File does not exist")

  expect_error(loadVariant(vcf_file, "vcf", keep_vcf_id = "not_logical"),
               "keep_vcf_id must be a single logical value")
})

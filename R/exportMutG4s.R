#' Export G4s and Associated Variant Information to File
#'
#' This function exports a \code{GRanges} object containing predicted G4
#' regions affected by variants to a TXT, CSV, or XLSX file. The output
#' includes both the original and mutated G-rich sequences, along with
#' detailed information about the associated variants.
#'
#' @param mut_G4 A \code{GRanges} object returned by \code{filterVarImpact}.
#' Both the \code{G4.info.sequence} and the \code{mutated.G4.seq} columns
#' must be included.
#' @param filename A \code{character} string specifying the output file path.
#' The file extension must be one of:
#' \code{.txt}, \code{.csv}, or \code{.xlsx}.
#' @param include_metadata A \code{logical} value.
#' Whether to include global metadata (e.g., G4 prediction parameters) in the
#' output.
#' Default is \code{TRUE}.
#' @param revcomp_G4_seq A \code{logical} value. Whether to
#' reverse-complement G4 sequences on the antisense (negative) strand.
#' Default is \code{TRUE}.
#' @param revcomp_mutG4_seq A \code{logical} value. Whether to
#' reverse-complement mutated G4 sequences on the antisense (negative) strand.
#' Default is \code{TRUE}.
#'
#' @return Invisibly returns a \code{data.frame} object.
#'
#' @export
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @importFrom tools file_ext
#' @importFrom data.table fwrite
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom S4Vectors metadata
#'
#' @examples
#' library(GenomicRanges)
#' fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
#' seq <- loadSequence(seq_path = fa_path)
#' # Predict G4s
#' G4_detected <- G4HunterDetect(seq)
#' # create mut granges object
#' mut <- data.frame(
#'   chr = c("seq1", "seq5"),
#'   pos = c(81, 11),
#'   ref = c("GGGTAGGG", "A"),
#'   alt = c("G", "AGGGGGGGGGGGGGGGG")
#' )
#'
#' mut <- GRanges(
#'   seqnames = mut$chr,
#'   ranges = IRanges(start = mut$pos, end = mut$pos),
#'   strand = "*",
#'   ref = mut$ref,
#'   alt = mut$alt
#' )
#'
#' mut_G4 <- G4VarImpact(G4_detected, mut, ref_col = "ref", alt_col = "alt")
#' filtered_mut_G4 <- filterVarImpact(mut_G4, score_diff_threshold = -0.2)
#' exportMutG4(filtered_mut_G4, "./result_mut.txt")
#' exportMutG4(filtered_mut_G4, "./result_mut.xlsx", include_metadata = FALSE)
#' exportMutG4(filtered_mut_G4, "./result_mut.csv", revcomp_mutG4_seq = FALSE)
#'
#' # remove all exported files
#' unlink("./result_mut.txt")
#' unlink("./result_mut.xlsx")
#' unlink("./result_mut.csv")
exportMutG4 <- function(mut_G4 = NULL,
                        filename = NULL,
                        include_metadata = TRUE,
                        revcomp_G4_seq = TRUE,
                        revcomp_mutG4_seq = TRUE) {

  validateExportMutG4Inputs(mut_G4,
                            filename,
                            include_metadata,
                            revcomp_G4_seq,
                            revcomp_mutG4_seq)

  df <- as.data.frame(mut_G4, stringsAsFactors = FALSE)

  if (revcomp_G4_seq) {
    if (!"G4.info.sequence" %in% names(df)) {
      stop("The 'mut_G4' object must contain a 'G4.info.sequence' column.")
    }

    df$G_rich_sequence <- df$G4.info.sequence

    antisense_idx <- which(df$strand == "-")
    if (length(antisense_idx) > 0) {
      df$G_rich_sequence[antisense_idx] <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(df$G4.info.sequence[antisense_idx])
        )
      )
    } else {
      message("All G4s are on the positive strand; ",
              "reverse complementation for G4 sequences was not applied.")
    }
  }

  if (revcomp_mutG4_seq) {
    if (!"mutated.G4.seq" %in% names(df)) {
      stop("The 'mut_G4' object must contain a 'mutated.G4.seq' column.")
    }

    df$G_rich_mut_sequence <- df$mutated.G4.seq

    antisense_idx <- which(df$strand == "-")
    if (length(antisense_idx) > 0) {
      df$G_rich_mut_sequence[antisense_idx] <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(df$mutated.G4.seq[antisense_idx])
        )
      )
    } else {
      message("All G4s are on the positive strand; ",
              "reverse complementation for mutated G4 ",
              "sequences was not applied.")
    }
  }

  if (include_metadata) {
    meta_cols <- names(metadata(mut_G4))
    if (length(meta_cols) > 0) {
      metadata_str <- paste(meta_cols,
                            metadata(mut_G4)[meta_cols],
                            sep = ":",
                            collapse = ";")
      metadata_str <- paste0("#", metadata_str)
    } else {
      message("[Note] No metadata columns available in 'mut_G4'.")
    }
  }

  ext <- tolower(file_ext(filename))

  if (ext %in% c("txt", "csv")) {

    if (exists("metadata_str")) {
      writeLines(metadata_str, filename)
      append_mode <- TRUE
    } else {
      append_mode <- FALSE
    }

    data.table::fwrite(
      df,
      file = filename,
      sep = ifelse(ext == "csv", ",", "\t"),
      quote = FALSE,
      col.names = TRUE,
      append = append_mode
    )
  } else {

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "G4SNVHunter_Output")

    if (exists("metadata_str")) {
      openxlsx::writeData(wb, sheet = "G4SNVHunter_Output",
                          metadata_str, startRow = 1, colNames = FALSE)
      data_start_row <- length(metadata_str) + 1
    } else {
      data_start_row <- 1
    }

    openxlsx::writeData(
      wb,
      sheet = "G4SNVHunter_Output",
      df,
      startRow = data_start_row,
      colNames = TRUE
    )

    openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
  }

  invisible(df)
}

#' validateExportMutG4Inputs
#'
#' This function validates the input arguments for the \code{exportMutG4} function.
#'
#' @param mut_G4 A \code{GRanges} object returned by \code{filterVarImpact}.
#' Both the \code{G4.info.sequence} and the \code{mutated.G4.seq} columns
#' must be included.
#' @param filename A \code{character} string specifying the output file path.
#' The file extension must be one of:
#' \code{.txt}, \code{.csv}, or \code{.xlsx}.
#' @param include_metadata A \code{logical} value.
#' Whether to include global metadata (e.g., G4 prediction parameters) in the
#' output.
#' Default is \code{TRUE}.
#' @param revcomp_G4_seq A \code{logical} value. Whether to
#' reverse-complement G4 sequences on the antisense (negative) strand.
#' Default is \code{TRUE}.
#' @param revcomp_mutG4_seq A \code{logical} value. Whether to
#' reverse-complement mutated G4 sequences on the antisense (negative) strand.
#' Default is \code{TRUE}.
#'
#' @return NULL
#' @importFrom tools file_ext
#' @noRd
validateExportMutG4Inputs <- function(mut_G4,
                                      filename,
                                      include_metadata,
                                      revcomp_G4_seq,
                                      revcomp_mutG4_seq) {

  if (!inherits(mut_G4, "GRanges")) {
    stop("'mut_G4' must be a GRanges object.")
  }

  if (length(mut_G4) == 0) {
    stop("'mut_G4' cannot be NULL. Please check your input.")
  }

  if (!is.character(filename) || length(filename) != 1) {
    stop("'filename' should be a single character string.")
  }

  valid_exts <- c("txt", "csv", "xlsx")
  ext <- tolower(file_ext(filename))
  if (!ext %in% valid_exts) {
    stop(sprintf("File extension must be one of: %s",
                 paste(valid_exts, collapse = ", ")))
  }

  if (!is.logical(include_metadata) || length(include_metadata) != 1) {
    stop("'include_metadata' must be a single logical value")
  }

  if (!is.logical(revcomp_G4_seq) || length(revcomp_G4_seq) != 1) {
    stop("'revcomp_G4_seq' must be a single logical value")
  }

  if (!is.logical(revcomp_mutG4_seq) || length(revcomp_mutG4_seq) != 1) {
    stop("'revcomp_mutG4_seq' must be a single logical value")
  }

}

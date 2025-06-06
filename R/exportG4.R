#' Export Predicted G4s to a File
#'
#' This function exports a \code{GRanges} object containing predicted G4s
#' generated by the \code{G4HunterDetect} function from the G4SNVHunter
#' package, to a file in TXT, CSV, or XLSX format.
#'
#' @param G4 A \code{GRanges} object returned by \code{G4HunterDetect}.
#' @param filename A \code{character} string specifying the output file path.
#' The file extension must be one of:
#' \code{.txt}, \code{.csv}, or \code{.xlsx}.
#' @param include_metadata A \code{logical} value.
#' Whether to include global metadata (G4 prediction parameters) in the output.
#' Default is \code{TRUE}.
#' @param revcomp_antisense A \code{logical} value. Whether to
#' reverse-complement sequences on the antisense (negative) strand.
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
#' fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
#' seq <- loadSequence(seq_path = fa_path)
#' # Predict G4s
#' G4_detected <- G4HunterDetect(seq)
#'
#' out_xlsx <- file.path(tempdir(), "results.xlsx")
#' out_txt  <- file.path(tempdir(), "results.txt")
#' out_csv  <- file.path(tempdir(), "results.csv")
#'
#' exportG4(G4_detected, out_xlsx)
#' exportG4(G4_detected, out_txt, include_metadata = FALSE)
#' exportG4(G4_detected, out_csv, revcomp_antisense = FALSE)
#'
#' unlink(c(out_xlsx, out_txt, out_csv))
exportG4 <- function(G4 = NULL,
                     filename = NULL,
                     include_metadata = TRUE,
                     revcomp_antisense = TRUE) {

  validateExportG4Inputs(G4,
                         filename,
                         include_metadata,
                         revcomp_antisense)

  df <- as.data.frame(G4, stringsAsFactors = FALSE)

  if (revcomp_antisense) {
    if (!"sequence" %in% names(df)) {
      stop("The 'G4' object must contain a 'sequence' column.")
    }

    df$G_rich_sequence <- df$sequence

    antisense_idx <- which(df$strand == "-")
    if (length(antisense_idx) > 0) {
      df$G_rich_sequence[antisense_idx] <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(df$sequence[antisense_idx])
        )
      )
    } else {
      message("All G4s are on the positive strand; ",
              "reverse complementation was not applied.")
    }
  }

  if (include_metadata) {
    meta_cols <- names(metadata(G4))
    if (length(meta_cols) > 0) {
      metadata_str <- paste(meta_cols,
                            metadata(G4)[meta_cols],
                            sep = ":",
                            collapse = ";")
      metadata_str <- paste0("#", metadata_str)
    } else {
      message("[Note] No metadata columns available in 'G4'.")
    }
  }

  file_ext <- tolower(file_ext(filename))

  if (file_ext %in% c("txt", "csv")) {

    if (exists("metadata_str")) {
      writeLines(metadata_str, filename)
      append_mode <- TRUE
    } else {
      append_mode <- FALSE
    }

    data.table::fwrite(
      df,
      file = filename,
      sep = ifelse(file_ext == "csv", ",", "\t"),
      quote = FALSE,
      col.names = TRUE,
      append = append_mode
    )
  } else {

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "G4SNVHunter-G4s")

    if (exists("metadata_str")) {
      openxlsx::writeData(wb, sheet = "G4SNVHunter-G4s",
                          metadata_str, startRow = 1, colNames = FALSE)
      data_start_row <- length(metadata_str) + 1
    } else {
      data_start_row <- 1
    }

    openxlsx::writeData(
      wb,
      sheet = "G4SNVHunter-G4s",
      df,
      startRow = data_start_row,
      colNames = TRUE
    )

    openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
  }

  invisible(df)
}

#' validateExportG4Inputs
#'
#' This function validates the input arguments for the \code{loadVariant}
#' function.
#'
#' @param G4 A \code{GRanges} object returned by \code{G4HunterDetect}.
#' @param filename A \code{character} string specifying the output file path.
#' The file extension must be one of:
#' \code{.txt}, \code{.csv}, or \code{.xlsx}.
#' @param include_metadata A \code{logical} value.
#' Whether to include global metadata (G4 prediction parameters) in the output.
#' Default is \code{TRUE}.
#' @param revcomp_antisense A \code{logical} value. Whether to
#' reverse-complement sequences on the antisense (negative) strand.
#' Default is \code{TRUE}.
#'
#' @return NULL
#' @importFrom tools file_ext
#' @noRd
validateExportG4Inputs <- function(G4,
                                   filename,
                                   include_metadata,
                                   revcomp_antisense) {

  if (!inherits(G4, "GRanges")) {
    stop("'G4' must be a GRanges object.")
  }

  if (length(G4) == 0) {
    stop("'G4' cannot be NULL. Please check your input.")
  }

  if (!is.character(filename) || length(filename) != 1) {
    stop("'filename' should be a single character string.")
  }

  valid_exts <- c("txt", "csv", "xlsx")
  file_ext <- tolower(file_ext(filename))
  if (!file_ext %in% valid_exts) {
    stop(sprintf("File extension must be one of: %s",
                 paste(valid_exts, collapse = ", ")), ".")
  }

  if (!is.logical(include_metadata) || length(include_metadata) != 1) {
    stop("'include_metadata' must be a single logical value.")
  }

  if (!is.logical(revcomp_antisense) || length(revcomp_antisense) != 1) {
    stop("'revcomp_antisense' must be a single logical value.")
  }

}

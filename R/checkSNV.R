#' Check the validity of SNVs
#'
#' This function is deprecated and will be removed in a future version.
#'
#' This function checks whether the user-provided SNVs are single
#' nucleotide substitutions.
#'
#' @param snv_gr A \code{GRanges} object containing SNV data.
#' @param mode A character string specifying the checks to be performed.
#' \code{w} for checking if all widths (w) are 1,
#' \code{r} for checking if all ref (r) values are A, T, C, or G,
#' \code{a} for checking if all alt (a) values are A, T, C, or G.
#' @param ref_col Column name for the ref bases in \code{snv_gr}.
#' Default is NULL.
#' @param alt_col Column name for the alt bases in \code{snv_gr}.
#' Default is NULL.
#'
#' @return A logical value indicating whether the user-provided SNVs
#' passed all checks.
#'
#' @export
#'
#' @examples
#'
#' if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'   install.packages("BiocManager")
#' }
#'
#' if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'   BiocManager::install("GenomicRanges")
#' }
#'
#' library(GenomicRanges)
#' gr1 <- GRanges("chr1", IRanges(start = 100, width = 1))
#' # check width ('w')
#' checkSNV(gr1, mode = "w")
#'
#' gr2 <- GRanges(
#'   seqnames = Rle("seq1"),
#'   ranges = IRanges(c(100, 200, 300), width = 1),
#'   ref = c("A", "C", "G"),
#'   alt = c("T", "T", "A")
#' )
#'
#' # check width ('w'), ref ('r'), and alt ('a')
#' checkSNV(gr2, mode = "wra", ref_col = "ref", alt_col = "alt")
#' # check width ('w') and alt ('a')
#' checkSNV(gr2, mode = "wa", alt_col = "alt")
#'
#' gr3 <- GRanges("chr1", IRanges(start = 100, width = 10))
#' # widths should be all one
#' checkSNV(gr3, mode = "w")
#'
#' gr4 <- GRanges(
#'   seqnames = Rle("seq1"),
#'   ranges = IRanges(start = 100, width = 1),
#'   ref = "AG",
#'   alt = "T"
#' )
#'
#' # ref should be all one
#' checkSNV(gr4, mode = "wr", ref_col = "ref")
#' @section Deprecated:
#' This function is no longer supported.
checkSNV <- function(snv_gr = NULL,
                     mode = "wra",
                     ref_col = NULL,
                     alt_col = NULL) {

  .Deprecated(msg = "checkSNV is deprecated.")

  if (!inherits(snv_gr, "GRanges")) {
    stop("'snv_gr' must be a GRanges object")
  }

  if (length(snv_gr) == 0) {
    stop("'snv_gr' is empty. Please provide a non-empty GRanges object.")
  }

  valid_modes <- c("w", "r", "a")
  if (any(!unlist(strsplit(mode, "")) %in% valid_modes)) {
    stop("Invalid 'mode' specified. Use a combination of 'w', 'r', 'a'.")
  }

  if (grepl("r", mode)) {
    if (is.null(ref_col)) {
      stop("Mode 'r' is specified but 'ref_col' is NULL. ",
           "Please provide a valid column name for 'ref'.")
    }
    if (!ref_col %in% colnames(mcols(snv_gr))) {
      stop("The 'ref_col' (", ref_col, ") does not exist in 'snv_gr'.")
    }
    if (all(is.na(mcols(snv_gr)[[ref_col]]))) {
      stop("The 'ref_col' (", ref_col, ") in 'snv_gr' contains only NAs.")
    }
  }

  if (grepl("a", mode)) {
    if (is.null(alt_col)) {
      stop("Mode 'a' is specified but 'alt_col' is NULL. ",
           "Please provide a valid column name for 'alt'.")
    }
    if (!alt_col %in% colnames(mcols(snv_gr))) {
      stop("The 'alt_col' (", alt_col, ") does not exist in 'snv_gr'.")
    }
    if (all(is.na(mcols(snv_gr)[[alt_col]]))) {
      stop("The 'alt_col' (", alt_col, ") in 'snv_gr' contains only NAs.")
    }
  }

  is_valid_nucleotide <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(FALSE)
    } else {
      return(all(x %in% c("A", "T", "C", "G")))
    }
  }

  check_results <- logical(3)

  if (grepl("w", mode)) {
    check_results[1] <- all(width(snv_gr) == 1)
    if (!check_results[1]) {
      message("Some SNVs do not have a width of 1.")
    }
  } else {
    check_results[1] <- TRUE
  }

  if (grepl("r", mode)) {
    ref_values <- mcols(snv_gr)[[ref_col]]
    check_results[2] <- is_valid_nucleotide(ref_values)
    if (!check_results[2]) {
      message("Invalid nucleotides found in ref column. Only (A, T, C, G) ",
              "are supported.")
    }
  } else {
    check_results[2] <- TRUE
  }

  if (grepl("a", mode)) {
    alt_values <- mcols(snv_gr)[[alt_col]]
    check_results[3] <- is_valid_nucleotide(alt_values)
    if (!check_results[3]) {
      message("Invalid nucleotides found in alt column. Only (A, T, C, G) ",
              "are supported.")
    }
  } else {
    check_results[3] <- TRUE
  }

  return(all(check_results))
}

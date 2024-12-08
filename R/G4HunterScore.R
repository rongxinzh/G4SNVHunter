#' Calculate the G4Hunter Score for a Given Sequence
#'
#' This function calculates the G4Hunter score for a given nucleotide sequence,
#' which reflects the ability of that sequence to form a G4 structure.
#'
#' @param seq A single character string representing the nucleotide sequence.
#' Must contain only the characters \code{A}, \code{T}, \code{C}, \code{G},
#' \code{U}, and \code{N}.
#' The length of the sequence should not be too short (e.g., less than 10 bp).
#'
#' @return A numeric value representing the G4Hunter score for the provided
#' sequence.
#'
#' @seealso \code{\link{G4HunterDetect}} for detecting the G4 sequences in a
#' given \code{DNAStringSet} object.
#'
#' @export
#'
#' @examples
#'
#' sequence <- "GGGTAAGGGATGGGTCGGG"
#' score <- G4HunterScore(sequence)
#' print(score)
#' # A negative value indicates that the G4 sequence
#' # is located on the reverse strand
#' sequence <- "GGGTAAGGGATGGGTCGGG"
#' score <- G4HunterScore(sequence)
#' print(score)

G4HunterScore <- function(seq = NULL) {
  if (!is.character(seq) || length(seq) != 1) {
    stop("'seq' must be a single string.")
  }

  valid_bases <- c("A", "T", "C", "G", "U", "N")
  seq_upper <- toupper(seq)
  if (any(!strsplit(seq_upper, NULL)[[1]] %in% valid_bases)) {
    stop("'seq' must only contain A, T, C, G, U, N characters.")
  }

  avg_score <- sum(G4HTranslate(seq_upper)) / length(G4HTranslate(seq_upper))

  return(signif(avg_score, 3))
}

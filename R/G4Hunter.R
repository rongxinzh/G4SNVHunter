#' Detect G4 Sequences Using the G4Hunter Algorithm
#'
#' This function detects G4 sequences from a given \code{DNAStringSet} object
#' using the G4Hunter algorithm.
#'
#' @param sequences A \code{DNAStringSet} object containing the input sequences
#' to be analyzed.
#' @param threshold A numeric value specifying the threshold for the G4Hunter
#' score (absolute value). Default is 1.5. It is not recommended to set the
#' threshold below 1.2.
#' @param window_size An integer specifying the window size (bp) used for
#' prediction. Default is 25. Another commonly used window size is 20.
#' However, 25 is generally preferred.
#' @param include_sequences A logical value (\code{TRUE}/\code{FALSE})
#' indicating whether to include the predicted G4 sequences in the output.
#' Default is \code{TRUE}. Setting this parameter to \code{FALSE} can reduce
#' memory usage, which may be beneficial for extremely large genomes. However,
#' we strongly recommend retaining the sequence information in the output, as
#' it is indispensable for subsequent analysis of the impact of variants on
#' G4 formation.
#' @param strands A character string specifying which strand(s) to consider:
#' \code{"b"} for both strands or \code{"p"} for the positive strand only.
#' Default is \code{"b"}.
#'
#' @return A \code{GRanges} object containing the predicted G4 sequences.
#' The \code{GRanges} object includes the following metadata columns:
#' \describe{
#'   \item{\code{score}}{The final G4Hunter score of the predicted G4 sequence
#'   after merging and pruning.}
#'   \item{\code{max_score}}{The maximum G4Hunter score observed across all
#'   sliding windows covering the G4.}
#'   \item{\code{sequence}}{The sequence of the predicted G4 (if
#'   \code{include_sequences = TRUE}).}
#' }
#' Additionally, the following parameters used during detection are stored
#' in the \code{metadata()} of the returned \code{GRanges} object:
#' \describe{
#'   \item{\code{threshold}}{The G4Hunter score threshold used.}
#'   \item{\code{window_size}}{The window size used.}
#'   \item{\code{include_sequences}}{Whether sequences were included.}
#'   \item{\code{strands}}{The strand(s) considered.}
#' }
#' If no G4 sequences are detected, an empty \code{GRanges} object is returned.
#'
#' @seealso \code{\link{loadSequence}} for loading genome sequences into a
#' \code{DNAStringSet} object.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import progress
#' @importFrom GenomeInfoDb seqlevels
#' @export
#'
#' @examples
#'
#' if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'   install.packages("BiocManager")
#' }
#'
#' if (!requireNamespace("Biostrings", quietly = TRUE)) {
#'   BiocManager::install("Biostrings")
#' }
#'
#' library(Biostrings)
#' sequences <- DNAStringSet(c(
#'   "AGTGAATGGGATGGGAGGAGGGACGGGGTAGTACAGCATAGCATG",
#'   "TAGGTAGCTACGACACCCTGCCCTACCCCTACCCCTATCTA"
#' ))
#' names(sequences) <- c("seq1", "seq2")
#'
#' G4s <- G4HunterDetect(sequences, threshold = 1.5, window_size = 25)
#' print(G4s)
#'
#' seq_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
#' G4s <- G4HunterDetect(loadSequence(seq_path = seq_path))
#' print(G4s)
#'
#' seq_path <- system.file("extdata", "seq.txt", package = "G4SNVHunter")
#' G4s <- G4HunterDetect(loadSequence(seq_path = seq_path))
#' print(G4s)
G4HunterDetect <- function(sequences = NULL,
                           threshold = 1.5,
                           window_size = 25,
                           include_sequences = TRUE,
                           strands = "b") {

  validateG4HunterDetectInputs(sequences,
                               threshold,
                               window_size,
                               include_sequences,
                               strands)

  message("Start processing...")

  if (length(sequences) > 2) {
    pb <- progress_bar$new(
      format = "Processing seq No.:current/:total [:bar] :percent in :elapsed",
      total = length(sequences),
      clear = FALSE, width = 60
    )
  } else {
    pb <- NULL
  }

  results <- lapply(seq_along(sequences), function(i) {
    seq <- toupper(as.character(sequences[[i]]))
    chr <- names(sequences)[i]
    if (length(sequences) <= 2) message("Now process: ", chr, ".")
    if (nchar(seq) < window_size) {
      message("The sequence length of ", chr,
              " is less than the window size. ",
              "Skipping this sequence.")
      result <- NULL
    } else {
      result <- G4HunterAnalysis(seq, chr, threshold, window_size,
                                 include_sequences, strands)
    }
    if (!is.null(pb)) {
      pb$tick()
    }

    return(result)
  })

  results <- Filter(Negate(is.null), results)

  if (length(results) == 0) {
    message("No G4 sequences were found.")

    G4_detected <- GRanges()
    metadata(G4_detected)$threshold <- threshold
    metadata(G4_detected)$window_size <- window_size
    metadata(G4_detected)$include_sequences <- include_sequences
    metadata(G4_detected)$strands <- strands

    message("Done!")

    return(G4_detected)
  } else {
    message("Merging predicted G4s...")
    seqlevels_union <- Reduce(union, lapply(results, seqlevels))
    results <- lapply(results, function(gr) {
      seqlevels(gr) <- seqlevels_union
      gr
    })
    G4_detected <- do.call(c, results)

    metadata(G4_detected)$threshold <- threshold
    metadata(G4_detected)$window_size <- window_size
    metadata(G4_detected)$include_sequences <- include_sequences
    metadata(G4_detected)$strands <- strands

    message("Done!")

    return(G4_detected)
  }
}

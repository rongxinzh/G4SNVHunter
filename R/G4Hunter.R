#' Detect G4 Sequences Using G4Hunter Algorithm
#'
#' This function detects G4 sequences for a given \code{DNAStringSet} object
#' using the G4Hunter algorithm.
#'
#' @param sequences A \code{DNAStringSet} object containing the sequences to be
#' analyzed.
#' @param threshold A numeric value specifying the threshold for the G4Hunter
#' score (absolute value). Default is 1.5. We do not recommend setting the
#' threshold below 1.2.
#' @param window_size An integer specifying the window size (bp) for
#' prediction. Default is 25. Another commonly used window size is 20.
#' However, 25 is generally preferred unless otherwise required.
#' @param include_sequences A logical value (\code{TRUE}/\code{FALSE})
#' indicating whether to include the predicted G4 sequences in the output.
#' Default is \code{TRUE}. Setting this parameter to \code{FALSE} can reduce
#' the memory overhead of the final output, which may be useful for extremely
#' large genomes. However, we strongly recommend retaining the sequence
#' information in the output, as it is indispensable for subsequent analysis of
#' the impact of variants on G4 formation.
#' @param strands A character string (\code{"b"} for both strands and
#' \code{"p"} for positive strand only) indicating which strand to be
#' considered. Default is \code{"b"}.
#'
#' @return A \code{GRanges} object containing the predicted G4 sequences.
#' The GRanges object includes the following metadata columns:
#' \describe{
#'   \item{\code{score}}{The actual G4Hunter score for the G4 sequence detected
#'   (after merging, pruning, etc.).}
#'   \item{\code{max_score}}{The maximum G4Hunter score for the window covering
#'   this G4.}
#'   \item{\code{sequence}}{The G4 sequence after merging, pruning, etc.}
#' }
#' If no G4 sequences were detected, it will return an empty \code{GRanges}
#' object.
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

  if (!inherits(sequences, "DNAStringSet")) {
    stop("'sequences' must be a DNAStringSet object.")
  }

  if (validateG4HunterParams(threshold, window_size,
                             include_sequences, strands)) {
    message("Start processing...")
  }

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
    message("Done!")

    return(GRanges())
  } else {
    message("Merging predicted G4s...")
    seqlevels_union <- Reduce(union, lapply(results, seqlevels))
    results <- lapply(results, function(gr) {
      seqlevels(gr) <- seqlevels_union
      gr
    })
    G4_detected <- do.call(c, results)
    message("Done!")

    return(G4_detected)
  }
}

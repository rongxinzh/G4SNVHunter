#' G4HTranslate
#'
#' Translates a DNA sequence into a numeric score vector based on G and C runs.
#'
#' @param sequence A string representing the DNA sequence.
#' @return A numeric vector with scores corresponding to \code{G} and \code{C}
#' runs in the sequence.
#' @import Rcpp
#' @noRd
#' @useDynLib G4SNVHunter
G4HTranslate <- function(sequence) {
  .Call(`_G4SNVHunter_G4HTranslate`, sequence)
}

#' findMaxVal
#'
#' Finds the maximum G4Hunter score (absolute) in a specified range of scores.
#'
#' @param start Integer, the start index of the range.
#' @param end Integer, the end index of the range.
#' @param scores Numeric vector of scores.
#' @return Numeric value representing the maximum absolute score within the
#' specified range.
#' @noRd
findMaxVal <- function(start, end, scores) {
  range_scores <- scores[start:end]
  max_value <- range_scores[which.max(abs(range_scores))]

  return(max_value)
}

#' maxScores
#'
#' Calculates the maximum G4Hunter scores in merged regions.
#'
#' @param gr A \code{GRanges} object containing the G4 ranges and their scores.
#' @param merged_gr A \code{GRanges} object with merged G4s.
#' @return Numeric vector of maximum G4Hunter scores for each merged G4.
#' @import GenomicRanges
#' @import IRanges
#' @noRd
maxScores <- function(gr, merged_gr) {

  overlaps <- findOverlaps(gr, merged_gr)

  gr_scores <- mcols(gr)$max_score

  result <- tapply(
    gr_scores[queryHits(overlaps)], subjectHits(overlaps),
    function(x) x[which.max(abs(x))]
  )

  output <- rep(NA_real_, length(merged_gr))
  output[as.integer(names(result))] <- result

  return(output)
}

#' extendAndTrimRanges
#'
#' Extends and trims genomic ranges based on sequence data and specific base
#' runs (G or C).
#'
#' @param gr A \code{GRanges} object with genomic ranges to be extended and
#' trimmed.
#' @param sequence Character string representing the sequence from which ranges
#' are derived.
#' @param include_sequences Logical value (\code{TRUE}/\code{FALSE}) indicating
#' whether to include the predicted G4 sequences in the output.
#' Default is \code{TRUE}.
#'
#' @return A \code{GRanges} object with extended and trimmed ranges.
#'
#' @import GenomicRanges
#' @import IRanges
#' @noRd
extendAndTrimRanges <- function(gr, sequence, include_sequences = FALSE) {
  seq_len <- nchar(sequence)
  sequence_vec <- strsplit(sequence, "")[[1]]

  pos_indices <- which(as.character(strand(gr)) == "+")
  neg_indices <- which(as.character(strand(gr)) == "-")

  extendLeft <- function(start_pos, base) {
    while (start_pos > 1 && sequence_vec[start_pos] == base &&
           sequence_vec[start_pos - 1] == base) {
      start_pos <- start_pos - 1
    }
    return(start_pos)
  }

  extendRight <- function(end_pos, base) {
    while (end_pos < seq_len && sequence_vec[end_pos] == base &&
           sequence_vec[end_pos + 1] == base) {
      end_pos <- end_pos + 1
    }
    return(end_pos)
  }

  trimLeft <- function(start_pos, base) {
    while (start_pos <= seq_len && sequence_vec[start_pos] != base) {
      start_pos <- start_pos + 1
    }
    return(start_pos)
  }

  trimRight <- function(end_pos, base) {
    while (end_pos > 0 && sequence_vec[end_pos] != base) {
      end_pos <- end_pos - 1
    }
    return(end_pos)
  }

  if (length(pos_indices) > 0) {
    pos_ranges <- ranges(gr)[pos_indices]

    pos_starts <- start(pos_ranges)
    pos_ends <- end(pos_ranges)

    new_pos_starts <- vapply(pos_starts, extendLeft,
                             base = "G", FUN.VALUE = numeric(1))
    new_pos_ends <- vapply(pos_ends, extendRight,
                           base = "G", FUN.VALUE = numeric(1))

    new_pos_starts <- vapply(new_pos_starts, trimLeft,
                             base = "G", FUN.VALUE = numeric(1))
    new_pos_ends <- vapply(new_pos_ends, trimRight,
                           base = "G", FUN.VALUE = numeric(1))

    extended_pos_ranges <- GRanges(
      seqnames = seqnames(gr)[pos_indices],
      ranges = IRanges(start = new_pos_starts, end = new_pos_ends),
      strand = strand(gr)[pos_indices],
      max_score = mcols(gr)$max_score[pos_indices]
    )
  } else {
    extended_pos_ranges <- GRanges()
  }

  if (length(neg_indices) > 0) {
    neg_ranges <- ranges(gr)[neg_indices]

    neg_starts <- start(neg_ranges)
    neg_ends <- end(neg_ranges)

    new_neg_starts <- vapply(neg_starts, extendLeft,
                             base = "C", FUN.VALUE = numeric(1))
    new_neg_ends <- vapply(neg_ends, extendRight,
                           base = "C", FUN.VALUE = numeric(1))

    new_neg_starts <- vapply(new_neg_starts, trimLeft,
                             base = "C", FUN.VALUE = numeric(1))
    new_neg_ends <- vapply(new_neg_ends, trimRight,
                           base = "C", FUN.VALUE = numeric(1))

    extended_neg_ranges <- GRanges(
      seqnames = seqnames(gr)[neg_indices],
      ranges = IRanges(start = new_neg_starts, end = new_neg_ends),
      strand = strand(gr)[neg_indices],
      max_score = mcols(gr)$max_score[neg_indices]
    )
  } else {
    extended_neg_ranges <- GRanges()
  }

  extended_ranges <- c(extended_pos_ranges, extended_neg_ranges)

  valid_ranges <- extended_ranges[start(extended_ranges) <=
                                    end(extended_ranges)]

  extended_sequences <- mapply(function(start, end) {
    substr(sequence, start, end)
  }, start(valid_ranges), end(valid_ranges))

  new_scores <- vapply(extended_sequences, function(seq) {
    # mean(G4HTranslate(seq))
    G4HunterScore(seq)
  }, FUN.VALUE = numeric(1))

  mcols(valid_ranges)$score <- new_scores

  if (include_sequences) {
    mcols(valid_ranges)$sequence <- extended_sequences
  }

  return(valid_ranges)
}

#' G4HunterAnalysis
#'
#' Perform G4Hunter analysis on a given DNA sequence to identify G4 regions.
#'
#' @param sequence A character string representing the DNA sequence.
#' @param chr A character string indicating the sequence identifier.
#' @param threshold A numeric value specifying the threshold for the G4Hunter
#' score (absolute value). Default is 1.5.
#' @param window_size An integer specifying the window size for G4Hunter
#' prediction. Default is 25.
#' @param include_sequences A logical value (\code{TRUE}/\code{FALSE})
#' indicating whether to include the predicted G4 sequences in the output.
#' Default is \code{TRUE}.
#' @param strands A character string (\code{"b"} for both strands and
#' \code{"p"} for positive strand only) indicating which strand to consider.
#' Default is \code{"b"}.
#'
#' @return A \code{GRanges} object with detected G4 regions and their scores.
#'
#' @import GenomicRanges
#' @import RcppRoll
#' @noRd
G4HunterAnalysis <- function(sequence, chr, threshold = 1.5, window_size = 25,
                             include_sequences, strands) {
  score <- G4HTranslate(sequence)
  running_mean_score <- roll_mean(score, window_size,
                                  fill = NA, align = "left")

  results <- list()

  if (strands == "b") {
    regions_above_threshold <- which(abs(running_mean_score) > threshold)
  } else {
    regions_above_threshold <- which(running_mean_score > threshold)
  }

  # Merge
  if (length(regions_above_threshold) > 0) {
    start_positions <- c(
      regions_above_threshold[1],
      regions_above_threshold[which(diff(regions_above_threshold) > 1) + 1]
      )
    end_positions <- c(
      regions_above_threshold[which(diff(regions_above_threshold) > 1)],
      regions_above_threshold[length(regions_above_threshold)]
      )

    max_values <- vapply(seq_along(start_positions), function(i) {
      start <- start_positions[i]
      end <- end_positions[i]
      findMaxVal(start, end, running_mean_score)
    }, numeric(1))

    gr <- GRanges(seqnames = chr,
                  ranges = IRanges(start = start_positions,
                                   end = end_positions + window_size - 1),
                  max_score = max_values,
                  strand = ifelse(max_values > 0, "+", "-"))
    merged_gr <- reduce(gr)

    max_scores_values <- maxScores(gr, merged_gr)

    merged_gr_with_scores <- GRanges(
      seqnames = seqnames(merged_gr),
      ranges = ranges(merged_gr),
      strand = strand(merged_gr),
      max_score = max_scores_values
    )

    G4_detect <- extendAndTrimRanges(merged_gr_with_scores,
                                     sequence, include_sequences)
    if (include_sequences) {
      mcols(G4_detect) <-
        mcols(G4_detect)[, c("score", "max_score", "sequence")]
    } else {
      mcols(G4_detect) <- mcols(G4_detect)[, c("score", "max_score")]
    }

    return(G4_detect)
  } else {
    return(NULL)
  }
}

#' validateG4HunterDetectInputs
#'
#'
#' This function validates the inputs for the G4Hunter analysis.
#'
#' @param sequences A \code{DNAStringSet} object containing the input sequences
#' to be analyzed.
#' @param threshold A numeric value specifying the threshold for the G4Hunter
#' score (absolute value). Default is 1.5.
#' @param window_size An integer specifying the window size for G4Hunter
#' prediction. Default is 25.
#' @param include_sequences A logical value (\code{TRUE}/\code{FALSE})
#' indicating whether to include the predicted G4 sequences in the output.
#' Default is \code{TRUE}.
#' @param strands A character string (\code{"b"} for both strands and
#' \code{"p"} for positive strand only) indicating which strand to consider.
#' Default is \code{"b"}.
#'
#' @return NULL
#'
#' @noRd
validateG4HunterDetectInputs <- function(sequences,
                                         threshold,
                                         window_size,
                                         include_sequences,
                                         strands) {

  if (!inherits(sequences, "DNAStringSet")) {
    stop("'sequences' must be a DNAStringSet object.")
  }

  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'threshold' must be a single numeric value.")
  }

  if (threshold > 4 || threshold < 1) {
    stop("'threshold' is not reasonable. It should be between 1 and 4.")
  }

  if (!is.numeric(window_size) || length(window_size) != 1) {
    stop("'window_size' must be a single numeric value.")
  }

  if (window_size %% 1 != 0 || window_size < 10 || window_size > 50) {
    stop("'window_size' must be an integer between 10 and 50.")
  }

  if (!is.logical(include_sequences)) {
    stop("'include_sequences' must be a logical value (TRUE or FALSE).")
  }

  if (!(strands %in% c("b", "p"))) {
    stop("'strands' should be 'b' or 'p', but got '", strands, "' instead.")
  }

}

#' Load Genome Sequences
#'
#' This function loads genomic sequences from various sources, including a
#' FASTA file, a text file with sequence identifiers and corresponding
#' sequences, or a data frame object.
#'
#' @param genome_seq A data frame containing sequence identifiers and
#' corresponding sequences.
#' @param seq_path A character string specifying the file path to a FASTA file
#' (suffixed with \code{.fa}, \code{.fna} or \code{.fasta}) or a text file with
#' sequence identifiers and corresponding sequences (file headers should not be
#' provided). Ignored if \code{genome_seq} is provided.
#'
#' @return A \code{DNAStringSet} object containing the genome sequences.
#'
#' @importFrom data.table fread
#' @importFrom stats setNames
#' @import Biostrings GenomeInfoDb
#' @export
#'
#' @examples
#'
#' # File path for sequences in fasta format
#' fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
#' seq <- loadSequence(seq_path = fa_path)
#' print(seq)
#'
#' # Another example
#' # Load sequences from data.frame
#' seq_df <- data.frame(
#'   chr = c("seq1", "seq2"),
#'   sequence = c(
#'     paste0(rep("G", 100), collapse = ""),
#'     paste0(rep("A", 100), collapse = "")
#'   )
#' )
#' seq <- loadSequence(genome_seq = seq_df)
#' print(seq)

loadSequence <- function(genome_seq = NULL, seq_path = NULL) {
  if (!is.null(genome_seq) && !is.null(seq_path)) {
    stop("Both 'genome_seq' and 'seq_path' cannot be provided ",
         "simultaneously. Please specify only one.")
  }

  if (!is.null(genome_seq)) {
    if (is.data.frame(genome_seq)) {
      if (ncol(genome_seq) != 2) {
        stop("The data frame must have two columns: 'chr' and 'sequence'.")
      }
      if (nrow(genome_seq) < 1) {
        stop("The data frame must have at least one record.")
      }
      colnames(genome_seq) <- c("chr", "sequence")
      dna_set <- DNAStringSet(setNames(genome_seq$sequence, genome_seq$chr))
      seqlengths(dna_set) <- nchar(genome_seq$sequence)
      return(dna_set)
    } else {
      stop("'genome_seq' must be a data.frame.")
    }
  } else if (!is.null(seq_path)) {
    if (file.exists(seq_path)) {
      if (grepl("\\.fasta$", seq_path) ||
          grepl("\\.fa$", seq_path) ||
          grepl("\\.fna$", seq_path)) {
        dna_set <- readDNAStringSet(seq_path)
        seqlengths(dna_set) <- width(dna_set)
        return(dna_set)
      } else {
        seq_df <- fread(seq_path, header = FALSE)
        if (ncol(seq_df) != 2) {
          stop("The file must have two columns: 'chr' and 'sequence'.")
        }
        colnames(seq_df) <- c("chr", "sequence")
        dna_set <- DNAStringSet(setNames(seq_df$sequence, seq_df$chr))
        seqlengths(dna_set) <- nchar(seq_df$sequence)

        return(dna_set)
      }
    } else {
      stop("The file specified in 'seq_path' does not exist.")
    }
  } else {
    stop("Either 'genome_seq' or 'seq_path' must be provided.")
  }
}

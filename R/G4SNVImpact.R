#' Evaluate the Impact of SNVs on G4 Sequences
#'
#' This function is deprecated and will be removed in a future version.
#'
#' This function evaluates the impact of SNVs on G4 formation.
#'
#' @param G4 A \code{GRanges} object representing the G4 regions. This object
#' must include metadata columns for G4 sequence and G4Hunter scores.
#' It should come from the results returned by the
#' \code{G4HunterDetect} function.
#' @param snvs A \code{GRanges} object representing the SNVs. This object must
#' include metadata columns for SNV alternative alleles and, optionally,
#' SNV IDs.
#' @param alt_col A character string specifying the name of the column in
#' \code{snvs} that contains the alternative alleles for the SNVs.
#' The default is \code{"alt"}.
#' @param mode A character string indicating the mode of operation.
#' Set to \code{"s"} to evaluate the impact of individual SNVs on G4 regions
#' one at a time.
#' Set to \code{"m"} to assess the combined impact of multiple SNVs that
#' overlap the same G4 region within a sample.
#' If using \code{"m"} mode, you must specify the \code{sampleid_col} and
#' \code{snvid_col} parameters.
#' @param sampleid_col A character string specifying the name of the column in
#' \code{snvs} that contains the sample IDs.
#' This parameter is required when \code{mode} is \code{"m"}, and is ignored if
#' \code{mode} is \code{"s"}.
#' @param snvid_col A character string specifying the column name in
#' \code{snvs} that contains SNV IDs.
#' Required when \code{mode} is \code{"m"}.
#' Ignored if \code{mode} is \code{"s"}.
#'
#' @return A \code{GRanges} object depending on the \code{mode} parameter:
#' \describe{
#'   \item{Mode "s":}{
#'     For \code{mode = "s"}, the returned \code{GRanges} object includes the
#'     following metadata columns:
#'     \itemize{
#'       \item \code{seqnames} Identifiers for SNVs.
#'       \item \code{ranges} Position of the SNVs.
#'       \item \code{strand} Strand of the SNVs.
#'       \item \code{SNV.info.*} Metadata columns related to each SNV.
#'       \item \code{G4.info.*} Metadata columns from the original G4 object.
#'       \item \code{mut.G4.seq} The mutated G4 sequence after applying the
#'       SNV change.
#'       \item \code{mut.G4.anno.seq} The mutated G4 sequence, with the
#'       mutated bases annotated using square brackets.
#'       \item \code{mut.score} The G4Hunter score of the mutated sequence.
#'       \item \code{score.diff} The difference between the mutated G4Hunter
#'       score and the original G4Hunter score.
#'       The value is calculated as the absolute value of the mutated G4Hunter
#'       score minus the absolute value of the original G4Hunter score.
#'     }
#'   }
#'   \item{Mode "m":}{
#'     For \code{mode = "m"}, the returned \code{GRanges} object includes the
#'     following metadata columns:
#'     \itemize{
#'       \item \code{seqnames} Identifiers for G4 sequences.
#'       \item \code{ranges} Position of the G4 sequences (start and end).
#'       \item \code{strand} Strand of the G4 sequences.
#'       \item \code{G4.info.*} Metadata columns from the original G4 object.
#'       \item \code{snv.ids} Concatenated SNV IDs for all SNVs affecting the
#'       G4 region.
#'       \item \code{sample.ids} A semicolon-separated list of sample IDs
#'       overlapping each G4 region.
#'       \item \code{mut.G4.seq} The mutated G4 sequence after applying the
#'       combined SNV changes.
#'       \item \code{mut.G4.anno.seq} The mutated G4 sequence, with the
#'       mutated bases annotated using square brackets.
#'       \item \code{mut.score} The G4Hunter score of the mutated sequence.
#'       \item \code{score.diff} The difference between the mutated G4Hunter
#'       score and the original G4Hunter score. This value is calculated as the
#'       absolute value of the mutated G4Hunter score minus the absolute value
#'       of the original G4Hunter score.
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{G4HunterDetect}} for detecting the G4 sequences in a
#'                given \code{DNAStringSet} object.
#'          \code{\link{G4HunterScore}} for calculating the G4Hunter scores
#'                for a given sequence.
#'          \code{\link{filterSNVImpact}} for filtering out SNVs with
#'                significant impact.
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom IRanges start
#' @importFrom Biostrings DNAString DNAStringSet replaceAt BString
#' @export
#'
#' @examples
#'
#' if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'   install.packages("BiocManager")
#' }
#'
#' if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   BiocManager::install("GenomicRanges")
#' }
#'
#' if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'   BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#' }
#'
#' library(GenomicRanges)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Load sequence for chromosome 21 (hg19)
#' hg19 <- BSgenome.Hsapiens.UCSC.hg19
#' chr21_seq <- DNAStringSet(hg19$chr21)
#' # Chromosome name is needed
#' names(chr21_seq) <- "chr21"
#'
#' # Detect G4s in human chromosome 21
#' G4 <- G4HunterDetect(chr21_seq)
#'
#' # 's' mode
#' # Load SNPs
#' data(snp_gr)
#'
#' # Obtain SNPs that overlap with G4 regions and assess their impact on G4.
#' # In variant-centric mode ('s'), evaluating each SNP individually.
#' res_snp <- SNVImpactG4(G4, snp_gr, alt_col = "alt")
#' print(res_snp)
#'
#' # 'm' mode
#' # Load SNVs
#' data(snv_gr)
#'
#' # Obtain SNVs that overlap with G4 regions and assess their impact on G4.
#' # In sample-centric mode ('m'), evaluate the combined impact of SNVs on G4s.
#' # Grouped by the sample IDs specified in 'sampleid_col'.
#' res_snv <- SNVImpactG4(G4, snv_gr,
#'   alt_col = "alt",
#'   mode = "m", sampleid_col = "sampleid", snvid_col = "snv_id"
#' )
#' print(res_snv)
#' @section Deprecated:
#' This function is no longer supported.
#' Use \code{\link{G4VarImpact}} instead.
SNVImpactG4 <- function(G4 = NULL,
                        snvs = NULL,
                        alt_col = "alt",
                        mode = "s",
                        sampleid_col = NULL,
                        snvid_col = NULL) {

  .Deprecated("G4VarImpact")

  if (is.null(G4) || !inherits(G4, "GRanges") || length(G4) == 0) {
    stop("'G4' must be a non-null, non-empty GRanges object.")
  }

  if (!"sequence" %in% colnames(mcols(G4))) {
    stop("'G4' must contain a 'sequence' column. Please use the ",
         "'G4HunterDetect' function with 'include_sequences = TRUE'.")
  }

  if (is.null(snvs) || !inherits(snvs, "GRanges") || length(snvs) == 0) {
    stop("'snvs' must be a non-null, non-empty GRanges object.")
  }

  if (!mode %in% c("s", "m")) {
    stop("Invalid mode. Choose 's' for variant-centric mode or ",
         "'m' for sample-centric mode.")
  }

  if (!alt_col %in% colnames(mcols(snvs))) {
    stop("snvs must contain a column named '", alt_col, "' in mcols.")
  }

  if (mode == "m" &&
      (is.null(sampleid_col) || !sampleid_col %in% colnames(mcols(snvs)))) {
    stop("When mode is 'm', 'sampleid_col' must be specified and present ",
      "in 'snvs'.")
  }

  if (mode == "s" && !is.null(sampleid_col)) {
    stop("When mode is 's', 'sampleid_col' should not be specified.")
  }

  if (mode == "m" &&
      (is.null(snvid_col) || !snvid_col %in% colnames(mcols(snvs)))) {
    stop("When mode is 'm', 'snvid_col' must be specified and present in ",
         "'snvs'.")
  }

  if (mode == "s" && !is.null(snvid_col)) {
    stop("When mode is 's', 'snvid_col' should not be specified.")
  }

  if (!checkSNV(snvs, mode = "wa", alt_col = alt_col)) {
    stop("Your SNV data did not pass the check.")
  }

  message("Start processing...")

  overlaps <- findOverlaps(snvs, G4)

  snv_hits <- queryHits(overlaps)
  g4_hits <- subjectHits(overlaps)

  if (length(snv_hits) == 0) {
    message("No overlaps found between SNVs and G4 regions")
    return(GRanges())
  }

  if (mode == "s") {
    snv_G4 <- GRanges(snvs[snv_hits],
                      SNV.info = data.frame(mcols(snvs[snv_hits])),
                      G4.info = data.frame(G4[g4_hits]))
    snv_pos_in_g4 <- start(snvs[snv_hits]) - start(G4[g4_hits]) + 1

    snv_seq_replaced <- mapply(function(seq, pos, alt) {
      new_seq <- replaceAt(DNAString(seq), IRanges(pos, pos), DNAString(alt))
      new_score <- G4HunterScore(as.character(new_seq))
      annotated_seq <- replaceAt(BString(seq), IRanges(pos, pos),
                                 paste0("[", alt, "]"))
      list(sequence = as.character(new_seq),
           annotated_sequence = as.character(annotated_seq),
           new_score = new_score)
    },
    mcols(G4[g4_hits])$sequence,
    snv_pos_in_g4,
    mcols(snvs[snv_hits])[[alt_col]],
    SIMPLIFY = FALSE)

    snv_seq_replaced_df <- do.call(rbind,
                                   lapply(snv_seq_replaced, function(res) {
      data.frame(sequence = res$sequence,
                 annotated_sequence = res$annotated_sequence,
                 new_score = res$new_score,
                 stringsAsFactors = FALSE)
    }))

    mcols(snv_G4)$mut.G4.seq <- snv_seq_replaced_df$sequence
    mcols(snv_G4)$mut.G4.anno.seq <- snv_seq_replaced_df$annotated_sequence
    mcols(snv_G4)$mut.score <- snv_seq_replaced_df$new_score
    mcols(snv_G4)$score.diff <-
      signif(
        abs(mcols(snv_G4)$mut.score) - abs(mcols(snv_G4)$G4.info.score), 3)
  } else if (mode == "m") {
    combined_data <- data.frame(
      snv_hits = snv_hits,
      g4_hits = g4_hits,
      sample_id = mcols(snvs[snv_hits])[[sampleid_col]],
      stringsAsFactors = FALSE
    )

    sample_g4_groups <- split(combined_data,
                              list(combined_data$sample_id,
                                   combined_data$g4_hits),
                              drop = TRUE)

    result_list <- lapply(sample_g4_groups, function(df) {
      sample_id <- unique(df$sample_id)
      g4_idx <- unique(df$g4_hits)
      base_g4 <- G4[g4_idx]
      base_seq <- mcols(base_g4)$sequence
      snvs_subset <- snvs[df$snv_hits]
      snv_pos_in_g4 <- start(snvs_subset) - start(base_g4) + 1
      snv_alts <- mcols(snvs_subset)[[alt_col]]

      mutated_dna <- replaceAt(DNAString(base_seq),
                               IRanges(snv_pos_in_g4, snv_pos_in_g4),
                               DNAStringSet(snv_alts))
      anno_mutated_dna <- replaceAt(BString(base_seq),
                                    IRanges(snv_pos_in_g4, snv_pos_in_g4),
                                    paste0("[", snv_alts, "]"))
      final_seq <- as.character(mutated_dna)
      final_anno_seq <- as.character(anno_mutated_dna)
      final_score <- G4HunterScore(final_seq)
      score_diff <- signif(abs(final_score) - abs(mcols(base_g4)$score), 3)

      snv_ids <- paste(mcols(snvs_subset)[[snvid_col]], collapse = ";")

      GRanges(base_g4,
        G4.info = mcols(base_g4),
        snv.ids = snv_ids,
        sample.ids = sample_id,
        mut.G4.seq = final_seq,
        mut.G4.anno.seq = final_anno_seq,
        mut.score = final_score,
        score.diff = score_diff
      )
    })

    snv_G4 <- Reduce(function(x, y) c(x, y), result_list)
  }

  message("Processing finished!")

  return(snv_G4)
}

#' Filter SNV Impact GRanges Object Based on User-Defined Thresholds
#'
#' This function is deprecated and will be removed in a future version.
#'
#' This function filters the SNV Impact \code{GRanges} object returned by the
#' \code{SNVImpactG4} function based on user-defined thresholds for the
#' \code{G4.info.score}, \code{mut.score}, and \code{score.diff} parameters.
#' This function filters SNVs that may significantly impair the formation of G4
#' structures using customizable filtering criteria. You are not required to
#' specify all three threshold parameters. However, at least one threshold
#' parameter must be provided.
#'
#' @param gr A \code{GRanges} object returned by the \code{SNVImpactG4}
#' function.
#' @param raw_score_threshold A positive numeric value no greater than 4 used
#' as the threshold for the absolute value of \code{G4.info.score}. G4s with
#' an absolute G4Hunter score exceeding this threshold will be retained.
#' If \code{NULL}, this threshold is not applied.
#' @param mut_score_threshold A positive numeric value no greater than 4 used
#' as the threshold for the absolute value of \code{mut.score}. Mutated G4s
#' with an absolute G4Hunter score below this threshold will be retained.
#' If \code{NULL}, this threshold is not applied.
#' @param score_diff_threshold A negative numeric value no less than -4 used as
#' the threshold for \code{score.diff}. G4s with a decrease in G4Hunter score
#' greater than this threshold after variation will be retained. If
#' \code{NULL}, this threshold is not applied.
#'
#' @seealso \code{\link{SNVImpactG4}} for assessing the impact of SNVs on
#' G4 formation.
#'
#' @return A filtered \code{GRanges} object, containing only the records that
#' meet the specified threshold criteria.
#' @import GenomicRanges
#' @export
#'
#' @examples
#'
#' if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'   install.packages("BiocManager")
#' }
#'
#' if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   BiocManager::install("GenomicRanges")
#' }
#'
#' if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'   BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#' }
#'
#' library(GenomicRanges)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' hg19 <- BSgenome.Hsapiens.UCSC.hg19
#' chr21_seq <- DNAStringSet(hg19$chr21)
#' # Chromosome name is needed
#' names(chr21_seq) <- "chr21"
#'
#' G4 <- G4HunterDetect(chr21_seq)
#'
#' data(snp_gr)
#'
#' res_snp <- SNVImpactG4(G4, snp_gr, alt_col = "alt")
#' filtered_snv_eff <- filterSNVImpact(res_snp,
#'   mut_score_threshold = 1.2,
#'   score_diff_threshold = -0.2
#' )
#' print(filtered_snv_eff)
#' @section Deprecated:
#' This function is no longer supported.
#' Use \code{\link{filterVarImpact}} instead.
filterSNVImpact <- function(gr,
                            raw_score_threshold = NULL,
                            mut_score_threshold = NULL,
                            score_diff_threshold = NULL) {

  .Deprecated("filterVarImpact")

  if (!inherits(gr, "GRanges")) {
    stop("The input object 'gr' must be a GRanges object.")
  }

  if (length(gr) == 0) {
    stop("Input GRanges object 'gr' is empty.")
  }

  if (is.null(raw_score_threshold) &&
      is.null(mut_score_threshold) &&
      is.null(score_diff_threshold)) {
    stop("Error: All thresholds ",
         "('raw_score_threshold', 'mut_score_threshold',",
         " and 'score_diff_threshold') are NULL. ",
         "Please provide at least one.")
  }

  if (!is.null(raw_score_threshold) &&
      (!is.numeric(raw_score_threshold) ||
       raw_score_threshold < 0 ||
       raw_score_threshold > 4)) {
    stop("'raw_score_threshold' must be a numeric value between 0 and 4.")
  }

  if (!is.null(mut_score_threshold) &&
      (!is.numeric(mut_score_threshold) ||
       mut_score_threshold < 0 ||
       mut_score_threshold > 4)) {
    stop("'mut_score_threshold' must be a numeric value between 0 and 4.")
  }

  if (!is.null(score_diff_threshold) &&
      (!is.numeric(score_diff_threshold) ||
       score_diff_threshold >= 0 ||
       score_diff_threshold < -4)) {
    stop("'score_diff_threshold' must be a numeric value between -4 and 0.")
  }

  if (!is.null(raw_score_threshold)) {
    gr <- gr[abs(gr$G4.info.score) >= raw_score_threshold]
  }

  if (!is.null(mut_score_threshold)) {
    gr <- gr[abs(gr$mut.score) <= mut_score_threshold]
  }

  if (!is.null(score_diff_threshold)) {
    gr <- gr[gr$score.diff <= score_diff_threshold]
  }

  return(gr)
}

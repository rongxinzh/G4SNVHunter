#' Evaluate the Impact of Variants (SNVs, Indels, MNVs) on G4 Sequences
#'
#' This function evaluates the impact of variants (SNVs, indels, and MNVs) on
#' G4 formation.
#'
#' @param G4 A \code{GRanges} object representing the G4 regions. This object
#' must include a \code{sequence} metadata column containing the G4 sequences.
#' It should come from the output of the \code{G4HunterDetect} function with
#' \code{include_sequences = TRUE}.
#' @param variants A \code{GRanges} object representing the variants. This
#' object must include metadata columns for reference and alternative alleles.
#' @param ref_col A \code{character} string specifying the name of the column
#' in \code{variants} that contains the reference alleles.
#' Default is \code{"ref"}.
#' @param alt_col A \code{character} string specifying the name of the column
#' in \code{variants} that contains the alternative alleles.
#' Default is \code{"alt"}.
#' @param mode A \code{character} string indicating the mode of operation.
#' Set to \code{"s"} to evaluate the impact of individual variants on G4
#' regions one at a time (single-variant mode).
#' Set to \code{"m"} to assess the combined impact of multiple variants that
#' overlap the same G4 region within a sample (multi-variant Mode).
#' If using \code{"m"} mode, you must specify \code{sampleid_col}.
#' @param sampleid_col A \code{character} string specifying the name of the
#' column in \code{variants} that contains the sample IDs.
#' Required when \code{mode} is \code{"m"};
#' ignored if \code{mode} is \code{"s"}.
#'
#' @return A \code{GRanges} object with variant impact results:
#' \describe{
#'   \item{Mode "s" (single-variant mode )}{
#'     For each variant-G4 overlap:
#'     \itemize{
#'       \item Original G4 metadata (\code{G4.info.*})
#'       \item Variant information (\code{variant.info.*})
#'       \item Mutated sequence (\code{mutated.G4.seq})
#'       \item Annotated mutation sequence (\code{mutated.G4.anno.seq})
#'       \item New G4Hunter max_score (\code{mutated.max_score})
#'       \item Score difference (\code{score.diff})
#'     }
#'   }
#'   \item{Mode "m" (multi-variant mode)}{
#'     For each sample-G4 combination:
#'     \itemize{
#'       \item Original G4 metadata (\code{G4.info.*})
#'       \item Combined variant information (\code{variant.info.*})
#'       \item Mutated sequence with all variants incorporated
#'       (\code{mutated.G4.seq})
#'       \item Annotated mutation sequence (\code{mutated.G4.anno.seq})
#'       \item New G4Hunter max_score (\code{mutated.max_score})
#'       \item Score difference (\code{score.diff})
#'     }
#'   }
#' }
#'
#' @seealso \code{\link{G4HunterDetect}} for detecting the G4 sequences in a
#'                given \code{DNAStringSet} object.
#'          \code{\link{filterVarImpact}} for filtering out variants with
#'                significant impact.
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom IRanges start findOverlaps
#' @export
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
#' # Load variants
#' data(snv_gr)
#'
#' # 's' mode; single-variant mode ('s')
#' # evaluating each variant individually.
#' res_snv_s <- G4VarImpact(G4,
#'                          snv_gr,
#'                          ref_col = "ref",
#'                          alt_col = "alt")
#' print(res_snv_s)
#'
#' # 'm' mode; multi-variant mode ('m')
#' # evaluating the combined impact of variants on G4s.
#' # Grouped by the sample IDs specified in the 'sampleid_col' column.
#' res_snv_m <- G4VarImpact(G4,
#'                          snv_gr,
#'                          ref_col = "ref",
#'                          alt_col = "alt",
#'                          mode = "m",
#'                          sampleid_col = "sampleid"
#' )
#' print(res_snv_m)
G4VarImpact <- function(G4 = NULL,
                        variants = NULL,
                        ref_col = NULL,
                        alt_col = NULL,
                        mode = "s",
                        sampleid_col = NULL) {

  validateG4VarImpactInputs(G4,
                            variants,
                            ref_col,
                            alt_col,
                            mode,
                            sampleid_col)

  message("Start processing...")

  mcols(variants)[[alt_col]] <- vapply(
    mcols(variants)[[alt_col]],
    function(x) gsub(" ", "", paste0(as.character(x), collapse = ",")),
    character(1)
  )

  if (mode == "s") {

    variants <- splitAltVariants(variants, alt_col)

    overlaps <- findOverlaps(G4, variants)
    if (length(overlaps) == 0) {
      message("No overlaps found between variants and G4 regions")
      return(GRanges())
    }
    g4_hits <- queryHits(overlaps)
    variant_hits <- subjectHits(overlaps)

    var_G4 <- variantImpactG4ModeS(G4,
                                   variants,
                                   ref_col = ref_col,
                                   alt_col = alt_col,
                                   g4_hits, variant_hits)

  } else if (mode == "m") {

    checkSampleID(variants, sampleid_col)
    variants <- assignVarIDs(variants, ref_col, sampleid_col)
    variants <- splitAltVariants(variants, alt_col)

    overlaps <- findOverlaps(G4, variants)
    if (length(overlaps) == 0) {
      message("No overlaps found between variants and G4 regions")
      return(GRanges())
    }
    g4_hits <- queryHits(overlaps)
    variant_hits <- subjectHits(overlaps)

    var_G4 <- variantImpactG4ModeM(G4,
                                   variants,
                                   ref_col = ref_col,
                                   alt_col = alt_col,
                                   sampleid_col = sampleid_col,
                                   variantid_col = "G4SNVHunter_Var_ID",
                                   g4_hits, variant_hits)

  }

  metadata(var_G4) <-  metadata(G4)

  message("Processing finished!")

  return(var_G4)
}

#' validateG4VarImpactInputs
#'
#' This function validates the input arguments for the \code{G4VarImpact}
#' function.
#'
#' @param G4 A \code{GRanges} object representing the G4 regions. This object
#' must include a \code{sequence} metadata column containing the G4 sequences.
#' It should come from the output of the \code{G4HunterDetect} function with
#' \code{include_sequences = TRUE}.
#' @param variants A \code{GRanges} object representing the variants. This
#' object must include metadata columns for reference and alternative alleles.
#' @param ref_col A \code{character} string specifying the name of the column
#' in \code{variants} that contains the reference alleles.
#' Default is \code{"ref"}.
#' @param alt_col A \code{character} string specifying the name of the column
#' in \code{variants} that contains the alternative alleles.
#' Default is \code{"alt"}.
#' @param mode A \code{character} string indicating the mode of operation.
#' Set to \code{"s"} to evaluate the impact of individual variants on G4
#' regions one at a time (single-variant mode).
#' Set to \code{"m"} to assess the combined impact of multiple variants that
#' overlap the same G4 region within a sample (multi-variant mode).
#' If using \code{"m"} mode, you must specify \code{sampleid_col}.
#' @param sampleid_col A \code{character} string specifying the name of the
#' column in \code{variants} that contains the sample IDs.
#' Required when \code{mode} is \code{"m"};
#' ignored if \code{mode} is \code{"s"}.
#'
#' @return NULL
#'
#' @noRd
validateG4VarImpactInputs <- function(G4,
                                      variants,
                                      ref_col,
                                      alt_col,
                                      mode,
                                      sampleid_col) {

  if (is.null(G4) || !inherits(G4, "GRanges") || length(G4) == 0) {
    stop("'G4' must be a non-null, non-empty GRanges object.")
  }

  if (!"sequence" %in% colnames(mcols(G4))) {
    stop("'G4' must contain a 'sequence' column. Please use the ",
         "'G4HunterDetect' function with 'include_sequences = TRUE'.")
  }

  if (is.null(variants) ||
      !inherits(variants, "GRanges") ||
      length(variants) == 0) {
    stop("'variants' must be a non-null, non-empty GRanges object.")
  }

  if (is.null(ref_col) || nchar(ref_col) == 0) {
    stop("'ref_col' cannot be NULL or empty")
  }

  if (is.null(alt_col) || nchar(alt_col) == 0) {
    stop("'alt_col' cannot be NULL or empty.")
  }

  if (!mode %in% c("s", "m")) {
    stop("Invalid mode. Choose 's' for single-variant mode or ",
         "'m' for multi-variant mode.")
  }

  if (!alt_col %in% colnames(mcols(variants))) {
    stop("variants must contain a column named '", alt_col, "' in mcols.")
  }

  if (mode == "m" &&
      (is.null(sampleid_col) ||
       !sampleid_col %in% colnames(mcols(variants)))) {
    stop("When mode is 'm', 'sampleid_col' must be specified and present ",
         "in 'variants'.")
  }

  if (mode == "s" && !is.null(sampleid_col)) {
    stop("When mode is 's', 'sampleid_col' should not be specified.")
  }

}

#' assignVarIDs
#'
#' Assign IDs to the \code{GRanges} variants
#'
#' @param gr A \code{GRanges} object containing variants.
#' @param ref_col Name of reference allele column.
#' @param sampleid_col A string specifying the column for sample identifier.
#'
#' @return A \code{GRanges} object with an additional ID metadata column
#' \code{G4SNVHunter_Var_ID}.
#'
#' @noRd
assignVarIDs <- function(gr, ref_col = NULL, sampleid_col = NULL) {

  gr_df <- as.data.frame(gr)

  sample_id_vec <- gr_df[[sampleid_col]]
  ref_vec <- gr_df[[ref_col]]

  gr_df$variant_key <- paste(sample_id_vec, gr_df$seqnames,
                             gr_df$start, ref_vec, sep = "_")

  unique_keys <- unique(gr_df$variant_key)
  key_to_id <- setNames(paste0("ID_", seq_along(unique_keys)), unique_keys)
  gr_df$id <- key_to_id[gr_df$variant_key]

  mcols(gr)$G4SNVHunter_Var_ID <- gr_df$id

  return(gr)
}


#' variantImpactG4ModeS
#'
#' Processes variant impacts on G4 formation in single-variant mode.
#'
#' @param G4  A \code{GRanges} object containing G4s.
#' @param variants  A \code{GRanges} object containing variants.
#' @param ref_col Name of reference allele column.
#' @param alt_col Name of alternate allele column.
#' @param g4_hits An integer vector of row indices in G4 object that overlap
#' with variants.
#' @param variant_hits An integer vector of row indices in variant object that
#' overlap with G4s.
#'
#' @return A \code{GRanges} object with mutated G4 sequences and associated
#' variant information.
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @noRd
variantImpactG4ModeS <- function(G4 = NULL,
                                 variants = NULL,
                                 ref_col = NULL,
                                 alt_col = NULL,
                                 g4_hits, variant_hits) {

  G4_filtered <- G4[g4_hits]

  results <- lapply(seq_along(G4_filtered), function(i) {
    G4_item <- G4_filtered[i]
    affected_var <- variants[variant_hits[i]]
    mut_result <- applyVariantsToG4(G4_item, affected_var, ref_col, alt_col)

    mut_G4_pred <- mutG4MaxScore(
      mut_seq = as.character(mut_result$mut_g4_seq),
      G4_orig_start = start(G4_item),
      threshold = metadata(G4_item)$threshold,
      window_size = metadata(G4_item)$window_size,
      strands = metadata(G4_item)$strands
    )

    GRanges(
      G4_item,
      G4.info = mcols(G4_item),
      variant.info.ranges.start = start(affected_var),
      variant.info.ranges.end = end(affected_var),
      variant.info.ranges.width = width(affected_var),
      variant.info = mcols(affected_var),
      mutated.G4.seq = mut_G4_pred$mut_G4_seq,
      mutated.G4.anno.seq = as.character(mut_result$mut_g4_anno_seq),
      mutated.max_score = mut_G4_pred$max_score,
      score.diff =
        signif(abs(mut_G4_pred$max_score) - abs(G4_item$max_score), 3)
    )

  })

  final_var_G4 <- do.call(c, results)

  return(final_var_G4)
}

#' variantImpactG4ModeM
#'
#' Processes variant impacts on G4 formation in multiple variant mode.
#'
#' @param G4 A \code{GRanges} object containing G4s.
#' @param variants A \code{GRanges} object containing variants.
#' @param ref_col Name of reference allele column.
#' @param alt_col Name of alternate allele column.
#' @param sampleid_col Name of sample ID column.
#' @param variantid_col  Name of variant ID column.
#' @param g4_hits An integer vector of row indices in G4 object that overlap
#' with variants.
#' @param variant_hits An integer vector of row indices in variant object that
#' overlap with G4s.
#'
#' @return A \code{GRanges} object with mutated G4 sequences and associated
#' variant information.
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @noRd
variantImpactG4ModeM <- function(G4 = NULL,
                                 variants = NULL,
                                 ref_col = NULL,
                                 alt_col = NULL,
                                 sampleid_col = NULL,
                                 variantid_col = NULL,
                                 g4_hits, variant_hits) {

  G4_filtered <- G4[unique(g4_hits)]

  results <- lapply(seq_along(G4_filtered), function(i) {
    G4_item <- G4_filtered[i]
    idx <- unique(g4_hits)[i]

    affected_var <- variants[variant_hits[g4_hits == idx]]

    #var_groups <- split(affected_var, mcols(affected_var)[[sampleid_col]])
    var_groups <- processVariants(affected_var,
                                  alt_col,
                                  sampleid_col,
                                  variantid_col)

    sample_results <- lapply(names(var_groups), function(sample) {
      var_s <- var_groups[[sample]]

      mut_result <- applyVariantsToG4(G4_item, var_s, ref_col, alt_col)

      mut_G4_pred <- mutG4MaxScore(
        mut_seq = as.character(mut_result$mut_g4_seq),
        G4_orig_start = start(G4_item),
        threshold = metadata(G4_item)$threshold,
        window_size = metadata(G4_item)$window_size,
        strands = metadata(G4_item)$strands
      )

      #meta_cols <- setdiff(colnames(mcols(var_s)),
      #                     c(sampleid_col, variantid_col))

      meta_cols <- setdiff(colnames(mcols(var_s)), c(variantid_col))

      variant_info <- do.call(DataFrame, c(
        list(variant.info.ranges.start = paste(as.character(start(var_s)),
                                               collapse = ";"),
             variant.info.ranges.end = paste(as.character(end(var_s)),
                                               collapse = ";"),
             variant.info.ranges.width = paste(as.character(width(var_s)),
                                               collapse = ";")),
        setNames(lapply(meta_cols, function(col) paste(mcols(var_s)[[col]],
                                                       collapse = ";")),
                 paste0("variant.info.", meta_cols))
      ))

      tmp_gr <- GRanges(
        G4_item,
        G4.info = mcols(G4_item)#,
        #variant.ids = paste(mcols(var_s)[[variantid_col]], collapse = ";"),
        #sample.ids = sub("_[^_]+$", "", sample)
      )

      mcols(tmp_gr) <- c(mcols(tmp_gr), variant_info, DataFrame(
        mutated.G4.seq = mut_G4_pred$mut_G4_seq,
        mutated.G4.anno.seq = as.character(mut_result$mut_g4_anno_seq),
        mutated.max_score = mut_G4_pred$max_score,
        score.diff =
          signif(abs(mut_G4_pred$max_score) - abs(G4_item$max_score), 3)
      ))

      sample_colname <- paste0("variant.info.", sampleid_col)
      sample_id_all <- mcols(tmp_gr)[[sample_colname]]
      cleaned_sample_id <- paste(unique(
        strsplit(sample_id_all, ";", fixed = TRUE)[[1]]), collapse = ";")
      mcols(tmp_gr)[[sample_colname]] <- cleaned_sample_id

      return(tmp_gr)

    })

    do.call(c, sample_results)
  })

  final_var_G4 <- do.call(c, results)

  return(final_var_G4)
}

#' applyVariantsToG4
#'
#' Introduces variants into a G4 sequence, and generate both the G4 mutated
#' sequence and the annotated sequence indicating the mutations.
#'
#' @param g4 A \code{GRanges} object containing a G4 sequence.
#' @param variants A \code{GRanges} object containing variants overlapping the
#' G4.
#' @param ref_col Name of reference allele column.
#' @param alt_col Name of alternate allele column.
#'
#' @return A \code{list} containing mutated sequence and annotated sequence
#' with mutations in \code{[ref>alt]} format.
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @noRd
applyVariantsToG4 <- function(g4,
                              variants,
                              ref_col = NULL,
                              alt_col = NULL) {

  end(variants) <- start(variants) + nchar(mcols(variants)[[ref_col]]) - 1

  g4_seq <- g4$sequence[1]
  g4_start <- start(g4)[1]
  g4_end <- end(g4)[1]

  variants <- variants[order(start(variants)), ]
  ref_values <- mcols(variants)[[ref_col]]
  #alt_values <- mcols(variants)[[alt_col]]
  alt_values <- gsub(" ", "", mcols(variants)[[alt_col]])
  is_insert <- ref_values == "-"

  rel_starts <- pmax(1, start(variants) - g4_start + 1)
  rel_ends <- pmin(end(variants) - g4_start + 1, nchar(g4_seq))

  cut_points <- sort(unique(c(1, rel_starts, rel_ends + 1, nchar(g4_seq) + 1)))
  seg_starts <- cut_points[-length(cut_points)]
  seg_ends <- cut_points[-1] - 1
  segments <- substring(g4_seq, seg_starts, seg_ends)

  is_var_segment <- rowSums(outer(seg_starts, rel_starts, `>=`) &
                              outer(seg_starts, rel_ends, `<=`)) > 0
  match_idx <- max.col(outer(seg_starts, rel_starts, `>=`) &
                         outer(seg_starts, rel_ends, `<=`),
                       "first")
  match_idx[!is_var_segment] <- NA

  var_idx <- match_idx[is_var_segment]
  #pos_refs <-
  #  str_sub(g4_seq, seg_starts[is_var_segment], seg_starts[is_var_segment])
  pos_refs <-
    substr(g4_seq, seg_starts[is_var_segment], seg_starts[is_var_segment])
  is_insert_segment <- is_insert[var_idx]

  modified_segments <- segments
  modified_segments[is_var_segment] <- ifelse(
    is_insert_segment,
    paste0(pos_refs, alt_values[var_idx]),
    alt_values[var_idx]
  )

  mut_anno_seq <- segments
  mut_anno_seq[is_var_segment] <- ifelse(
    is_insert_segment,
    paste0(pos_refs, "[->", alt_values[var_idx], "]"),
    paste0("[", ref_values[var_idx], ">",
           alt_values[var_idx], "]")
  )

  list(
    mut_g4_seq =
      gsub("-", "", paste0(modified_segments, collapse = ""), fixed = TRUE),
    mut_g4_anno_seq = paste0(mut_anno_seq, collapse = "")
  )
}

#' mutG4MaxScore
#'
#' Recalculates the maximum G4Hunter score for a mutated G4 sequence.
#'
#' @param mut_seq A \code{character} string of the mutated G4 sequence.
#' @param G4_orig_start An \code{integer} specifying the original start
#' position of the G4.
#' @param threshold A \code{numeric} value specifying the G4Hunter score
#' threshold used for prediction.
#' @param window_size An \code{integer} specifying the window size used in
#' G4Hunter prediction.
#' @param strands A \code{character} string indicating the strand(s) to
#' predict on.
#'
#' @return A \code{list} containing the maximum G4Hunter score identified
#' after mutation and the G4-forming sequence corresponding to that score.
#'
#' @import S4Vectors
#'
#' @noRd
mutG4MaxScore <- function(mut_seq = "",
                          G4_orig_start = 25,
                          threshold = 1.5,
                          window_size = 25,
                          strands = "b"
                             ) {

  if (nchar(mut_seq) < window_size) {
    diff_nuc <- window_size - nchar(mut_seq)
    left_pad <- ifelse(ceiling(diff_nuc / 2) > G4_orig_start,
                       G4_orig_start - 1, ceiling(diff_nuc / 2))
    right_pad <- diff_nuc - left_pad
  } else {
    left_pad <- 0
    right_pad <- 0
  }

  tmp_G4 <- G4HunterAnalysis(paste0(strrep("N", left_pad),
                                    as.character(mut_seq),
                                    strrep("N", right_pad)),
                             "CHR", threshold, window_size,
                             include_sequences = T, strands)

  if (length(tmp_G4) == 0) {
    # the G4 score fell below the threshold after the mutation
    # Re-predict at threshold 0 to find the most high-scoring G4s (unstable)
    tmp_G4 <- G4HunterAnalysis(paste0(strrep("N", left_pad),
                                      as.character(mut_seq),
                                      strrep("N", right_pad)),
                               "CHR", 0, window_size,
                               include_sequences = T, strands)
  }

  tmp_pred <- list()
  if (length(tmp_G4) == 0) {
    # No any Gs/Cs can be found
    tmp_pred$max_score <- 0
    tmp_pred$mut_G4_seq <- mut_seq
  } else {
    idx <- which.max(abs(mcols(tmp_G4)$max_score))[1]
    tmp_G4 <- tmp_G4[idx]
    tmp_pred$max_score <- tmp_G4$max_score
    tmp_pred$mut_G4_seq <- tmp_G4$sequence
  }

  return(tmp_pred)
}

#' splitAltVariants
#'
#' Split multi-allelic variants into separate ones.
#'
#' @param gr A \code{GRanges} object containing variant data.
#' @param alt_col A \code{character} string specifying the column name for
#' alternative alleles.
#'
#' @return A \code{GRanges} object, where each record represents a
#' non-multiallelic variant
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest
#' @importFrom data.table :=
#' @importFrom magrittr %>%
#'
#' @noRd
splitAltVariants <- function(gr, alt_col = NULL) {

  df <- as.data.frame(gr)

  if (!any(grepl(",", df[[alt_col]], fixed = TRUE))) {
    return(gr)
  } else {
    df_split <- df %>%
      mutate(!!alt_col := strsplit(as.character(.data[[alt_col]]), ",")) %>%
      unnest(!!alt_col) %>% mutate(!!alt_col := as.character(.data[[alt_col]]))

    metadata_cols <- df_split[, !colnames(df_split) %in%
                                c("seqnames", "start", "end",
                                  "width", "strand")]
    new_gr <- GRanges(
      seqnames = df_split$seqnames,
      ranges = IRanges(start = df_split$start, end = df_split$end),
      strand = df_split$strand,
      metadata_cols)

    return(new_gr)
  }
}

#' checkSampleID
#'
#' Validate if sample IDs are missing.
#'
#' @param gr A \code{GRanges} object containing variant data.
#' @param id_col A \code{character} string specifying the column name for
#' sample ID.
#'
#' @return No return value. Stops execution if invalid IDs
#' (missing or empty) are found.
#'
#' @noRd
checkSampleID <- function(gr, id_col = NULL) {
  id_values <- mcols(gr)[[id_col]]

  bad_idx <- which(is.na(id_values) | id_values == "")

  if (length(bad_idx) > 0) {
    message("Invalid sample ID first found at:\n")
    print(gr[bad_idx[1]])
    stop("Sample ID contains missing (NA) or empty values.")
  }
}

#' processVariants
#'
#' Generate all possible multi-allelic combinations within the G4 region for
#' each sample.
#'
#' @param gr A \code{GRanges} object containing variant data.
#' @param alt_col A \code{character} string specifying the name of the column
#' in \code{gr} that contains the alternative alleles.
#' @param sampleid_col A \code{character} string specifying the name of the
#' column in \code{gr} that contains the sample IDs.
#' @param variantid_col A \code{character} string specifying the name of the
#' column in \code{gr} that contains variant IDs.
#'
#' A \code{GRangesList}, where each element represents a combination of
#' multiallelic variants for a specific sample.
#'
#' @importFrom data.table CJ
#' @importFrom GenomicRanges GRangesList
#'
#' @noRd
processVariants <- function(gr, alt_col, sampleid_col, variantid_col) {

  samples <- split(gr, mcols(gr)[[sampleid_col]])

  result_list <- lapply(names(samples), function(sample_name) {
    sample_data <- samples[[sample_name]]

    variant_options <- split(mcols(sample_data)[[alt_col]],
                             mcols(sample_data)[[variantid_col]])
    variant_options <- lapply(variant_options, unique)

    all_combinations <- do.call(CJ, c(variant_options, unique = TRUE))

    gr_list <- lapply(seq_len(nrow(all_combinations)), function(i) {
      new_gr <-
        GRangesList(lapply(names(variant_options), function(variant_id) {
        #alt_value <- all_combinations[i, ..variant_id]
        alt_value <- all_combinations[i, get(variant_id)]
        pos_data <-
          sample_data[mcols(sample_data)[[variantid_col]] == variant_id][1]
        mcols(pos_data)[[alt_col]] <- as.character(alt_value)
        pos_data
      }))
      unlist(new_gr, use.names = FALSE)
    })

    names(gr_list) <- paste0(sample_name, "_", seq_along(gr_list))

    gr_list
  })

  gr_split <- do.call(c, result_list)

  return(gr_split)
}

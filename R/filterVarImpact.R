#' Filter Variant Impact G4 Object Based on User-Defined Thresholds
#'
#' This function filters the G4 (\code{GRanges}) object returned by the
#' \code{G4VarImpact} function based on user-defined thresholds for the
#' \code{G4.info.max_score}, \code{mutated.max_score}, and \code{score.diff}
#' columns.
#' Users are not required to specify all three thresholds;
#' however, at least one must be provided.
#'
#' @param mut_G4 A \code{GRanges} object returned by the \code{G4VarImpact}
#' function.
#' @param raw_score_threshold A positive numeric value (no greater than 4)
#' used as the threshold for the absolute value of \code{G4.info.max_score}.
#' G4s with an absolute G4Hunter max_score exceeding this threshold will be
#' retained.
#' If \code{NULL}, this threshold is not applied.
#' @param mut_score_threshold A positive numeric value (no greater than 4)
#' used as the threshold for the absolute value of \code{mutated.max_score}.
#' Mutated G4s with an absolute G4Hunter max_score below this threshold will be
#' retained.
#' If \code{NULL}, this threshold is not applied.
#' @param score_diff_threshold A negative numeric value (no less than -4)
#' used as the threshold for \code{score.diff}. G4s with a decrease in
#' G4Hunter score greater than this threshold after mutation will be retained.
#' If \code{NULL}, this threshold is not applied.
#'
#' @seealso \code{\link{G4VarImpact}} for assessing the impact of variants on
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
#' res_snp <- G4VarImpact(G4, snp_gr, ref_col = "ref", alt_col = "alt")
#' filtered_var_eff <- filterVarImpact(res_snp,
#'   mut_score_threshold = 1.2,
#'   score_diff_threshold = -0.2
#' )
#' print(filtered_var_eff)
filterVarImpact <- function(mut_G4,
                            raw_score_threshold = NULL,
                            mut_score_threshold = NULL,
                            score_diff_threshold = NULL) {

  if (!inherits(mut_G4, "GRanges")) {
    stop("The input object 'mut_G4' must be a GRanges object.")
  }

  if (length(mut_G4) == 0) {
    stop("Input GRanges object 'mut_G4' is empty.")
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
    mut_G4 <- mut_G4[abs(mut_G4$G4.info.max_score) >= raw_score_threshold]
  }

  if (!is.null(mut_score_threshold)) {
    mut_G4 <- mut_G4[abs(mut_G4$mutated.max_score) <= mut_score_threshold]
  }

  if (!is.null(score_diff_threshold)) {
    mut_G4 <- mut_G4[mut_G4$score.diff <= score_diff_threshold]
  }

  if (!is.null(raw_score_threshold)) {
    metadata(mut_G4)$raw_score_threshold <- raw_score_threshold
  }

  if (!is.null(mut_score_threshold)) {
    metadata(mut_G4)$mut_score_threshold <- mut_score_threshold
  }

  if (!is.null(score_diff_threshold)) {
    metadata(mut_G4)$score_diff_threshold <- score_diff_threshold
  }

  return(mut_G4)
}

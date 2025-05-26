#' Plot the Impact of Variants on G4 Formation
#'
#' This function generates two panels to visualize the impact of variants on
#' G4 formation:
#' \enumerate{
#'   \item A 2D density plot of absolute maximum G4Hunter scores before and
#'   after mutation.
#'   \item A density plot showing the distribution of score differences.
#' }
#'
#' @param gr A \code{GRanges} object containing G4Hunter scores and
#' associated metadata.
#' Required metadata columns include:
#' \itemize{
#'   \item \code{G4.info.max_score}: Numeric vector of
#'   original maximum G4Hunter scores.
#'   \item \code{mutated.max_score}: Numeric vector of
#'   maximum G4Hunter scores after mutation.
#'   \item \code{score.diff}: Numeric vector of score differences.
#' }
#'
#' @param p_colors A character vector of five colors used for plotting.
#' The vector must include:
#' \describe{
#'   \item{\code{"fill_2d"}}{Fill Color for the 2D density plot.}
#'   \item{\code{"diagonal"}}{Color for the diagonal line in the 2D
#'   density plot.}
#'   \item{\code{"density_fill"}}{Fill color for the density plot.}
#'   \item{\code{"density_line"}}{Line color for the density plot.}
#'   \item{\code{"vline"}}{Color for the vertical line at zero
#'   in the density plot.}
#' }
#'
#' @return A \code{ggplot} object combining two subplots:
#'   \enumerate{
#'     \item 2D density plot of the absolute maximum G4Hunter scores
#'     before and after mutation.
#'     \item Density plot of the maximum G4Hunter score differences.
#'   }
#'
#' @seealso \code{\link{G4VarImpact}} for computing the effects of variants
#' on G4 formation,
#' and \code{\link{filterVarImpact}} for screening G4s whose formation
#' propensity may be impaired by variants.
#'
#' @import ggplot2
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom ggplot2 scale_color_viridis_c
#' @importFrom ggdensity geom_hdr
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
#' # Load SNVs
#' data(snp_gr)
#'
#' res_snp <- G4VarImpact(G4, snp_gr, ref_col = "ref", alt_col = "alt")
#' plotVarImpact(res_snp)
plotVarImpact <- function(gr,
                          p_colors = c(
                            fill_2d = "#376597",
                            diagonal = "#b22d2d",
                            density_fill = "#6ca4d6",
                            density_line = "#2d69b0",
                            vline = "#1f77b4"
                          )) {
  required_names <- c("fill_2d", "diagonal", "density_fill",
                      "density_line", "vline")
  if (!all(required_names %in% names(p_colors))) {
    stop("'p_colors' must be a named vector with elements: ",
         paste(required_names, collapse = ", "))
  }

  df <- as.data.frame(gr)

  max_range <- max(abs(df$G4.info.max_score), abs(df$mutated.max_score))
  min_range <- min(abs(df$G4.info.max_score), abs(df$mutated.max_score))

  p1 <- ggplot(df, aes(x = abs(.data$G4.info.max_score),
                       y = abs(.data$mutated.max_score))) +
    geom_hdr(fill = p_colors["fill_2d"]) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.8,
                linetype = "dashed", color = p_colors["diagonal"]) +
    scale_x_continuous(limits = c(min_range, max_range), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min_range, max_range), expand = c(0, 0)) +
    labs(
      title = "2D Density Plot of\nMax G4Hunter Scores",
      x = "Original Score (absolute)\n",
      y = "Mutated Score (absolute)"
    ) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )

  p2 <- ggplot(df, aes(x = .data$score.diff)) +
    geom_density(fill = p_colors["density_fill"],
                 color = p_colors["density_line"],
                 alpha = 0.5, linewidth = 0.8) +
    geom_vline(xintercept = 0, color = p_colors["vline"],
               linewidth = 0.8, linetype = "dashed") +
    labs(
      title = "Distribution of Changes in\nMax G4Hunter Scores",
      x = "Change in Max G4Hunter Score\nabs(Mutated) - abs(Original)",
      y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )

  combined_plot <- plot_grid(p1, p2, ncol = 2,
                             labels = c("a", "b"),
                             label_size = 14)

  title <- ggdraw() +
    draw_label(
      "Mutational Impact on G4 Formation",
      fontface = 'bold',
      size = 14,
      x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(title, combined_plot,
                          ncol = 1, rel_heights = c(0.1, 1))

  return(final_plot)
}

#' Sequence Logo of the G4 Sequence and Its Mutant Counterpart
#'
#' This function visualizes and compares the G4 and its mutant counterpart
#' using sequence logos.
#'
#' @param mut_G4 A \code{GRanges} object containing the mutant G4.
#' If multiple G4s are provided, only the first is plotted.

#' @param keep_gstrand Logical. If \code{TRUE} (default), negative
#' strand (C-rich) G4 will be displayed in its complementary reversed form.
#'
#' @return A ggplot2-based object containing the G4 sequence logo.
#'
#' @import ggplot2
#' @import ggseqlogo
#' @importFrom grid arrow unit
#' @importFrom cowplot plot_grid
#' @importFrom GenomicRanges strand
#'
#' @examples
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
#' # Load SNVs
#' data(snp_gr)
#'
#' res_snp <- G4VarImpact(G4, snp_gr, ref_col = "ref", alt_col = "alt")
#' plotImpactedG4(res_snp[1])
#'
#' @export
plotImpactedG4 <- function(mut_G4 = NULL, keep_gstrand = TRUE) {

  if (!inherits(mut_G4, "GRanges")) {
    stop("'mut_G4' must be a GRanges object")
  }

  if (length(mut_G4) == 0) {
    stop("Input GRanges object 'mut_G4' is empty.")
  }

  if (length(mut_G4) > 1) {
    warning("Multiple records detected in the 'mut_G4' object. ",
            "Only the first record will be plotted.\n",
            "Please use a loop function to plot all records.")
  }

  mut_G4 <- mut_G4[1]

  version =
    ifelse("mutated.max_score" %in% colnames(mcols(mut_G4)),
           "new", "initial")

  if (version == "initial") {
    stop("This function requires G4SNVHunter version > 1.0.0.
         Please upgrade your package.")
  }

  original_score <- mut_G4$G4.info.max_score[1]
  mutated_score <- mut_G4$mutated.max_score[1]

  label <- paste0("Original max G4Hunter score: ", original_score,
                  "\nMutated max G4Hunter score: ", mutated_score)

  mutated_seq <- toupper(mut_G4$mutated.G4.anno.seq[1])
  mutated_seqs <- parseVarSeq(mutated_seq)

  if (as.character(strand(mut_G4)) == "-" & keep_gstrand) {
    mutated_seqs <- c(comRevBase(mutated_seqs[1]), comRevBase(mutated_seqs[2]))
  }

  original_seq_char <- strsplit(mutated_seqs[1], "")[[1]]
  mutated_seq_char <- strsplit(mutated_seqs[2], "")[[1]]

  aln_ori_mat <- data.frame(char = original_seq_char,
                            group = rep(paste0("Ref seq\n(Max score:",
                                               ifelse(keep_gstrand,
                                                      abs(original_score),
                                                      original_score), ")"),
                                        length(original_seq_char)),
                            position = seq_along(original_seq_char),
                            mut = rep('no', length(original_seq_char)))

  aln_mut_mat <- data.frame(char = mutated_seq_char,
                            group = rep(paste0("Mutated seq\n(Max score:",
                                               ifelse(keep_gstrand,
                                                      abs(mutated_score),
                                                      mutated_score), ")"),
                                        length(mutated_seq_char)),
                            position = seq_along(mutated_seq_char),
                            mut = rep('no', length(mutated_seq_char)))

  diff_positions <- which(original_seq_char != mutated_seq_char)
  aln_ori_mat$mut[diff_positions] <- "yes"
  aln_mut_mat$mut[diff_positions] <- "yes"

  aln_mat <- rbind(aln_ori_mat, aln_mut_mat)
  fig <- plotSeqlogo(mutated_seqs, label = label, aln_mat = aln_mat)

  return(fig)

}

#' parseVarSeq
#'
#' This function splits a variant sequence into original and
#' mutated sequences.
#'
#' @param variant_seq A \code{character} string representing a G4 sequence
#' with embedded mutation annotation(s)
#'
#' @return A \code{character} vector of length 2 containing the original and
#' mutated sequences.
#' @noRd
parseVarSeq <- function(variant_seq) {
  seq1 <- character(0)
  seq2 <- character(0)

  tokens <- regmatches(variant_seq,
                       gregexpr("\\[[^]]*\\]|[^[]+",
                                variant_seq))[[1]]

  for (token in tokens) {
    if (grepl("^\\[.*\\]$", token)) {
      content <- substr(token, 2, nchar(token)-1)

      if (grepl(">", content)) {
        parts <- strsplit(content, ">")[[1]]
        left <- parts[1]
        right <- parts[2]

        len_left <- nchar(left)
        len_right <- nchar(right)

        if (len_left > len_right) {
          right <-
            paste0(right, paste(rep("-", len_left - len_right), collapse=""))
        } else if (len_right > len_left) {
          left <-
            paste0(left, paste(rep("-", len_right - len_left), collapse=""))
        }

        seq1 <- c(seq1, left)
        seq2 <- c(seq2, right)
      }

    } else {
      seq1 <- c(seq1, token)
      seq2 <- c(seq2, token)
    }
  }

  result1 <- paste(seq1, collapse="")
  result2 <- paste(seq2, collapse="")

  return(c(result1, result2))
}

#' comRevBase
#'
#' This function returns the reverse complement of a nucleotide sequence.
#'
#' @param seq A \code{character} string representing a DNA or RNA sequence.
#'
#' @return A \code{character} string of the reverse complement sequence.
#'
#' @noRd
comRevBase <- function(seq) {

  if (grepl("U", seq, ignore.case = TRUE)) {
    tmp_seq <- chartr("AUGCTaugct", "UACGAuacga", seq)  # RNA
  } else {
    tmp_seq <- chartr("ATGCUatgcu", "TACGAtacga", seq)  # DNA
  }

  tmp_seq <- paste(rev(strsplit(tmp_seq, "")[[1]]), collapse = "")
  return(tmp_seq)
}

#' plotSeqlogo
#'
#' This function generates a sequence logo plot, with mutations highlighted.
#'
#' @param seqs A \code{character} vector of original and mutated G4 sequences.
#' @param label A \code{character} string of the seqplot title.
#' @param aln_mat A \code{data.frame} containing variant
#' annotation information.
#' @param alphabet A \code{character} vector of nucleotide types to be counted.
#'
#' @return A ggplot object containing the sequence logo plot.
#' @noRd
plotSeqlogo <- function(seqs,
                        label = "",
                        aln_mat = "",
                        alphabet = c("A", "T", "G", "C", "N", "U", "-")) {
  seq_mat <- do.call(rbind, strsplit(seqs, split = ""))

  ppm <- vapply(
    X = seq_len(ncol(seq_mat)),
    FUN = function(i) {
      table(factor(seq_mat[, i], levels = alphabet))
    },
    FUN.VALUE = numeric(length(alphabet))
  )
  rownames(ppm) <- alphabet
  ppm <- ppm / nrow(seq_mat)

  num_nonzero_per_col <-
    vapply(seq_len(ncol(ppm)), function(i) sum(ppm[, i] > 0), numeric(1))
  mut_cols <- which(num_nonzero_per_col > 1)

  if (length(mut_cols) > 0) {
    mut_blocks <- split(mut_cols, cumsum(c(1, diff(mut_cols) != 1)))
  } else {
    mut_blocks <- list()
  }

  p <- ggseqlogo(ppm, method = "custom",
                 seq_type = ifelse(any(grepl("U", seqs, ignore.case = TRUE)),
                                   "rna", "dna")) +
    ylab("probability") +
    labs(title = label) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))

  for (block in mut_blocks) {
    p <- p + annotate("rect",
                      xmin = min(block) - 0.5,
                      xmax = max(block) + 0.5,
                      ymin = -0.05,
                      ymax = 1.05,
                      alpha = 0.15,
                      col='black',
                      fill = "yellow")
  }

  segment_df <- data.frame(
    x = mut_cols,
    xend = mut_cols,
    y = -0.2,
    yend = -0.1
  )

  p <- p + geom_segment(
    data = segment_df,
    aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
    color = '#b22d2d',
    linewidth = 0.8,
    arrow = arrow(type = "closed", length = unit(0.1, "inches")))

  p2 <- ggplot(aln_mat, aes(.data$position, .data$group)) +
    geom_text(aes(label = .data$char, color = .data$mut, size = .data$mut)) +
    scale_x_continuous(breaks = seq_len(10), expand = c(0.105, 0)) +
    xlab('') + ylab('') +
    scale_color_manual(values = c('black', '#b22d2d')) +
    scale_size_manual(values = c(5, 6)) +
    theme_logo() +
    theme(legend.position = 'none', axis.text.x = element_blank())

  return(plot_grid(p, p2, ncol = 1, align = 'v'))
}

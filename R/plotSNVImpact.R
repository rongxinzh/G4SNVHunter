#' Plot the Impact of SNVs on G4Hunter Scores
#'
#' This function generates two plots for visualizing the impact of SNVs on
#' G4 stability:
#' 1. A scatter plot with density shading comparing the original G4Hunter score
#' and the mutant G4Hunter score.
#' 2. A density plot showing the distribution of score changes between the
#' original and mutant G4 sequences.
#'
#' @param gr A \code{GRanges} object containing the G4Hunter scores and
#' associated metadata. The object should include the following columns:
#' \itemize{
#'   \item \code{G4.info.score}: Numeric vector of the original G4Hunter
#'   scores.
#'   \item \code{mut.score}: Numeric vector of the G4Hunter scores after
#'   mutation.
#'   \item \code{score.diff}: Numeric vector of the differences between the
#'   original and mutant G4Hunter scores.
#' }
#'
#' @param p_colors A vector of four colors used for plotting.
#'
#' @return A combined plot (using \code{plot_grid}) containing two subplots:
#'         - Density scatter plot comparing original vs mutant G4Hunter scores.
#'         - Density plot of the G4Hunter score differences.
#'
#' @seealso \code{\link{SNVImpactG4}} for evaluating the impact of SNVs on G4
#' stability, and \code{\link{filterSNVImpact}} for filtering G4s that are
#' significantly affected by SNVs.
#'
#' @import ggplot2
#' @import ggpointdensity
#' @import viridis
#' @import cowplot
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
#' res_snp <- SNVImpactG4(G4, snp_gr, alt_col = "alt")
#' plotSNVImpact(res_snp)

plotSNVImpact <- function(gr,
                          p_colors = c("#b22d2d", "#6ca4d6",
                                       "#2d69b0", "#1f77b4")) {
  if (length(p_colors) != 4) {
    stop("'p_colors' must be a vector of exactly four colors.")
  }

  df <- as.data.frame(gr)

  max_range <- max(abs(df$G4.info.score), abs(df$mut.score))
  min_range <- min(abs(df$G4.info.score), abs(df$mut.score))

  p1 <- ggplot(data = df, mapping = aes(x = abs(.data$G4.info.score),
                                        y = abs(.data$mut.score))) +
    geom_pointdensity() +
    scale_color_viridis() +
    geom_abline(intercept = 0, slope = 1, linewidth = 1,
                linetype = "dashed", color = p_colors[1]) +
    scale_x_continuous(limits = c(min_range, max_range)) +
    scale_y_continuous(limits = c(min_range, max_range)) +
    labs(
      title = "Density Plot of G4Hunter Score",
      x = "Original G4 Score",
      y = "Mutant G4 Score"
    ) +
    coord_fixed() +
    theme_bw()

  p2 <- ggplot(df, aes(x = .data$score.diff)) +
    geom_density(fill = p_colors[2], alpha = 0.5, color = p_colors[3]) +
    geom_vline(xintercept = 0, color = p_colors[4],
               linewidth = 1, linetype = "dashed") +
    labs(
      title = "Distribution of G4Hunter Score Changes",
      x = "Score Changes",
      y = "Density Estimate"
    ) +
    theme_bw()

  combined_plot <- plot_grid(p1, p2, ncol = 2, labels = "AUTO", align = 'v')

  title <- ggdraw() +
    draw_label(
      "Changes in G4Hunter Score",
      fontface = 'bold',
      x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(title, combined_plot,
                          ncol = 1, rel_heights = c(0.1, 1))

  return(final_plot)
}

#' Visualize the variants in G4 sequence
#'
#' This function plot sequence logos to visualize sequence variants caused by
#' SNVs or SNPs, with the location of the variants highlighted by rectangles
#' and arrows.
#'
#' @param filtered_gr A \code{GRanges} object containing sequence data and
#' G4Hunter scores. The object must have metadata columns named
#' \code{G4.info.score}, \code{mut.score}, \code{G4.info.sequence}, and
#' \code{mut.G4.seq}.
#' @param ncol An integer specifying the number of columns in the output plot
#' grid. Default is 1.
#'
#' @return A plot that displays the grid of sequence logos, showing the
#' differences between the original and mutated sequences.
#'
#' @seealso \code{\link{SNVImpactG4}} for evaluating the impact of SNVs on G4
#' stability, and \code{\link{filterSNVImpact}} for filtering G4s that are
#' significantly affected by SNVs.
#'
#' @import ggplot2
#' @import ggseqlogo
#' @import cowplot
#' @export
#'
#' @examples
#'
#' if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   BiocManager::install("GenomicRanges")
#' }
#'
#' library(GenomicRanges)
#'
#' seq <- data.frame(chr = c("seq1", "seq2"),
#'                   seq = c("ATTTGGGAGGGAGGGAGGGATGATGAAAATTTTATTTATTTTATTTA",
#'                           "TTTATACTATTCCCTTACCCTCCCATCCCCATACGGCATCTAGATC"))
#' seq_gr <- loadSequence(seq)
#' G4 <- G4HunterDetect(seq_gr)
#'
#' snv_gr <- GRanges(seqnames = c("seq1", "seq2"),
#'                   ranges = IRanges(start = c(18, 23), width = 1),
#'                   ref = c("G", "C"),
#'                   alt = c("C", "G"))
#'
#' effect <- SNVImpactG4(G4, snv_gr, alt_col = "alt")
#' plotImpactSeq(effect, ncol = 2)

plotImpactSeq <- function(filtered_gr,
                          ncol = 1) {

  if (!inherits(filtered_gr, "GRanges")) {
    stop("'filtered_gr' must be a GRanges object")
  }

  if (length(filtered_gr) == 0) {
    stop("Input GRanges object 'filtered_gr' is empty.")
  }

  if (!is.numeric(ncol) || ncol <= 0 || floor(ncol) != ncol) {
    stop("'ncol' must be a positive integer.")
  }

  fig_list <- lapply(seq_along(filtered_gr), function(i) {

    original_score <- filtered_gr$G4.info.score[i]
    mutated_score <- filtered_gr$mut.score[i]

    label <- paste0("Original G4Hunter score: ", original_score,
                    " \nMutated G4Hunter score: ", mutated_score)

    original_seq <- filtered_gr$G4.info.sequence[i]
    mutated_seq <- filtered_gr$mut.G4.seq[i]

    original_seq_char <- strsplit(original_seq, "")[[1]]
    mutated_seq_char <- strsplit(mutated_seq, "")[[1]]

    diff_positions <- which(original_seq_char != mutated_seq_char)

    aln_ori_mat <- data.frame(char = original_seq_char,
                              group = rep(paste0("Original\n(Score:",
                                                 original_score, ")"),
                                          length(original_seq_char)),
                              position = seq_along(original_seq_char),
                              mut = rep('no', length(original_seq_char)))

    aln_mut_mat <- data.frame(char = mutated_seq_char,
                              group = rep(paste0("Mutated\n(Score:",
                                                 mutated_score, ")"),
                                          length(mutated_seq_char)),
                              position = seq_along(mutated_seq_char),
                              mut = rep('no', length(mutated_seq_char)))

    aln_mut_mat$mut[diff_positions] <- "yes"

    aln_mat <- rbind(aln_ori_mat, aln_mut_mat)

    p1 <- ggseqlogo(setNames(list(c(original_seq, mutated_seq)), label),
                    method = 'prob') +
       annotate('rect', xmin = diff_positions - 0.5,
                xmax = diff_positions + 0.5, ymin = -0.05, ymax = 1.05,
                alpha = .2, col='black', fill='yellow') +
       geom_segment(aes(x = diff_positions, xend = diff_positions, y = -0.2,
                        yend = -0.1),
                    color = '#b22d2d', size = 0.8,
                    arrow = arrow(type = "closed",
                                  length = unit(0.2, "inches"))) +
       theme(axis.text.x = element_blank())

    p2 <- ggplot(aln_mat, aes(.data$position, .data$group)) +
      geom_text(aes(label = .data$char, color = .data$mut, size = .data$mut)) +
      scale_x_continuous(breaks = seq_len(10), expand = c(0.105, 0)) +
      xlab('') + ylab('') +
      scale_color_manual(values = c('black', '#b22d2d')) +
      scale_size_manual(values = c(5, 6)) +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())

    return(plot_grid(p1, p2, ncol = 1, align = 'v'))
  })

  combined_plot <- plot_grid(plotlist = fig_list, labels = "AUTO",
                             ncol = ncol, align = 'v')

  return(combined_plot)
}

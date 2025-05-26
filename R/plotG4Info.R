#' Plot Basic Statistics of G4s Detected by the G4HunterDetect Function
#'
#' This function generates a series of plots to visualize basic statistics of
#' G4 sequences predicted by the G4Hunter algorithm.
#' The function produces the following plots:
#' \itemize{
#'   \item Distribution of max scores (absolute values).
#'   \item Distribution of max scores, split by strand.
#'   \item Distribution of scores (absolute values).
#'   \item Distribution of scores, split by strand.
#'   \item Distribution of sequence lengths (with lengths greater than 50 bp
#'   grouped into a single bin).
#'   \item Distribution of sequence lengths, split by strand (with lengths
#'   greater than 50 bp grouped into a single bin).
#' }
#'
#' @param G4 A GRanges object containing G4 sequences. The object must include
#' at least the following metadata columns:
#' \itemize{
#'   \item \code{score}: Numeric vector of scores for G4 sequences.
#'   \item \code{max_score}: Numeric vector of maximum scores for G4 sequences.
#'   \item \code{sequence}: Character vector of G4 sequence strings.
#' }
#' The GRanges object must include \code{strand} information.
#' It should come from the results returned by the
#' \code{G4HunterDetect} function.
#'
#' @param p_colors A vector of colors to be used for the plots. It should
#' contain nine color values corresponding to the different plot elements.
#'
#' @return A combined plot displays all the generated plots arranged in a
#' grid layout.
#'
#' @seealso \code{\link{G4HunterDetect}} for detecting the G4 sequences in a
#' given \code{DNAStringSet} object.
#'
#' @import ggplot2
#' @import cowplot
#' @import GenomicRanges
#' @export
#'
#' @examples
#'
#' seq_path <- system.file("extdata", "seq.txt", package = "G4SNVHunter")
#' G4 <- G4HunterDetect(loadSequence(seq_path = seq_path))
#' plotG4Info(G4)
plotG4Info <- function(G4,
                       p_colors = c("#6B9ECC", "#D91117", "#0E619C",
                                    "#58AC7B", "#D91117", "#0E619C",
                                    "#F39C40", "#D91117", "#0E619C")) {
  if (length(p_colors) != 9) {
    stop("'p_colors' must be a vector of exactly nine colors.")
  }

  if (!inherits(G4, "GRanges")) {
    stop("The input object G4 must be a GRanges object.")
  }

  df <- as.data.frame(G4)

  p1 <- ggplot(df, aes(x = abs(.data$max_score))) +
    geom_histogram(binwidth = 0.1, color = "black", fill = p_colors[1]) +
    labs(title = "Distribution of max_score (Absolute)",
         x = "max_score", y = "Frequency") +
    theme_minimal()

  p2 <- ggplot(df, aes(x = .data$max_score, fill = .data$strand)) +
    geom_histogram(binwidth = 0.1, position = "identity",
                   alpha = 0.6, color = "black") +
    labs(title = "Distribution of max_score (By Strand)",
         x = "max_score", y = "Frequency") +
    theme_minimal() +
    scale_fill_manual(values = c("+" = p_colors[2], "-" = p_colors[3]))

  p3 <- ggplot(df, aes(x = abs(.data$score))) +
    geom_histogram(binwidth = 0.1, color = "black", fill = p_colors[4]) +
    labs(title = "Distribution of score (Absolute)",
         x = "score", y = "Frequency") +
    theme_minimal()

  p4 <- ggplot(df, aes(x = .data$score, fill = .data$strand)) +
    geom_histogram(binwidth = 0.1, position = "identity",
                   alpha = 0.6, color = "black") +
    labs(title = "Distribution of score (By Strand)",
         x = "score", y = "Frequency") +
    theme_minimal() +
    scale_fill_manual(values = c("+" = p_colors[5], "-" = p_colors[6]))

  df$length_group <- cut(df$width,
    breaks = c(seq(0, 50, by = 5), Inf),
    labels = c(paste0(seq(0, 45, by = 5), "-", seq(5, 50, by = 5)), ">50"),
    right = FALSE
  )

  p5 <- ggplot(df, aes(x = .data$length_group)) +
    geom_bar(color = "black", fill = p_colors[7]) +
    labs(title = "Distribution of Sequence Length",
         x = "Length", y = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  df$strand_length_group <- interaction(df$strand,
                                        df$length_group,
                                        drop = TRUE)

  p6 <- ggplot(df, aes(x = .data$strand_length_group, fill = .data$strand)) +
    geom_bar(position = "dodge", color = "black", alpha = 0.6) +
    labs(title = "Distribution of Sequence Length (By Strand)",
         x = "Length", y = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("+" = p_colors[8], "-" = p_colors[9]))

  # Combine all plots

  combined_plot <- plot_grid(p1, p2, p3, p4, p5, p6,
                             ncol = 2, labels = "AUTO", align = 'v')

  title <- ggdraw() +
    draw_label(
      "Basic Statistics of G4s Predicted by G4Hunter Algorithm",
      fontface = 'bold',
      x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(title, combined_plot,
                          ncol = 1, rel_heights = c(0.1, 1))

  return(final_plot)
}

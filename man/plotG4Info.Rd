% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotG4Info.R
\name{plotG4Info}
\alias{plotG4Info}
\title{Plot Basic Statistics of G4s Detected by the G4HunterDetect Function}
\usage{
plotG4Info(
  G4,
  p_colors = c("#6B9ECC", "#D91117", "#0E619C", "#58AC7B", "#D91117", "#0E619C",
    "#F39C40", "#D91117", "#0E619C")
)
}
\arguments{
\item{G4}{A GRanges object containing G4 sequences. The object must include
at least the following metadata columns:
\itemize{
  \item \code{score}: Numeric vector of scores for G4 sequences.
  \item \code{max_score}: Numeric vector of maximum scores for G4 sequences.
  \item \code{sequence}: Character vector of G4 sequence strings.
}
The GRanges object must include \code{strand} information.
It should come from the results returned by the
\code{G4HunterDetect} function.}

\item{p_colors}{A vector of colors to be used for the plots. It should
contain nine color values corresponding to the different plot elements.}
}
\value{
A combined plot displays all the generated plots arranged in a
grid layout.
}
\description{
This function generates a series of plots to visualize basic statistics of
G4 sequences predicted by the G4Hunter algorithm.
The function produces the following plots:
\itemize{
  \item Distribution of max scores (absolute values).
  \item Distribution of max scores, split by strand.
  \item Distribution of scores (absolute values).
  \item Distribution of scores, split by strand.
  \item Distribution of sequence lengths (with lengths greater than 50 bp
  grouped into a single bin).
  \item Distribution of sequence lengths, split by strand (with lengths
  greater than 50 bp grouped into a single bin).
}
}
\examples{

seq_path <- system.file("extdata", "seq.txt", package = "G4SNVHunter")
G4 <- G4HunterDetect(loadSequence(seq_path = seq_path))
plotG4Info(G4)
}
\seealso{
\code{\link{G4HunterDetect}} for detecting the G4 sequences in a
given \code{DNAStringSet} object.
}

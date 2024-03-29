% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAMAlignR_functions.R
\name{AdjustStartingPosition}
\alias{AdjustStartingPosition}
\title{Adjust the Starting Position of the read based on CIGAR String}
\usage{
AdjustStartingPosition(read_pos, read_CIGAR)
}
\arguments{
\item{read_pos}{Integer representing the original leftmost mapping position of the first matching base for the read of interest}

\item{read_CIGAR}{CIGAR string associated to the read of interest}
}
\value{
A new value (integer) for the adjusted starting position of the read mapping on the reference genome
}
\description{
This function changes the value POS of the read of interest,
according to the presence of soft-clipped bases (S in the CIGAR string) and insertions (I) at the beginning of the read.
The "S" operation in the CIGAR string indicates that the corresponding bases are present in the read but not included in the alignment to the reference.
Soft clipping is often used to represent regions of the read that may contain sequencing artifacts, adapter contamination, or regions with low-quality bases.
The aim of the function is to define the new starting position for extracting the region of interest from the reference genome.
}
\examples{
read_pos <- 100
read_CIGAR <- "3S2M"
adjusted_pos <- AdjustStartingPosition(read_pos, read_CIGAR)
cat("Adjusted Starting Position:", adjusted_pos, "\n")

}
\references{
\url{https://samtools.github.io/hts-specs/SAMv1.pdf}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAMAlignR_functions.R
\name{CigarAnalyzerRead}
\alias{CigarAnalyzerRead}
\title{Analyze Read Alignment with CIGAR string}
\usage{
CigarAnalyzerRead(read_CIGAR, read_seq_char, read_flag)
}
\arguments{
\item{read_CIGAR}{CIGAR string representing the alignment operations for the read}

\item{read_seq_char}{String representing the read sequence}

\item{read_flag}{FLAG value associated to the read of interest}
}
\value{
A list containing a string representing the aligned sequence of the read and its length, and a message indicating whether the read is mapping on the forward or reverse strand
}
\description{
This function analyses the alignment of a read using its CIGAR string and associated sequence. It generates the aligned sequence of the read, taking into account various
CIGAR operations (M, I, D, =, S, X, N).In particular, it adds a '-' into the sequence if the CIGAR operations are D (deletion from the reference),
or N (skipped region from the reference), in such a way that the read correctly aligns with the reference sequence.
Furthermore, it takes as input the read FLAG value to use the \link{CheckStrand} function to print out if the read is mapping on the forward strand or reverse strand.
The function returns information about the aligned sequence and its length.
}
\examples{
read_CIGAR <- "2S5M10D2M"
read_seq_char <- "ACGACTGCA"
read_flag <- 0
aligned_seq_info <- CigarAnalyzerRead(read_CIGAR, read_seq_char, read_flag)
cat("Aligned Sequence:", aligned_seq_info$aligned_sequence, "\n")
cat("Length of Aligned Sequence:", aligned_seq_info$length_aligned_seq, "\n")

}
\references{
\url{https://samtools.github.io/hts-specs/SAMv1.pdf}
}
\seealso{
\code{\link{CheckStrand}}
}

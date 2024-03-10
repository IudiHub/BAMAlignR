#' Extract information from a specific read and perform some checks
#'
#' @description
#' This function takes as input a BAM file path and a specific read ID and it extracts some relevant information for that read,
#' i.e. the sequence (chromosome or contig) within the reference genome where the read is aligned (RNAME), the leftmost mapping position of the first matching base within the read (POS),
#' the CIGAR string, the actual read sequence, the FLAG value, the strand to which the read maps.
#' Additionally, it verifies that the input parameters are correct and that the BAM file contains the proper information about the read of interest.
#'
#' @references \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#'
#' @param bamPath Path to the BAM file
#' @param readID A character string representing the ID of the read to extract information for
#'
#' @return A list containing the following information about the read of interest:
#' \item{read_rname}{The reference sequence name to which the read is aligned.}
#' \item{read_pos}{The leftmost mapping position of the first matching base within the read.}
#' \item{read_CIGAR}{The CIGAR string representing the alignment pattern of the read.}
#' \item{read_seq_char}{The sequence of the read, stored as a character string.}
#' \item{read_flag}{The flag value associated with the read.}
#'
#' #' @examples
#' bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package = "BAMAlignR")
#' readID <- "SRR062635.8528185"
#' readInfo <- ExtractReadInfoAndCheck(bamPath, readID)
#'
#' @import Rsamtools
#' @import GenomicRanges
#'
#' @export
ExtractReadInfoAndCheck <- function(bamPath, readID) {
  #bamPath <- system.file("extdata", "bamName.bam", package="ReadAligner")

  if (!file.exists(bamPath)) stop("The BAM file is not present in the extdata folder")

  # scan the BAM file to retrieve all the reads and their information
  reads <- scanBam(bamPath)[[1]]

  # get the index of a specific read in 'reads' list
  read_index <- which(reads$qname == readID)

  if (!is.character(readID)) stop("ReadID must be a character!")

  if (!(readID %in% reads$qname)) stop("That read ID is not present in the QNAME field of the BAM file. Try with another read ID.")

  if (readID == "*") stop("The information about the read ID is unavailable")

  # extract a specific read in the list based on its index
  read_of_interest <- lapply(reads, function(xx) xx[read_index])

  # if read ID of interest is associated with more than one value in the different fields of the list, take just the first one for simplicity
  if (length(read_of_interest[[1]]) > 1) {
    read_of_interest <- lapply(read_of_interest, function(x) x[1])
  }

  # extract the reference sequence name to which the read is aligned
  read_rname <- read_of_interest$rname # read_rname is a factor
  #read_rname_char <- as.character(read_rname)

  if (read_rname == "*") stop("The read is unmapped or not aligned to any specific sequence in the reference genome")

  # extract the leftmost mapping position of the first matching base within the read
  read_pos <- read_of_interest$pos

  if (read_pos == 0) stop("The read is unmapped without coordinate")

  # extract the CIGAR string for the read
  read_CIGAR <- read_of_interest$cigar

  if (read_CIGAR == "*") stop("The CIGAR string is unavailable")

  # extract the sequence of the read
  read_seq_width <- read_of_interest$seq  #read_seq_width is a DNAStringSet object of length 1 containing the width of the sequence and the sequence itself
  # convert the sequence into a character
  read_seq_char <- toString(read_seq_width[[1]]) #or as.character(read_seq_width[[1]])

  if (read_seq_char == "*") stop("The sequence of the read is not stored")

  # extract the flag value of the read
  read_flag <- read_of_interest$flag

  read_info <- list(read_rname = read_rname, read_pos = read_pos, read_CIGAR = read_CIGAR , read_seq_char = read_seq_char, read_flag = read_flag)
  return(read_info)
}
#' Check the strand orientation of the mapping read from the FLAG field of the BAM file
#'
#' @description
#' This function checks if the read is mapping on the forward strand or on the reverse strand by taking as input the FLAG value from the BAM file.
#' The 5th bit (16) of the FLAG field indicates whether the read is reverse complemented or not.
#' This bit is set (1) if the read is reverse complemented, and it is not set (0) if it is not reverse complemented.
#' It performs a bitwise AND operation ('&') between the flag field and the value 16 (00010000).
#' The strand orientation of the mapping read is printed out.
#'
#' @param read_flag FLAG value associated to the read of interest
#'
#' @return A message indicating whether the read is mapped to the forward or reverse strand.
#'
#'
#' @examples
#' CheckStrand(16) # Output: "The read is mapping to the reverse strand"
#' CheckStrand(0)  # Output: "The read is mapping to the forward strand"
#'
#' @export
CheckStrand <- function(read_flag) {
  # check if the 5th value is set
  is_reverse_strand <- (bitwAnd(read_flag, 16) != 0)
  if (is_reverse_strand) {
    print("The read is mapping to the reverse strand")
  } else {
    print("The read is mapping to the forward strand")
  }
}
#' Analyze Read Alignment with CIGAR string
#'
#' @description
#' This function analyses the alignment of a read using its CIGAR string and associated sequence. It generates the aligned sequence of the read, taking into account various
#' CIGAR operations (M, I, D, =, S, X, N).In particular, it adds a '-' into the sequence if the CIGAR operations are D (deletion from the reference),
#' or N (skipped region from the reference), in such a way that the read correctly aligns with the reference sequence.
#' Furthermore, it takes as input the read FLAG value to use the \link{CheckStrand} function to print out if the read is mapping on the forward strand or reverse strand.
#' The function returns information about the aligned sequence and its length.
#'
#' @references \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#'
#' @param read_CIGAR CIGAR string representing the alignment operations for the read
#' @param read_seq_char String representing the read sequence
#' @param read_flag FLAG value associated to the read of interest
#'
#' @return A list containing a string representing the aligned sequence of the read and its length, and a message indicating whether the read is mapping on the forward or reverse strand
#'
#' @seealso \code{\link{CheckStrand}}
#'
#' @examples
#' read_CIGAR <- "2S5M10D2M"
#' read_seq_char <- "ACGACTGCA"
#' read_flag <- 0
#' aligned_seq_info <- CigarAnalyzerRead(read_CIGAR, read_seq_char, read_flag)
#' cat("Aligned Sequence:", aligned_seq_info$aligned_sequence, "\n")
#' cat("Length of Aligned Sequence:", aligned_seq_info$length_aligned_seq, "\n")
#'
#' @import stringr
#'
#' @export
CigarAnalyzerRead <- function(read_CIGAR, read_seq_char, read_flag) {

  # Extract the integers from the CIGAR string
  integers <- as.numeric(str_extract_all(read_CIGAR, "\\d+")[[1]])

  # Extract the characters from the CIGAR string
  characters <- str_replace_all(read_CIGAR, "\\d+", "")
  cigar_operations <- strsplit(characters, "")[[1]]

  # Initialize the aligned sequence of the read
  aligned_sequence <- ""

  # Iterate over the integers and the cigar operations (i.e., the characters of the CIGAR string) simultaneously
  for (i in seq_along(integers)) {
    op_type <- cigar_operations[i]

    # if the CIGAR operations are M, I, =, S, X, rewrite the read sequence as it is
    if (op_type %in% c("M", "I", "=", "S", "X")) {
      aligned_sequence <- paste0(aligned_sequence, substring(read_seq_char, 1, integers[i]))
      read_seq_char <- substring(read_seq_char, integers[i] + 1, nchar(read_seq_char))
      # else if the CIGAR operations are D, N, add a "-"
    } else if (op_type %in% c("D", "N")) {
      aligned_sequence <- paste0(aligned_sequence, strrep("-", integers[i]))
    }
  }
  CheckStrand(read_flag)
  # length of the aligned sequence
  length_aligned_seq <- nchar(aligned_sequence)
  aligned_seq_info <- list(aligned_sequence = aligned_sequence, length_aligned_seq = length_aligned_seq)
  return(aligned_seq_info)
}

#' Adjust the Starting Position of the read based on CIGAR String
#'
#' @description
#' This function changes the value POS of the read of interest,
#' according to the presence of soft-clipped bases (S in the CIGAR string) and insertions (I) at the beginning of the read.
#' The "S" operation in the CIGAR string indicates that the corresponding bases are present in the read but not included in the alignment to the reference.
#' Soft clipping is often used to represent regions of the read that may contain sequencing artifacts, adapter contamination, or regions with low-quality bases.
#' The aim of the function is to define the new starting position for extracting the region of interest from the reference genome.
#'
#' @references \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#'
#' @param read_pos Integer representing the original leftmost mapping position of the first matching base for the read of interest
#' @param read_CIGAR CIGAR string associated to the read of interest
#'
#' @return A new value (integer) for the adjusted starting position of the read mapping on the reference genome
#'
#' @examples
#' read_pos <- 100
#' read_CIGAR <- "3S2M"
#' adjusted_pos <- AdjustStartingPosition(read_pos, read_CIGAR)
#' cat("Adjusted Starting Position:", adjusted_pos, "\n")
#'
#' @import stringr
#'
#' @export
AdjustStartingPosition <- function(read_pos, read_CIGAR) {

  # Extract the integers from the CIGAR string
  integers <- as.numeric(str_extract_all(read_CIGAR, "\\d+")[[1]])

  # Extract the characters from the CIGAR string
  characters <- str_replace_all(read_CIGAR, "\\d+", "")
  cigar_operations <- strsplit(characters, "")[[1]]

  for (i in seq_along(integers)) {
    op_type <- cigar_operations[i]
    integers_0 <- integers[1]

    if (op_type %in% c("M", "D", "N", "=", "X")) {
      break
    } else {
      # if the first operation of the CIGAR string is S (or I), update the starting position
      # by subtracting from the original position the associated number of operations
      read_pos <- read_pos - integers[i]
    }
  }
  return(read_pos)
}

#' Extract the reference sequence
#'
#' @description
#' This function takes as input the sequence within the reference genome where the read is aligned (i.e. RNAME in the BAM file, that is extracted in the function \link{ExtractReadInfoAndCheck}),
#' the length of the reference sequence you need to extract (i.e., the length of the aligned read, that is computed in the function \link{CigarAnalyzerRead}),
#' the position of the nucleotide in the reference genome from which the extracted sequence starts (i.e. POS in the BAM file, that is the adjusted with the function
#' \link{AdjustStartingPosition} and the CIGAR string.
#' From this information, the function extracts the reference sequence to which the read is aligned.
#'
#' @param read_rname String representing the reference sequence name to which the read is aligned
#' @param length_aligned_seq Integer representing the length of the aligned read
#' @param read_pos Integer representing the starting position of the mapping read on the reference sequence
#' @param read_CIGAR CIGAR string associated to the read of interest
#'
#' @return A string containing the extracted reference sequence
#'
#' @seealso \code{\link{AdjustStartingPosition}}
#'
#' @examples
#' read_rname <- 11
#' length_aligned_seq <- 100
#' read_pos <- 100
#' read_CIGAR <- "100M"
#' ref_sequence <- ExtractReferenceSequence(read_rname, length_aligned_seq, read_pos, read_CIGAR)
#'
#' cat("Extracted Reference Sequence:", ref_sequence, "\n")
#'
#' @import GenomicRanges
#' @import IRanges
#' @import BSgenome.Hsapiens.NCBI.GRCh38
#'
#' @export
ExtractReferenceSequence <- function(read_rname, length_aligned_seq, read_pos, read_CIGAR) {
  read_pos <- AdjustStartingPosition(read_pos, read_CIGAR)
  read_rname <- as.character(read_rname)
  # define the genomic range of the reference sequence of interest
  gr <- GRanges(read_rname, IRanges(read_pos, width=length_aligned_seq))
  ref_seq <- getSeq(Hsapiens, gr)          # DNAStringSet object of length 1 containing the sequence and its width
  ref_string <- toString(ref_seq[[1]])
  return(ref_string)
}
#' Analyze the CIGAR string and print the aligned reference sequence
#'
#' @description
#' This function uses the \link{ExtractReferenceSequence} function to extract the reference sequence of interest, and then based on the CIGAR string, it adds a '-' into the sequence
#' if the CIGAR operation is I (i.e., insertion to the reference sequence), in such a way that the reference sequence correctly aligns with the read.
#' It uses the same logic as \link{CigarAnalyzerRead} function, but inverting insertions and deletions/skipped regions.
#'
#' @references \url{https://samtools.github.io/hts-specs/SAMv1.pdf}
#'
#' @param read_CIGAR CIGAR string associated to the read of interest
#' @param read_rname String representing the reference sequence name to which the read is aligned
#' @param length_aligned_seq Integer representing the length of the aligned read
#' @param read_pos Integer representing the starting position of the mapping read on the reference sequence
#'
#' @return A string representing the reference sequence that aligns with the read
#'
#' @seealso \code{\link{ExtractReferenceSequence}}
#'
#' @examples
#' read_CIGAR <- "3M2I5M"
#' read_seq_char <- "TCTAACGC"
#' read_rname <- 11
#' length_aligned_seq <- 10
#' read_pos <- 100
#' aligned_ref_seq <- CigarAnalyzerReference(read_CIGAR, read_rname, length_aligned_seq, read_pos)
#' cat("Aligned Reference Sequence:", aligned_ref_seq, "\n")
#'
#' @import stringr
#'
#' @export
CigarAnalyzerReference <- function(read_CIGAR, read_rname, length_aligned_seq, read_pos) {
  ref_string <- ExtractReferenceSequence(read_rname, length_aligned_seq, read_pos, read_CIGAR)
  # Extract integers from the CIGAR string
  integers <- as.numeric(str_extract_all(read_CIGAR, "\\d+")[[1]])

  # Extract characters from the CIGAR string
  characters <- str_replace_all(read_CIGAR, "\\d+", "")
  cigar_operations <- strsplit(characters, "")[[1]]

  # Initialize aligned_sequence
  aligned_ref_seq <- ""

  # Iterate over the integers and cigar_operations simultaneously
  for (i in seq_along(integers)) {
    op_type <- cigar_operations[i]

    # if the CIGAR operations are M, =, S, X, D, N, rewrite the reference sequence as it is
    if (op_type %in% c("M", "=", "S", "X", "D", "N")) {
      aligned_ref_seq <- paste0(aligned_ref_seq, substring(ref_string, 1, integers[i]))
      ref_string <- substring(ref_string, integers[i] + 1, nchar(ref_string))
      # else if the CIGAR operations is I add a "-"
    } else if (op_type %in% c("I")) {
      aligned_ref_seq <- paste0(aligned_ref_seq, strrep("-", integers[i]))
    }
  }
  return(aligned_ref_seq)
}
#' Show the read alignment and the CIGAR string
#'
#' @description
#' Core function of the package: it shows the precise alignment of the read to the reference genome according to the mapped
#' read's CIGAR string. It first extracts relevant information using the \link{ExtractReadInfoAndCheck} function. It then
#' utilizes the \link{CigarAnalyzerRead} and \link{CigarAnalyzerReference} functions to process the CIGAR string and generate
#' aligned sequences for both the read and the reference. The function then prints the CIGAR string, the aligned reference
#' sequence, and the aligned read sequence.
#'
#'
#' @param bamPath Path to the BAM file
#' @param readID A character string representing the ID of the read
#'
#'
#' @return A representation of the  precise alignment of the read to the reference genome according to the mapped read's CIGAR string
#'
#' @seealso \code{\link{ExtractReadInfoAndCheck}}, \code{\link{CigarAnalyzerRead}}, \code{\link{CigarAnalyzerReference}}
#'
#' @examples
#' bamPath <- system.file(
#' "extdata",
#' "HG00096.chrom11.ILLUMINA.bwa.bam",
#' package = "BAMAlignR"
#' )
#' readID <- "SRR062635.8528185"
#' ShowAlignment(bamPath, readID)
#'
#' @export
ShowAlignment <- function(bamPath, readID) {
  read_info <- ExtractReadInfoAndCheck(bamPath, readID)
  read_flag <- read_info$read_flag
  read_seq_char <- read_info$read_seq_char
  read_CIGAR <- read_info$read_CIGAR
  read_pos <- read_info$read_pos
  read_rname <- read_info$read_rname
  aligned_seq_info <- CigarAnalyzerRead(read_CIGAR, read_seq_char, read_flag)
  aligned_sequence <- aligned_seq_info$aligned_sequence
  length_aligned_seq <- aligned_seq_info$length_aligned_seq
  aligned_ref_seq <- CigarAnalyzerReference(read_CIGAR, read_rname, length_aligned_seq, read_pos)
  cat("CIGAR:", read_CIGAR, "\n", "Ref.:", aligned_ref_seq, "\n", "Read:", aligned_sequence, "\n")
}



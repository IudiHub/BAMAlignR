## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BAMAlignR)
library(stringr)
library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(IRanges)

## -----------------------------------------------------------------------------
# Specify the BAM file path and read ID
library(BAMAlignR)
bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package = "BAMAlignR")
readID <- "SRR062635.8528185"

# Display alignment information using the core function
ShowAlignment(bamPath, readID)

## -----------------------------------------------------------------------------
ExtractReadInfoAndCheck <- function(bamPath, readID) {

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

bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package = "BAMAlignR")
readID <- "SRR062635.8528185"
print(ExtractReadInfoAndCheck(bamPath, readID))

## -----------------------------------------------------------------------------
CheckStrand <- function(read_flag) {
  # check if the 5th value is set
  is_reverse_strand <- (bitwAnd(read_flag, 16) != 0)
  if (is_reverse_strand) {
    print("The read is mapping to the reverse strand")
  } else {
    print("The read is mapping to the forward strand")
  }
}

CheckStrand(16)
CheckStrand(0)

## -----------------------------------------------------------------------------
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

read_CIGAR <- "2S5M10D2M"
read_seq_char <- "ACGACTGCA"
read_flag <- 0
aligned_seq_info <- CigarAnalyzerRead(read_CIGAR, read_seq_char, read_flag)
cat("Aligned Sequence:", aligned_seq_info$aligned_sequence, "\n")
cat("Length of Aligned Sequence:", aligned_seq_info$length_aligned_seq, "\n")

## -----------------------------------------------------------------------------
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

read_pos <- 100
read_CIGAR <- "3S2M"
adjusted_pos <- AdjustStartingPosition(read_pos, read_CIGAR)
cat("Adjusted Starting Position:", adjusted_pos, "\n")

## -----------------------------------------------------------------------------
ExtractReferenceSequence <- function(read_rname, length_aligned_seq, read_pos, read_CIGAR) {
  read_pos <- AdjustStartingPosition(read_pos, read_CIGAR)
  read_rname <- as.character(read_rname)
  # define the genomic range of the reference sequence of interest
  gr <- GRanges(read_rname, IRanges(read_pos, width=length_aligned_seq))
  ref_seq <- getSeq(Hsapiens, gr)          # DNAStringSet object of length 1 containing the sequence and its width
  ref_string <- toString(ref_seq[[1]])
  return(ref_string)
}

read_rname <- 11
length_aligned_seq <- 100
read_pos <- 100
read_CIGAR <- "100M"
ref_sequence <- ExtractReferenceSequence(read_rname, length_aligned_seq, read_pos, read_CIGAR)

## -----------------------------------------------------------------------------
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

read_CIGAR <- "6M3I2M"
read_seq_char <- "TCTAACGC"
read_rname <- 11
length_aligned_seq <- 100
read_pos <- 100
aligned_ref_seq <- CigarAnalyzerReference(read_CIGAR, read_rname, length_aligned_seq, read_pos)
cat("Aligned Reference Sequence:", aligned_ref_seq, "\n")

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
sessionInfo()


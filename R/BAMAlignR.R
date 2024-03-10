#' BAMAlignR package
#'
#' This package is designed to read in a BAM/SAM file, allow users to choose a specific read entry
#' and display how that read's sequence aligns to the reference genome using its CIGAR string.
#'
#' @docType package
#'
#' @author Cristina Iudica \email{cristina.iudica@mail.polimi.it}
#'
#' @name BAMAlignR
#' @export ExtractReadInfoAndCheck
#' @export CheckStrand
#' @export CigarAnalyzerRead
#' @export AdjustStartingPosition
#' @export ExtractReferenceSequence
#' @export CigarAnalyzerReference
#' @export ShowAlignment
#' @import Rsamtools
#' @import GenomicRanges
#' @import stringr
#' @import IRanges
#' @import BSgenome.Hsapiens.NCBI.GRCh38
NULL

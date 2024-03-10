library(Rsamtools)
library(GenomicRanges)

test_that("ReadID must be a character", {
  bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package="BAMAlignR")
  readID <- 289382

  expect_error(ExtractReadInfoAndCheck(bamPath, readID), "ReadID must be a character!")
})


test_that("Read ID not present in BAM file", {
  bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package="BAMAlignR")
  readID <- 'non-existent read'
  expect_error(ExtractReadInfoAndCheck(bamPath, readID), "That read ID is not present in the QNAME field of the BAM file.")
})


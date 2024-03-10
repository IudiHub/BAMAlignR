test_that("No error occurs during execution of the core ShowAlignment function", {
  expect_no_error({
    bamPath <- system.file("extdata", "HG00096.chrom11.ILLUMINA.bwa.bam", package="BAMAlignR")
    readID <- "SRR062635.8528185"
    ShowAlignment(bamPath, readID)
  })
})

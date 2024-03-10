test_that("The output of ExtractReferenceSequence is a character", {
  read_rname <- 11
  length_aligned_seq <- 9
  read_pos <- 100
  read_cigar <- "3M2D4M"

  ref_string <- ExtractReferenceSequence(read_rname, length_aligned_seq, read_pos, read_cigar)

  # Check if the output is a DNAString object of the expected length
  expect_is(ref_string, "character")
})

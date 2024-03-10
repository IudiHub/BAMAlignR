test_that("The length of the aligned read sequence is equal to the length of the aligned reference sequence", {
  read_cigar <- "5M2D4M"
  read_rname <- 11
  read_length <- 11
  read_pos <- 100

  aligned_ref_seq <- CigarAnalyzerReference(read_cigar, read_rname, read_length, read_pos)

  # Check if the length of the output matches the expected length
  expect_equal(nchar(aligned_ref_seq), read_length)
})

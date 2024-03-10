test_that("Aligned sequence length matches expected", {
  read_cigar <- "3M2D4M"
  read_seq_char <- "ATCGTAC"
  read_flag <- 0

  aligned_seq_info <- CigarAnalyzerRead(read_cigar, read_seq_char, read_flag)
  expect_equal(aligned_seq_info$length_aligned_seq, 9)
})


test_that("Aligned sequence contains proper characters", {
  read_cigar <- "3M2D4M"
  read_seq_char <- "ATCGTAC"
  read_flag <- 0

  aligned_seq_info <- CigarAnalyzerRead(read_cigar, read_seq_char, read_flag)
  expect_equal(aligned_seq_info$aligned_sequence, "ATC--GTAC")
})


test_that("Aligned sequence contains proper characters", {
  read_cigar <- "8M"
  read_seq_char <- "ATCGTACG"
  read_flag <- 0

  aligned_seq_info <- CigarAnalyzerRead(read_cigar, read_seq_char, read_flag)
  expect_equal(aligned_seq_info$aligned_sequence, "ATCGTACG")
})

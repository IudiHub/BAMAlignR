test_that("Starting position adjustment with soft-clipped operations", {
  read_pos <- 100
  read_cigar <- "5S3M2D4M"

  read_pos <- AdjustStartingPosition(read_pos, read_cigar)
  expect_equal(read_pos, 95)
})

test_that("Starting position adjustment without soft-clipped operations", {
  read_pos <- 100
  read_cigar <- "3M2D4M"

  read_pos <- AdjustStartingPosition(read_pos, read_cigar)
  expect_equal(read_pos, 100)
})

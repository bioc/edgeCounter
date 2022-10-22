test_that("bamEdges() working", {
  test.bam <- test_path("fixtures", "sample.bam")
  test.edges <- test_path("fixtures", "sample.bamEdges.rds")
  expect_true(
    all(bamEdges(test.bam) == readRDS(test.edges))
  )
})

test_that("countEdges() working", {
  test.bam <- test_path("fixtures", "sample.bam")
  test.range <- readRDS(test_path("fixtures", "sample.testGRanges.rds"))
  expect_equal(
    countEdges(test.range, test.bam)@elementMetadata@listData$edge.counts,
    c(2, 4, 0)
  )
})

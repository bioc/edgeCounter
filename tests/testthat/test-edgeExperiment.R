test_that("edgeExtract() works", {
    bam.path <- test_path("fixtures", "sample.bam")
    edges <- edgeExtract(c(a = bam.path, b = bam.path))
    expect.edges <- readRDS(test_path("fixtures", "sample.bamEdges.rds"))

    expect_true(
        all(names(edges) == c("a", "b"))
    )
    expect_true(
        all(
            # edges == expect.edges gives a LogicalList of length 2;
            # each list element is a logical vector
            #   ... indicating whether individual ranges are the same
            sapply(edges == expect.edges, function(x) all(x))
        )
    )
})

test_that("edgeExperimentFromCounts() works with GRanges and GRangesList", {
    # Get edges by edgeExtract
    bam.path <- test_path("fixtures", "sample.bam")
    edges <- edgeExtract(c(a = bam.path, b = bam.path))

    # Get a GRanges from an example in test data
    range.path <- test_path("fixtures", "sample.testGRanges.rds")
    range <- readRDS(range.path)

    # First, test with GRanges
    ranges <- range
    res1 <- edgeExperimentFromCounts(edges, ranges)
    #   Count data should be expected
    res1.expected <- matrix(c(2, 4, 0, 2, 4, 0), ncol = 2)
    colnames(res1.expected) <- c("a", "b")
    expect_equal(
        res1@assays@data@listData$counts,
        res1.expected
    )
    #   rowRanges should be expected
    expect_identical(SummarizedExperiment::rowRanges(res1), breakGRanges(ranges))


    # Next, test with GRangesList
    ranges <- GenomicRanges::GRangesList(
        list(range, range),
        compress = TRUE
    )
    res2 <- edgeExperimentFromCounts(edges, ranges)
    #   Count data should be expected
    res2.expected <- matrix(c(6, 6, 6, 6), ncol = 2)
    colnames(res2.expected) <- c("a", "b")
    expect_equal(
        res2@assays@data@listData$counts,
        res2.expected
    )
    #   rowRanges should be expected
    expect_identical(SummarizedExperiment::rowRanges(res2), ranges)
})

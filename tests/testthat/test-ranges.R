test_that("consensusGRanges() works with minimum.width defined", {
    # Testing by an artifical GRanges.
    #   NOTE - if the artifical GRanges described here changes,
    #     We need to save the `.rds` file again!
    gr <- makeGRanges(
        tibble::tribble(
            ~start, ~end, ~seqnames,
            1, 100, "chr2L",
            80, 120, "chr2L",
            200, 300, "chr2L",
            280, 320, "chr2L",
            360, 460, "chr2L",
            360, 400, "chr2L",
            100, 300, "chrX"
        )
    )
    expect_true({
        expect <- readRDS(test_path("fixtures", "consensusGRanges.example.rds"))
        all(expect == consensusGRanges(gr, minimum.width = 100))
    })
})


test_that("consensusGRanges() works with minimum.width NULL mode", {
    # Fully coded input and output here; no saved RDS file
    gr <- makeGRanges(
        tibble::tribble(
            ~start, ~end, ~seqnames,
            100, 200, "chr2L",
            150, 240, "chr2L",
            300, 500, "chr2L",
            490, 510, "chr2L",
            100, 200, "chrX",
            300, 400, "chrX"
        )
    )
    expect_true({
        expect <- makeGRanges(
            tibble::tribble(
                ~start, ~end, ~seqnames,
                100, 240, "chr2L",
                300, 510, "chr2L",
                100, 200, "chrX",
                300, 400, "chrX"
            )
        )
        all(expect == consensusGRanges(gr, minimum.width = NULL))
    })
})


test_that("readNarrowPeaks() works for single input", {
    # Testing by a single narrowPeak file provided.
    #   NOTE - the `.narrowPeak` and the `.rds` are linked.
    #     If we use a another `.narrowPeak`, we need to manually save `.rds` again!
    example.path <- test_path(
        "fixtures",
        "readNarrowPeaks.example.narrowPeak"
    )
    example.res.path <- test_path(
        "fixtures",
        "readNarrowPeaks.example.rds"
    )
    grs <- readNarrowPeaks(example.path)

    expect_equal(length(grs), 1)
    expect_equal(names(grs), example.path)
    expect_true(
        all(grs[[1]] == readRDS(example.res.path))
    )
})

test_that("readNarrowPeaks() works for multiple inputs with names", {
    # Testing by a single narrowPeak file provided, reading in two times.
    #   NOTE - the `.narrowPeak` and the `.rds` are linked.
    #     If we use a another `.narrowPeak`, we need to manually save `.rds` again!
    example.path <- test_path(
        "fixtures",
        "readNarrowPeaks.example.narrowPeak"
    )
    example.res.path <- test_path(
        "fixtures",
        "readNarrowPeaks.example.rds"
    )
    example.names <- c("f1", "f2")
    grs <- readNarrowPeaks(rep(example.path, times = 2),
        names = example.names
    )

    expect_equal(length(grs), 2)
    expect_true(
        all(names(grs) == example.names)
    )
    expect_true(
        all(grs[[1]] == readRDS(example.res.path))
    )
    expect_true(
        all(grs[[2]] == readRDS(example.res.path))
    )
})

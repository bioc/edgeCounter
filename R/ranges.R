# How to establish the intended GRanges object for counting?
#   The following methods are suggested:
#     1. Read MACS2 narrowPeaks or other BED-like files
#     2. Genomic features from packages, e.g., `GenomicFeatures::TxDb`.
#   This package provides a set of helper functions to deal with #1.
#   To process #2 data, refer to `inter-range-methods` of `GenomicRanges`.

#' Read in MACS2 narrowPeaks from one/multiple files into GRanges objects
#'
#' @description `readNarrowPeaks()` reads in MACS2 narrowPeaks from files
#' into a list of GRanges object(s).
#'
#'
#' @param paths A character vector of MACS2 narrowPeak file paths
#' @param names (optional) A character vector for names of the GRanges
#' @returns A named list of GRanges objects. Names are the input paths.
#' If `names` is provided, `names` is used
#' as list names instead of input paths.
#'
#' @examples
#' np.path <- system.file("extdata", "readNarrowPeaks.example.narrowPeak",
#'     package = "edgeCounter"
#' )
#' readNarrowPeaks(np.path, names = "sample") # Returns a list of GRanges
#' @export
readNarrowPeaks <- function(paths, names = NULL) {
    # Extra columns of narrowPeak file
    extraCols <- list(
        names = c(
            "name", "score", "strand", "foldChange",
            "-log10Pvalue", "-log10Qvalue", "summitPos"
        ),
        types = "cicdddi"
    )

    # Read in each narrowPeak file using the `importBED` of this package
    grs <- lapply(paths, function(path) {
        importBED(path = path, extra.cols = extraCols)
    })

    # Assign list names
    if (is.null(names)) {
        names(grs) <- paths
    } else {
        stopifnot(length(names) == length(paths))
        names(grs) <- names
    }

    grs
}

#' Compute 'consensus' GRanges from a list of GRanges objects
#'
#' @description `consensusGRanges()` computes consensus GRanges from a list of
#' GRanges objects. Definition of 'consensus' is explained in Details section.
#'
#'
#' @param grs A single GRanges or a list of GRanges, could be length 1 or more.
#' @param minimum.width Minimum width of the consensus ranges. Could be either
#' positive integer or NULL (which will trigger two computation modes, see
#' detail).
#' @returns A GRanges object consisting the consensus ranges.
#'
#' @details
#' # Details
#'
#' This section describes behavior of this function in detail. Two modes are
#' supported by this function, which are defined by value of `minimum.width`.
#'
#' ## Consensus algorithm (when `minimum.width` is NOT `NULL`)
#'
#' The following steps are adopted to compute a single `GRanges` object from
#' a list of input `GRanges`.
#' Inter-range methods of the `GenomicRanges` package are used for computation.
#'
#' 1. Concatenate all ranges from all GRanges in list into a single GRanges.
#' 2. Get a collection of non-overlapping ranges by `disjoin`.
#' 3. Only retain non-overlapping ranges that overlaps with at least two ranges
#' from the original concatenated GRanges.
#' 4. Merge connecting (i.e., distance is zero bp) ranges by `reduce`.
#' 5. Extend ranges to >= `minimum.width`bp in length (or keep unchanged).
#' 6. Merge connecting ranges by `reduce`.
#' 7. Those ranges are considered to be 'consensus'.
#'
#' ### Minimum width of the consensus ranges
#'
#' When `minimum.width` is a defined integer, the consensus ranges computed by
#' the algorithm above is enforced to be at least of `minimum.width` in width.
#' If ranges are not wide enough, they are extended to `minimum.width` width by
#' `GenomicRanges::resize(..., fix = "center")`.
#'
#' ## Consensus algorithm (when `minimum.width` is `NULL`)
#'
#' When minimum.width is set to be `NULL`, the function tries to get "widest"
#' and "most" possible consensus ranges. The following steps are adopted:
#'
#' 1. Concatenate all ranges from all GRanges in list into a single GRanges.
#' 2. Get a collection of non-overlapping ranges by `disjoin`.
#' 3. Merge connecting (i.e., distance is zero bp) ranges by `reduce`.
#' 4. Those ranges are considered to be 'consensus'.
#'
#' Therefore, under this mode even ranges that occur only once will be reported
#' as a "consensus range".
#'
#' @seealso
#' The following inter-range methods from the `GenomicRanges` package are used:
#'
#' * [GenomicRanges::disjoin()]
#' * [GenomicRanges::reduce()]
#' * [GenomicRanges::countOverlaps()]
#'
#' For `countOverlaps`, this functions always adopt parameter `type = "any"`.
#'
#' @examples
#' grs <- makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         100, 200, "chr2L",
#'         150, 250, "chr2L",
#'         200, 400, "chrX",
#'         250, 500, "chrX"
#'     )
#' )
#' consensusGRanges(grs) # With default, a minimum width is enforced
#'
#' grs <- makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         100, 240, "chr2L",
#'         300, 510, "chr2L",
#'         100, 200, "chrX"
#'     )
#' )
#' consensusGRanges(grs, minimum.width = NULL)
#'
#' @export
consensusGRanges <- function(grs, minimum.width = 100) {
    # Sanity check of minimum.width
    stopifnot((is.numeric(minimum.width) & minimum.width > 0) |
        is.null(minimum.width))
    # Convert into a single GRanges of peaks and strip metadata
    if (is.list(grs)) grs <- Reduce(c, grs)
    stopifnot(inherits(grs, "GRanges"))
    GenomicRanges::mcols(grs) <- NULL

    # Get disjoin set - those are potential consensus peaks
    grs.dj <- GenomicRanges::disjoin(grs)

    if (!is.null(minimum.width)) {
        # Disjoin ranges that overlap >1 times with the original
        grs.dj.count <- GenomicRanges::countOverlaps(grs.dj, grs, type = "any")
        consensus <- grs.dj[grs.dj.count > 1]
        consensus <- GenomicRanges::reduce(consensus)

        # For those width < minimum.width adjust their width with center fixed
        consensus.small.set <- GenomicRanges::width(consensus) < minimum.width
        consensus.small <- consensus[consensus.small.set]
        consensus.small <- GenomicRanges::resize(
            x = consensus.small,
            width = minimum.width,
            fix = "center"
        )
        consensus[consensus.small.set] <- consensus.small

        # Reduce the ranges
        consensus <- GenomicRanges::reduce(consensus, drop.empty.ranges = TRUE)
    } else {
        consensus <- GenomicRanges::reduce(grs.dj, drop.empty.ranges = TRUE)
    }

    consensus
}

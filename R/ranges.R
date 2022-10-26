# How to establish the intended GRanges object for counting?
#   The following methods are suggested:
#     1. Read MACS2 narrowPeaks or other BED-like files
#     2. Genomic features from established packages, e.g., `GenomicFeatures::TxDb`.
#   This package provides a set of helper functions to deal with #1.
#   To further process #2 data, refer to `inter-range-methods` of `GenomicRanges`.

#' Read in MACS2 narrowPeaks from one/multiple files into GRanges objects
#'
#' @description `readNarrowPeaks()` reads in MACS2 narrowPeaks from files
#' into a list of GRanges object(s).
#'
#'
#' @param paths A character vector of MACS2 narrowPeak file paths
#' @param names (optional) A character vector for names of the GRanges
#' @returns A named list of GRanges objects. Names are the input paths.
#' If `names` is provided, `names` is used as list names instead of input paths.
#'
#' @export
readNarrowPeaks <- function(paths, names = NULL){

  # Extra columns of narrowPeak file
  extraCols <- list(
    names = c("name", "score", "strand", "foldChange",
              "-log10Pvalue", "-log10Qvalue", "summitPos"),
    types = "cicdddi"
  )

  # Read in each narrowPeak file using the `importBED` of this package
  grs <- lapply(paths, function(path){
    importBED(path = path, extra.cols = extraCols)
  })

  # Assign list names
  if (is.null(names)){
    names(grs) <- paths
  } else{
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
#' @param grs A list of GRanges objects, could be length 1 or more.
#' @param minimum.width Minimum width of the consensus ranges.
#' @returns A GRanges object consisting the consensus ranges.
#'
#' @details
#' # Details
#'
#' This section describes behavior of this function in detail.
#'
#' ## Consensus algorithm
#'
#' The following steps are adopted to compute a single `GRanges` object from
#' a list of input `GRanges`. Inter-range methods of the `GenomicRanges` package
#' are used for computation.
#'
#' 1. Concatenate all ranges from all GRanges in the list into a single GRanges.
#' 2. Get a collection of non-overlapping ranges by `disjoin`.
#' 3. Only retain non-overlapping ranges that overlaps with at least two ranges
#' from the original concatenated GRanges.
#' 4. Merge connecting (i.e., distance is zero bp) ranges by `reduce`.
#' 5. Extend ranges to at least `minimum.width`bp in length (or keep unchanged).
#' 6. Merge connecting ranges by `reduce`.
#' 7. Those ranges are considered to be 'consensus'.
#'
#' ## Minimum width of the consensus ranges
#'
#' TODO
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
#' @export
consensusGRanges <- function(grs, minimum.width = 100){

  # Convert into a single GRanges of peaks and strip metadata
  if (is.list(grs)) grs <- Reduce(c, grs)
  stopifnot(class(grs) == "GRanges")
  GenomicRanges::mcols(grs) <- NULL

  # Get disjoin set - those are potential consensus peaks
  grs.dj <- GenomicRanges::disjoin(grs)

  # Consensus range PART I - disjoin ranges that overlap >1 times with the original
  grs.dj.count <- GenomicRanges::countOverlaps(grs.dj, grs, type = "any")
  consensus.p1 <- grs.dj[grs.dj.count > 1]
  consensus.p1 <- GenomicRanges::reduce(consensus.p1)

  # TODO - DO WE WANT THIS?
  # Consensus range PART II - ranges that overlap with nobody else (original)
  #grs.self.count <- GenomicRanges::countOverlaps(grs, grs, type = "any")
  #consensus.p2 <- grs[grs.self.count == 1]

  # Build consensus
  #consensus <- c(consensus.p1, consensus.p2)
  consensus <- consensus.p1

  # For those width < minimum.width adjust their width with center fixed
  consensus.small <- consensus[GenomicRanges::width(consensus) < minimum.width]
  consensus.small <- GenomicRanges::resize(x = consensus.small,
                                           width = minimum.width,
                                           fix = "center")
  consensus[GenomicRanges::width(consensus) < minimum.width] <- consensus.small

  # Reduce the ranges
  consensus <- GenomicRanges::reduce(consensus, drop.empty.ranges = TRUE)

  consensus
}

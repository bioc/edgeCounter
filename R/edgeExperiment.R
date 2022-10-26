# Construct experiment which counts edge over a bunch of samples

# Given that you have the following:
#   1. A bunch of BAM files in which you want to count edges of the reads
#   2. A GRanges or GRangesList object containing the ranges that you want to count
# How to get an object of edge-counting results?


#' Extract edges from BAM files
#'
#' @description `edgeExtract()` is a simple wrapper of `bamEdges` to extract edges
#' of the BAM reads from multiple BAM files each of which represents one sample.
#' This step is IO and memory intensive. Do NOT support parallel computing.
#'
#' @param bam.paths  A named character vector containing paths to the BAM files.
#' Names are used as sample IDs.
#' @returns A compressed GRangesList containing edges from each of the BAM files.
#' IDs will be saved in `names()`.
#'
#' @export
edgeExtract <- function(bam.paths) {
  # Sanity check, especially given that this is an expensive operation
  #   require existence of sample IDs as names
  stopifnot(!is.null(names(bam.paths)))
  #   check file existence
  stopifnot(all(file.exists(bam.paths)))

  # Construct the GRangesList
  grl <- GenomicRanges::GRangesList(
    lapply(bam.paths, bamEdges), compress = TRUE
  )
  # Add names to the list
  names(grl) <- names(bam.paths)

  grl
}


#' Construct a RangedSummarizedExperiment with BAM paths and ranges of interest
#'
#' @description `edgeExperiment()` constructs a RangedSummarizedExperiment object containing
#' the edge-counting results given BAM file paths and a GRanges/GRangesList with
#' the target ranges for counting. You need to specify IDs of the samples.
#'
#' @param bam.paths A named character vector containing paths to the BAM files.
#' Names are used as sample IDs.
#' @param ranges A GRanges or GRangesList object containing the ranges for counting.
#' If GRanges is provided, each single range row is calculated individually.
#' If GRangesList is provided, each GRanges in the list is calculated individually
#' instead.
#' @returns A RangedSummarizedExperiment object containing the counting results.
#' One `assay` is saved, `counts`, containing the counting results.
#' `colData` will contain the sample IDs as row names but no column.
#' @seealso [SummarizedExperiment::RangedSummarizedExperiment]
#' @export
edgeExperiment <- function(bam.paths, ranges) {
  # Call edgeExtract to get a named list of counts first.
  # TODO
}


#' Construct a RangedSummarizedExperiment with edge-counting results.
#'
#' @description `edgeExperimentFromCounts()` constructs a RangedSummarizedExperiment
#' object containing the edge-counting results, given pre-extracted edges of
#' the samples and a GRanges/GRangesList with the target ranges for counting.
#'
#' @param edges Pre-extracted edges using [edgeExtract()].
#' @param ranges A GRanges or GRangesList object containing the ranges for counting.
#' If GRanges is provided, each single range row is calculated individually.
#' If GRangesList is provided, each GRanges in the list is calculated individually
#' instead.
#' @returns A RangedSummarizedExperiment object containing the counting results.
#' One `assay` is saved, `counts`, containing the counting results.
#' `colData` will contain the sample IDs as row names but no column.
#' @seealso [SummarizedExperiment::RangedSummarizedExperiment],
#' [edgeExtract()], [ecParallel()].
#' @details
#' # Details
#'
#' To extract edges of the samples use the `edgeExtract` function. This function
#' supports parallel computing as set up by [ecParallel()].
#' @export
edgeExperimentFromCounts <- function(edges, ranges) {
  # Only need to run `countEdgesCached` and not `bamEdges`.

  # Sanity check of edges - must be edgeExtract results or the same class
  stopifnot(inherits(edges, "CompressedGRangesList"))

  # count edges for each sample
  sampleIDs <- names(edges)
  counts <- ecParallelFunc(countEdgesCached,
                           ranges = list(ranges), edges = as.list(edges))
  counts <- lapply(counts, function(gr) GenomicRanges::mcols(gr)[["edge.counts"]])
  counts <- matrix(unlist(counts), ncol = length(counts))
  colnames(counts) <- sampleIDs

  # If `ranges` is a GRanges we will need to break it into a list
  if (inherits(ranges, "GRanges")) ranges <- breakGRanges(ranges)

  # construct the RSE object
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowRanges = ranges,
    colData = S4Vectors::DataFrame(sampleID = sampleIDs)
  )
}

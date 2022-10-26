# Main functions for counting edges

#' Get edges of BAM alignment reads
#'
#' @description
#' `bamEdges()` reads a BAM file and returns the 5' AND 3' ends of EACH
#' paired-end read PAIR as a `GRanges` object.
#'
#' @param bam.path Path to the BAM file. The BAM file is presumed to contain
#' read PAIRs.
#' @returns A `GRanges` object representing edges of BAM reads.
#'
#' @details
#' # Use cases
#'
#' * For ATAC sequencing, read edges are roughly transposition sites of the Tn5
#' transposase. To get the exact base-pair-precision transposition sites
#' you need to manually shift the read start/end positions.
#'
#' * Other cases?
#' @export
bamEdges <- function(bam.path){

  # Read in BAM file; only take necessary fields
  scanParam <- Rsamtools::ScanBamParam(
    what = c('rname', 'pos', 'isize'),
    # For a pair, only take one mate
    flag = Rsamtools::scanBamFlag(isMinusStrand = F)
    )
  bam <- Rsamtools::scanBam(bam.path, param = scanParam)[[1]]

  # Get GRange of the 5' and 3' ends for each read pair
  # each read PAIR will be converted to TWO entries in the GRange
  ends <- c(bam$pos, bam$pos + abs(bam$isize))
  # each entry is of width=1 reflecting the 5' and 3' ends of the PAIR
  GenomicRanges::GRanges(
    seqnames = c(bam$rname, bam$rname),
    ranges = IRanges::IRanges(
      start = ends,
      end = ends
      ),
    strand = NULL,
    seqinfo = NULL
  )
}


#' Count number of BAM alignment read edges for each range of a `GRanges/GRangesList` object
#'
#' @description `countEdgesCached()` counts number of BAM read pair edges for each
#' range of a given `GRanges` object. Compared to the exported `countEdges()`,
#' `countEdgesCached()` requires providing the pre-extracted edges using `bemEdges()`.
#'
#'
#' @param ranges A GRanges/GRangesList object containing the ranges to consider
#' @param edges A GRanges containing the edges. Use `bamEdges()` to get this.
#' @returns The input GRanges/GRangesList object with an appended metadata column
#' `edge.counts`.
#'
#' @seealso
#' For more detailed discussion on "edges"
#' of a BAM alignment file, see [bamEdges()].
#' Consult [GenomicRanges::GRanges] about `GRanges` objects.
countEdgesCached <- function(ranges, edges){

  counts <- GenomicRanges::countOverlaps(ranges, edges, type = "any")
  GenomicRanges::mcols(ranges)[["edge.counts"]] <- counts

  ranges
}


#' Count number of BAM alignment read edges for each range of a `GRanges/GRangesList` object
#'
#' @description `countEdges()` counts number of BAM read pair edges for each
#' range of a given `GRanges` object.
#'
#'
#' @param ranges A GRanges/GRangesList object containing the ranges to consider
#' @param bam.path Path to the BAM file.
#' @returns The input GRanges/GRangesList object with an appended metadata column
#' `edge.counts`
#'
#' @seealso
#' For more detailed discussion on "edges"
#' of a BAM alignment file, see [bamEdges()].
#' Consult [GenomicRanges::GRanges] about `GRanges` objects.
#' This function calls internal function [countEdgesCached()].
#' @export
countEdges <- function(ranges, bam.path){

  edges <- bamEdges(bam.path)
  countEdgesCached(ranges, edges)
}

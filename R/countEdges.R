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
#' # Usage cases
#'
#' * For ATAC sequencing, read edges are roughly transposition sites of the Tn5
#' transposase. However, to get the exact base-pair-precision
#' transposition sites you need to manually shift the read start/end positions.
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

#' Count number of BAM alignment read edges for each range of a `GRanges` object
#'
#' @description `countEdges()` counts number of BAM read pair edges for each
#' range of a given `GRanges` object.
#'
#'
#' @param ranges A GRanges object containing the ranges to consider
#' @param bam.path Path to the BAM file.
#' @returns The input GRanges object with an appended metadata column
#' `edge.counts`
#'
#' @seealso
#' For more detailed discussion on "edges"
#' of a BAM alignment file, see [bamEdges()].
#' Consult [GenomicRanges::GRanges] about `GRanges` objects.
#' @export
countEdges <- function(ranges, bam.path){
  # Given ranges as provided in `peaks` (GRanges object)
  # Count how many 5' and 3' ends are within each of the peak
  #   ... returns a numeric vector of length = length(peaks)
  edges <- bamEdges(bam.path)
  counts <- GenomicRanges::countOverlaps(ranges, edges, type = "any")
  GenomicRanges::mcols(ranges)[["edge.counts"]] <- counts

  ranges
}

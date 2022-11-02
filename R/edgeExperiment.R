# Construct experiment which counts edge over a bunch of samples

# Given that you have the following:
#   1. A bunch of BAM files in which you want to count edges of the reads
#   2. A GRanges or GRangesList object with the ranges that you want to count
# How to get an object of edge-counting results?


#' Extract edges from BAM files
#'
#' @description `edgeExtract()` is a simple wrapper of `bamEdges`
#' to extract edges of the BAM reads from multiple BAM files each of which
#' represents one sample. This step is IO and memory intensive.
#' Do NOT support parallel computing.
#'
#' @param bam.paths  A named character vector with paths to the BAM files.
#' Names are used as sample IDs.
#' @returns A compressed GRangesList with edges from each of the BAM files.
#' IDs will be saved in `names()`.
#' @examples
#' bam.paths <- system.file("extdata", "sample.near.Gapdh1.bam",
#'     package = "edgeCounter"
#' )
#' names(bam.paths) <- "sample"
#' edges <- edgeExtract(bam.paths) # Reads edges of one sample named "sample"
#' @export
edgeExtract <- function(bam.paths) {
    # Sanity check, especially given that this is an expensive operation
    #   require existence of sample IDs as names
    stopifnot(!is.null(names(bam.paths)))
    #   check file existence
    stopifnot(all(file.exists(bam.paths)))

    # Construct the GRangesList
    grl <- GenomicRanges::GRangesList(
        lapply(bam.paths, bamEdges),
        compress = TRUE
    )
    # Add names to the list
    names(grl) <- names(bam.paths)

    grl
}


#' Construct a RangedSummarizedExperiment with BAM paths and ranges of interest
#'
#' @description `edgeExperiment()` constructs a RangedSummarizedExperiment
#' object containing the edge-counting results given BAM file paths and
#' a GRanges/GRangesList with the target ranges for counting.
#' You need to specify IDs of the samples.
#'
#' @param bam.paths A named character vector containing paths to the BAM files.
#' Names are used as sample IDs.
#' @param ranges A GRanges or GRangesList with the ranges for counting.
#' If GRanges is provided, each single range row is calculated individually.
#' If GRangesList is provided, each GRanges in the list
#' is calculated individually instead.
#' @returns A RangedSummarizedExperiment containing the counting results.
#' One `assay` is saved, `counts`, containing the counting results.
#' `colData` will contain the sample IDs as row names but no column.
#' @seealso [SummarizedExperiment::RangedSummarizedExperiment]
#'
#' @examples
#' bam.paths <- system.file("extdata", "sample.near.Gapdh1.bam",
#'     package = "edgeCounter"
#' )
#' names(bam.paths) <- "sample"
#' ranges <- makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         7791000, 7794000, "chr2R", # inclusive of Gapdh1 locus
#'         80, 120, "chr2L", # random range
#'         100, 300, "chrX" # random range
#'     )
#' )
#' exp <- edgeExperiment(bam.paths, ranges) # SummarizedExperiment with counts
#' @export
edgeExperiment <- function(bam.paths, ranges) {
    # Sanity check
    stopifnot(!is.null(names(bam.paths)))
    stopifnot(inherits(ranges, "GRanges") | inherits(ranges, "GRangesList"))
    # Call edgeExtract to get a named list of counts first.
    edges <- edgeExtract(bam.paths)
    # Call edgeExperimentFromCounts to get the final result
    edgeExperimentFromCounts(edges, ranges)
}


#' Construct a RangedSummarizedExperiment with edge-counting results.
#'
#' @description `edgeExperimentFromCounts()` constructs a
#' RangedSummarizedExperiment object containing the edge-counting results,
#' given pre-extracted edges of the samples
#' and a GRanges/GRangesList with the target ranges for counting.
#'
#' @param edges Pre-extracted edges using [edgeExtract()].
#' @param ranges A GRanges or GRangesList containing the ranges for counting.
#' If GRanges is provided, each single range row is calculated individually.
#' If GRangesList is provided, each GRanges in the list
#' is calculated individually instead.
#' @returns A RangedSummarizedExperiment containing the counting results.
#' One `assay` is saved, `counts`, containing the counting results.
#' `colData` will contain the sample IDs as row names but no column.
#' @seealso [SummarizedExperiment::RangedSummarizedExperiment],
#' [edgeExtract()], [ecParallel()].
#' @details
#' # Details
#'
#' To extract edges of the samples use the `edgeExtract` function.
#' This function supports parallel computing as set up by [ecParallel()].
#'
#' ## colData
#'
#' sampleID = IDs provided as names of the `edges` input list.
#' edges.ir = Number of edges within the given ranges.
#' edges.total = Number of edges in the original `edges` input.
#'
#' @examples
#' bam.paths <- system.file("extdata", "sample.near.Gapdh1.bam",
#'     package = "edgeCounter"
#' )
#' names(bam.paths) <- "sample"
#' ranges <- makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         7791000, 7794000, "chr2R", # inclusive of Gapdh1 locus
#'         80, 120, "chr2L", # random range
#'         100, 300, "chrX" # random range
#'     )
#' )
#' edges <- edgeExtract(bam.paths) # Reads edges of one sample named "sample"
#' exp <- edgeExperimentFromCounts(edges, ranges) # You can use other ranges
#' @export
edgeExperimentFromCounts <- function(edges, ranges) {
    # Only need to run `countEdgesCached` and not `bamEdges`.

    # Sanity check of edges - must be edgeExtract results or the same class
    stopifnot(inherits(edges, "CompressedGRangesList"))

    # count edges for each sample
    sampleIDs <- names(edges)
    counts <- ecParallelFunc(countEdgesCached,
        ranges = list(ranges), edges = as.list(edges)
    )
    counts <- lapply(
        counts,
        function(gr) GenomicRanges::mcols(gr)[["edge.counts"]]
    )
    counts <- matrix(unlist(counts), ncol = length(counts))
    colnames(counts) <- sampleIDs

    # If `ranges` is a GRanges we will need to break it into a list
    #   TODO: THIS ARRANGEMENT SEEMS TO SLOW DOWN CODE
    #         SIGNIFICANTLY MORE THAN IT SHOULD BE.
    #   UPDATE: DOES NOT SEEM NECESSARY?
    # if (inherits(ranges, "GRanges")) ranges <- breakGRanges(ranges)

    # Gather normalization counts
    total.edges.ir <- apply(counts, 2, sum) # total edges in ranges
    total.edges <- vapply(edges, length, 0) # total edges

    # construct the RSE object
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = ranges,
        colData = S4Vectors::DataFrame(
            sampleID = sampleIDs,
            edges.ir = total.edges.ir,
            edges.total = total.edges
        )
    )
}


#' Plots variance ~ mean of range counts
#'
#' @description `varMeanPlot()` plots variance ~ mean of range counts given
#' edge-counting result `SummarizedExperiment` object.
#'
#' @param edgeExp Edge-counting results as constructed by [edgeExperiment()].
#' @param name Optional. Extra identifier of the dataset to be added on title.
#' @returns None. Side effect of generating a single plot.
#' @seealso [edgeExperiment()], [edgeExperimentFromCounts()].
#' @examples
#' character(0) # Function is only meaningful with real data; see vignette
#' @export
varMeanPlot <- function(edgeExp, name = NULL) {
    # Sanity check
    stopifnot(inherits(edgeExp, "SummarizedExperiment"))
    # Get counts data
    counts <- SummarizedExperiment::assay(edgeExp, "counts")
    # Retain ranges which have at least two NON-ZERO count values
    zeros <- apply(counts == 0, 1, sum)
    nonzeros <- ncol(counts) - zeros
    retain <- nonzeros >= 2
    message(sum(!retain), " ranges are removed.")
    counts <- counts[retain, ]
    # Compute mean and variance
    means <- apply(counts, 1, mean)
    vars <- apply(counts, 1, stats::var)
    # Plot the results
    graphics::plot.default(
        x = means, y = vars, main = paste("var ~ mean", name),
        xlab = "Mean", ylab = "Variance", log = "xy"
    )
    graphics::abline(a = 0, b = 1, untf = TRUE, col = "red", lty = 2)
}


#' Export edge-counting results
#'
#' @description `export.edgeExperiment()` exports the edge-counting results
#' as BED tracks for visualization. Counts data normalization method may be
#' specified to allow comparison among different samples.
#' Current implementation is relatively slow due to `IRanges::RleList`
#' performance issues.
#'
#' @param edgeExp A SummarizedExperiment object with counts assay data.
#' May be created by [edgeExperiment()] functions.
#' @param base.path Directory to save the exported files. Note that for
#' each sample one BED file will be exported as `{sampleID}.bed`.
#' @param norm Normalization method. Currently supported are
#' `epm_ir`, `epm` and `none`. Refer to Details for explanation.
#' @returns None. Only side effect of saving BED files.
#' @details
#' # Details
#'
#' Nature of the `edgeExperiment()` function family is to generate counts
#' data for downstream analysis. Apart from differential analysis, which is
#' supported by various packages such as `edgeR` and `DESeq2`, visualization
#' is important for getting intuition of the data and verifying the results.
#'
#' To support this, `export.edgeExperiment()` exports counts data into BED
#' tracks. Further, to allow comparison among different samples,
#' normalization methods as described below may be adopted.
#'
#' ## Normalization of counts
#'
#' Problem: different samples have different numbers of reads. Solution:
#' calculate 'edges per million in range (`epm_ir`)'. It divides all counts
#' by the total number of edges in the sample.
#'
#' This function also allows normalization against ALL edges in a sample
#' with the `epm` option. Refer to [edgeExperimentFromCounts()] for
#' more details.
#'
#' No normalization: `none`.
#'
#' ## BED columns
#'
#' Export follows standard 12-column BED format:
#'
#' chrom, chromStart, chromEnd, name, score, strand
#' thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
#'
#' name is `names` of the `rowRanges`
#' score is (optionally normalized) count values with one decimal place
#' strand is always "." (i.e., no strand)
#' itemRgb is always 255,0,0 (i.e., red color)
#' blockCount/Sizes/Starts will be set according to the ranges
#'
#' ## BED file track line
#'
#' One track line is written to first line of the BED file to
#' hopefully set up appropriate display settings. Refer to the
#' [UCSC BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and
#' [IGV forum](https://groups.google.com/g/igv-help/c/IyBzO5_j5E8).
#'
#' As of November 2022, we only found `igv.js` as a browser that supports all
#' features.
#'
#' @examples
#' \dontrun{
#' export.edgeExperiment(exp, ".", norm = "epm_ir")
#' }
#' @export
export.edgeExperiment <- function(edgeExp, base.path, norm = "epm_ir") {
    # Sanity check
    stopifnot(inherits(edgeExp, "SummarizedExperiment"))
    stopifnot(norm %in% c("epm_ir", "epm", "none"))
    # Create directory if not exist
    if (!dir.exists(base.path)) dir.create(base.path, recursive = TRUE)
    # Get counts (normalize if specified)
    counts <- SummarizedExperiment::assay(edgeExp, "counts")
    if (norm == "epm_ir") {
        totals <- SummarizedExperiment::colData(edgeExp)$edges.ir
        counts <- sweep(counts, 2, totals, "/") * 1E6 # per million
    }
    if (norm == "epm") {
        totals <- SummarizedExperiment::colData(edgeExp)$edges.total
        counts <- sweep(counts, 2, totals, "/") * 1E6 # per million
    }
    SummarizedExperiment::assay(edgeExp, "counts") <- counts
    # Export BEDs
    bed.cols <- rowRangesToBED(edgeExp)
    bed.standard.order <-
        c(
            "chrom", "chromStart", "chromEnd", "name", "score", "strand",
            "thickStart", "thickEnd", "itemRgb",
            "blockCount", "blockSizes", "blockStarts"
        )
    lapply(
        SummarizedExperiment::colnames(edgeExp),
        function(sampleID) {
            #  Add score column
            cols <- bed.cols
            score <- counts[, sampleID]
            score <- round(score, 1)
            cols$score <- score
            #  Rearrange order
            cols <- cols[bed.standard.order]
            #  Write "track line"
            tl <- getBEDtrackLine(name = sampleID, score = score)
            path <- file.path(base.path, paste0(sampleID, ".bed"))
            readr::write_lines(x = tl, file = path)
            #  Append BED data
            cols <- as.data.frame(cols)
            readr::write_delim(
                x = cols, file = path, delim = "\t",
                append = TRUE, col_names = FALSE
            )
        }
    )
    invisible()
}


#' Export edge distribution data for visualization
#'
#' @description `export.edgeExtract()` exports edges cached by the
#' `edgeExtract()` function for visualization. It generates bigwig files
#' for random access. This function relies on Bioconductor `rtracklayer`
#' and will NOT work on Windows.
#'
#' @param edgeEx A named list of GRanges object containing the BAM edges.
#' May be created by [edgeExtract()].
#' @param base.path Directory to save the exported files. Note that for
#' each sample one BIGWIG file will be exported as `{sampleID}.bw`.
#' @param bin.width Bin width. See Details.
#' @param bin.min Bin filter. See Details.
#' @param smooth.width If not NULL, smooth the edge positions. See details.
#' @param seqinfo Genome sequence info, a [GenomeInfoDb::Seqinfo] object.
#' @param norm Normalization method. Currently supported are
#' `epm_ir`, `epm` and `none`. Refer to Details for explanation.
#' @returns None. Only side effect of writing bigwig files.
#' @details
#' # Details
#'
#' Nature of the `edgeExperiment()` function family is to generate counts
#' data for downstream analysis. Apart from differential analysis, which is
#' supported by various packages such as `edgeR` and `DESeq2`, visualization
#' is important for getting intuition of the data and verifying the results.
#'
#' To support this, `export.edgeExtract()` exports BAM edges data to BIGWIG
#' tracks. The following sections explain exact procedure of this
#' function in order.
#'
#' ## Gathering genome size information
#'
#' To generate bigwig files genome lengths must be known. With
#' `seqinfo` parameter the user provides a [GenomeInfoDb::Seqinfo]
#' object containing the sequence lengths.
#'
#' ## Counting edges with fixed-width bins
#'
#' The whole genome is divided into fixed-width bins along each seqname.
#'
#' If `smooth.width` is set to a positive integer, edges expanded to
#' `width = smooth.width` centering at the original edge position.
#'
#' Then, edges are counted w.r.t. the bins using
#' [edgeExperimentFromCounts()].
#'
#' ## Normalization of counts
#'
#' Problem: different samples have different numbers of reads. Solution:
#' calculate 'edges per million (`epm`)'. It divides all counts
#' by the total number of edges counted in the ranges. Refer
#' to [edgeExperimentFromCounts()] for more details.
#'
#' No normalization: `none`.
#'
#' ## BIGWIG conversion with `rtracklayer`
#'
#' This function relies on `rtracklayer` for Bigwig generation. Refer to
#' [rtracklayer::BigWigFile] for more detail.
#'
#' @examples
#' \dontrun{
#' export.edgeExtract(exp, ".", norm = "epm_ir")
#' }
#' @export
export.edgeExtract <- function(edgeEx, base.path, bin.width = 20,
    bin.min = 1, smooth.width = NULL,
    seqinfo, norm = "epm") {
    # Sanity check - GRanges edgeEx, existent seqlengths, smooth
    stopifnot(
        inherits(edgeEx, "GRangesList"),
        all(vapply(edgeEx, function(x) inherits(x, "GRanges"), TRUE))
    )
    stopifnot(
        inherits(seqinfo, "Seqinfo"),
        all(!is.na(GenomeInfoDb::seqlengths(seqinfo)))
    )
    stopifnot(is.null(smooth.width) | smooth.width > 0)
    stopifnot(norm %in% c("epm", "none"))
    # Create directory if not exist
    if (!dir.exists(base.path)) dir.create(base.path, recursive = TRUE)
    # Making GRanges of fixed-width genomic bins
    bins <- GenomicRanges::tileGenome(seqinfo,
        tilewidth = bin.width,
        cut.last.tile.in.chrom = TRUE
    )
    # Smooth if specified
    if (!is.null(smooth.width)) {
        edgeEx <- S4Vectors::endoapply(edgeEx, function(gr) {
            GenomicRanges::resize(gr, smooth.width, fix = "center")
        })
    }
    # Counting with `edgeExperimentFromCounts`
    exp <- edgeExperimentFromCounts(edgeEx, bins)
    # Exporting bigwig with rtracklayer
    for (sampleID in SummarizedExperiment::colnames(exp)) {
        ranges <- bins
        score <- SummarizedExperiment::assay(exp, "counts")[, sampleID]
        # normalize score as specified
        if (norm == "epm") {
            factor <- SummarizedExperiment::colData(exp)$edges.total[sampleID]
            score <- score / factor * 1E6
        }
        GenomicRanges::mcols(ranges)$score <- score
        # filter ranges
        ranges <- ranges[score >= bin.min]
        # write bigwig using rtracklayer
        rtracklayer::export(
            object = ranges,
            con = file.path(base.path, paste0(sampleID, ".bw")),
            format = "bw"
        )
    }
}

# Useful utility functions that satisfy the following criteria:
# 1. Wide applicability for testing/internal usage purposes
# OR
# 2. Useful for end-user to perform certain tasks

#' Simple wrapper for making GRanges
#'
#' @description `makeGRanges()` makes an artificial GRanges based
#' on user input. Useful for testing purposes and vignette building.
#'
#' @param df A tibble/data.frame object with at least two columns: start, end.
#' Optionally provide seqnames (otherwise will take default value)
#' Optionally provide strand (otherwise no strand information, i.e., '*')
#'
#' @param seqnames.default What value should seqnames take
#' if `is.null(df$seqnames)`.
#' This is restricted to length = 1 character vector.
#'
#' @examples
#' makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         1, 100, "chr2L",
#'         80, 120, "chr2L",
#'         100, 300, "chrX"
#'     )
#' )
#'
#' makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end,
#'         1, 100,
#'         80, 120,
#'         100, 300
#'     )
#' )
#'
#' @returns A GRanges object
#' @export
makeGRanges <- function(df, seqnames.default = "chr2L") {
    stopifnot(is.data.frame(df)) # tibble is also a data.frame
    stopifnot(
        !is.null(df$start),
        !is.null(df$end)
    )
    suppressWarnings(
        result <-
            GenomicRanges::GRanges(
                seqnames = ifelse(rep(is.null(df$seqnames), times = nrow(df)),
                    rep(seqnames.default, times = nrow(df)),
                    df$seqnames
                ),
                ranges = IRanges::IRanges(
                    start = df$start,
                    end = df$end
                ),
                strand = df$strand
            )
    )

    result
}


#' Break down a `GRanges` into a `GRangesList` of single ranges.
#'
#' @description `breakGRanges()` breaks a `GRanges` object into a `GRangesList`
#' where each list element is a single range (i.e., row)
#' in the original `GRanges`.
#'
#' @param gr Input `GRanges`
#' @returns A `GRangesList` object
breakGRanges <- function(gr) {
    GenomicRanges::GRangesList(
        lapply(1:length(gr), function(i) gr[i, ]),
        compress = TRUE
    )
}


#' Reading BED and BED-like file with minimum dependency
#'
#' @description `importBED()` read in a BED file as a GRanges object.
#'
#' @param path Path to the BED file. See details for format of the BED file.
#' @param extra.cols A named list for information on extra columns. It should
#' contain two elements `names` and `types` for name and type of the columns.
#'
#' @details
#' # Details
#'
#' ## BED file format
#'
#' The standard BED format supported by this function is three-column. Refer to
#' \href{https://genome.ucsc.edu/FAQ/FAQformat.html}{UCSC FAQ} for details.
#'
#' ## Extra columns
#'
#' Apart from the standard BED format, you can supply arbitrary extra columns
#' by providing names and types of the columns in the `extra.cols`
#' optional argument.
#'
#' The extra column names must be a character vector.
#'
#' Types must follow the shorthand types used by [readr::read_delim()].
#'
#' ## Strand column
#'
#' According to UCSC documentation, `strand` is an optional column
#' at position 6. Therefore, 6th column of the file, if existent,
#' must specify strand info. Supported strand spec characters are
#' `.` (i.e., no strand), `+` and `-`.
#'
#' @returns A GRanges object
#' @examples
#' bed.path <- system.file("extdata", "importBED.example.bed",
#'     package = "edgeCounter"
#' )
#' importBED(bed.path)
#' @export
importBED <- function(path, extra.cols = NULL) {
    # Sanity checks of the extra column specs if provided
    if (!is.null(extra.cols)) {
        stopifnot(
            is.character(extra.cols$names),
            is.character(extra.cols$types)
        )
    }

    # Get the full column specs
    names <- c(
        c("seqnames", "start", "end"), # BED default
        extra.cols$names
    )
    types <- paste0(
        "cii", # BED default
        extra.cols$types
    )

    # Read in the BED file
    df <- readr::read_delim(
        file = path, delim = "\t",
        col_names = names, col_types = types
    )

    # Add strand information if provided
    if (ncol(df) >= 6) {
        # Extract strand
        strands <- df[[6]]
        # Sanity check
        stopifnot(
            is.character(strands),
            all(strands %in% c(".", "+", "-"))
        )
        # Change `.` (UCSC) to `*` (GenomicRanges)
        strands[strands == "."] <- "*"
        # Write strand information and make sure it is named `strand`.
        #   This `strand` name is used by `makeGRanges` function.
        df[[6]] <- strands
        names(df)[6] <- "strand"
    }

    # Convert into GRanges
    gr <- makeGRanges(df)

    # Adding the metadata columns
    df$start <- NULL
    df$end <- NULL
    df$seqnames <- NULL
    suppressWarnings({
        df$strand <- NULL # tibble warns if `strand` column is not existent.
    })
    GenomicRanges::mcols(gr) <- df

    gr
}

# Parallel computing supporting functions
#   Uses global options, refer to `edgeCounterOptions`


#' Set parallel computing options of edgeCounter
#'
#' @description `ecParallel()` sets global parallel computing options.
#'
#' @param PSOCK.names `name` of [parallel::makePSOCKcluster].
#' If `NULL` is provided, parallel computing feature will be turned off.
#'
#' @returns None.
#' @examples
#' ecParallel(6) # 6 workers, local machine
#' getOption("ecOptions")
#' ecParallel(NULL) # turn off PSOCK parallelization
#' getOption("ecOptions")
#' @export
ecParallel <- function(PSOCK.names) {
    # Sanity checking
    valid <- function(x) is.null(x) | ecOptions.valid()[["PSOCK.names"]](x)
    stopifnot(valid(PSOCK.names))

    # Set parallel options
    if (is.null(PSOCK.names)) {
        ecOptions(parallel = FALSE)
    } else {
        ecOptions(parallel = TRUE, PSOCK.names = PSOCK.names)
    }

    invisible()
}


#' Verify if an object is list-like
#'
#' @description `is.list.like()` checks if an object is 'list-like'.
#'
#' @param obj A single object
#'
#' @returns Boolean where TRUE means the object is 'list-like'.
#'
#' @details
#'
#' A 'list-like' object is anything that has the following operators:
#'
#' * `[[]]` double bracket subsetting, which should
#' return an object in the list. Integer subsetting must be supported.
#' * `length` which should return length of the list.
#'
#' @examples
#' gr <- makeGRanges(
#'     tibble::tribble(
#'         ~start, ~end, ~seqnames,
#'         1, 100, "chr2L",
#'         80, 120, "chr2L",
#'         100, 300, "chrX"
#'     )
#' )
#'
#' is.list.like(gr) # FALSE
#' is.list.like(GenomicRanges::GRangesList(list(gr, gr))) # TRUE
#'
#' @export
is.list.like <- function(obj) {
    tryCatch(
        {
            length(obj)
            obj[[length(obj)]]
            return(TRUE)
        },
        error = function(e) {
            return(FALSE)
        }
    )
}


#' Parallel execution of any function with arbitrary arguments.
#'
#' @description `ecParallelFunc()` runs any function with arbitrary arguments
#' in a parallel manner specified by `ecOptions()`.
#'
#' @param func Function to run, which must take arguments specified later.
#' @param benchmark If TRUE will message execution time. The time measured
#' includes all overheads such as PSOCK cluster initialization.
#' @param ... Arbitrary arguments of the function. All arguments must be
#' list-like' objects (see details).
#'
#' @returns A list of the execution results.
#' @details
#' # Arguments supplied
#'
#' `ecParallelFunc()` takes arbitrary arguments supplied as list-like objects
#' and feed consecutive combinations of the arguments.
#'
#' ## 'list-like' objects
#'
#' A list-like object must support the following operators.
#'
#' * `[[]]` double bracket subsetting, return an individual object in the list.
#' * `length` length accessor, return length of the list.
#'
#' # Behavior depending on `ecOptions`
#'
#' [ecOptions()] describes options of this package
#' including parallel computation.
#' If `parallel == FALSE`, `ecParallelFunc()` will be running
#' the combinations in serial. If `parallel == TRUE`,
#' `ecParallelFunc()` will initialize PSOCK cluster, execute computation,
#' stop the PSOCK cluster, and then return the results.
#'
#' @examples
#' ecParallelFunc(function(x, y) x + y, x = c(1, 2, 3), y = list(5))
#' @export
ecParallelFunc <- function(func, benchmark = FALSE, ...) {
    # If `benchmark` get time
    if (benchmark) time <- Sys.time()
    # Get the arguments
    args <- list(...)
    # Sanity checking
    stopifnot(is.function(func)) #   Must provide a function
    stopifnot(all(vapply(args, is.list.like, TRUE)))  #   list-like arguments
    #   Argument list lengths must be the same or 1
    arg.lens <- vapply(args, length, 1)
    arg.lens <- dimnames(table(arg.lens))[[1]]
    arg.lens <- arg.lens[arg.lens != "1"]
    stopifnot(length(arg.lens) %in% c(0, 1))
    # Decide on parallel computation
    par <- getOption("ecOptions")$parallel
    if (par) {
        # Start cluster
        par.names <- getOption("ecOptions")$PSOCK.names
        par.packages <- getOption("ecOptions")$parallel.packages
        cl <- parallel::makePSOCKcluster(par.names)
        # Load packages
        parallel::clusterCall(cl,
            function(packages) {
                lapply(packages, function(x) do.call(require, list(x)))
            },
            packages = par.packages
        )
        # Load function for computation
        parallel::clusterExport(cl, "func", envir = environment())
        pMApply <- function(...) {
            parallel::clusterMap(cl = cl, fun = func,
                                 .scheduling = "static", ...)}
    } else {
        pMApply <- function(...) mapply(FUN = func, SIMPLIFY = FALSE, ...)
    }
    # Apply function
    results <- pMApply(...)
    # Clean up
    if (par) {
        parallel::stopCluster(cl)
    }
    # If `benchmark` message elapsed time
    if (benchmark) {
        message(
            "Elapsed time ",
            format(round(difftime(Sys.time(), time, units = "secs"),
                digits = 1
            ))
        )
    }

    results
}

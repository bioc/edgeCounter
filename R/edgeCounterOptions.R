# Global options of edgeCounter package, used by various functions
#   The following aspects are considered currently:
#     parallel computation support (by PSOCK cluster)
#       `parallel` - BOOLEAN
#       `PSOCK.names` - integer or character, see `parallel::makeCluster`
#       The package always stop PSOCK clusters on exit of supported functions.
#   Check out zzz.R for initialization of the options.


#' Validation of global options
#'
#' @description `ecOptions.valid()` is internally used for global option
#' parameters sanity validation.
#'
#' @returns A named list where names are names of available global options and
#' values are functions for validating the option.
#' Validation functions must take in a single argument and return Boolean.
ecOptions.valid <- function() {
    list(
        parallel = function(x) is.logical(x) & length(x) == 1,
        PSOCK.names = function(x) is.numeric(x) | is.character(x),
        parallel.packages = function(x) is.character(x)
    )
}

#' Defaults of global options
#'
#' @description `ecOptions.default()` is internally used for global option
#' defaults
#'
#' @returns A named list of the default values for the global options.
ecOptions.default <- function() {
    list(
        parallel = FALSE,
        PSOCK.names = NA,
        parallel.packages = "edgeCounter"
    )
}


#' Global options for the edgeCounter package
#'
#' @description `ecOptions()` sets edgeCounter package global options.
#'
#' @param ... Parameter name-value pairs supported by the package. See details
#' for a list of global parameters. If no parameter is provided, do not change
#' current options.
#'
#' @returns Up-to-date edgeCounter options.
#'
#' @details
#'
#' Global options are accessible by `getOption("ecOptions")`.
#' # Supported options
#'
#' * `parallel`: a Boolean specifying whether to enable parallel computation.
#'
#' * `PSOCK.names`: PSOCK cluster `names` parameter.
#' See [parallel::makePSOCKcluster]. Not used if `parallel` is false.
#'
#' * `parallel.packages`: Which packages to load in PSOCK cluster nodes.
#'
#' ## Parallel computation
#'
#' This package always use PSOCK cluster as supported by the `parallel` package
#' for parallel computation. A 'diligent' approach for initialization
#' and clean-up of the PSOCK cluster is used, meaning that the PSOCK cluster is
#' initialized only when functions supporting parallel computation is called,
#' and the cluster is stopped whenever the called function exits.
#' @examples
#' ecOptions() # Returns current option
#' ecOptions(parallel = FALSE) # Turns off parallel computation
#' @export
ecOptions <- function(...) {
    # INTERNAL NOTE - WHEN INITIALIZING THE PACKAGE, FEED NO ARGUMENT IS FINE.

    # Fetch arguments into a list
    args <- list(...)

    # If arguments are provided, we validate it
    if (length(args) != 0) {
        # Sanity check
        supported <- ecOptions.valid()
        #   All arguments must be supported
        stopifnot(all(names(args) %in% names(supported)))
        #   All argument values must be valid
        stopifnot(all(vapply(
            names(args),
            function(n) supported[[n]](args[[n]]),
            TRUE
        )))
    }

    # Set the global options
    #   Fetch current options
    opt <- getOption("ecOptions", default = list())
    #   Change options
    for (n in names(args)) opt[[n]] <- args[[n]]
    #   Fill undefined options with defaults
    def <- ecOptions.default()
    for (n in names(def)) if (!(n %in% names(opt))) opt[[n]] <- def[[n]]
    #   Update options
    options(ecOptions = opt)

    opt
}

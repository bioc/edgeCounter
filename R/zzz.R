# Side effects of upon package loading and unloading

.onLoad <- function(libname, pkgname) {
    # Initialize global options
    ecOptions()

    invisible()
}

.onUnload <- function(libpath) {
    # Reverse side effects
    options(ecOptions = NULL)

    invisible()
}

require_suggested_package <- function(pkg, ver = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE) ||
      (!is.null(ver) && utils::packageVersion(pkg) < ver)) {
    stop("Please install ",
         if (!is.null(ver)) paste("at least version", ver, "of "),
         "the ", pkg, " package.", call. = FALSE)
  }
}

is_family <- function(x) {
  inherits(x, "family")
}

#' Check if family is binomial/bernoulli
#' @noRd
#' @param x User's specified `family` argument.
#' @return `TRUE` or `FALSE` if no error.
family_is_binomial <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  if (!is.character(x) && !is_family(x)) {
    stop("Model family must be a string or family object", call. = FALSE)
  }
  if (is_family(x)) {
    x <- x$family
  }
  x %in% c("binomial", "bernoulli")
}
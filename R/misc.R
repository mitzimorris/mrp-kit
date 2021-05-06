`%||%` <- function (a, b) {
  if (is.null(a) || isTRUE(is.na(a))) b else a
}

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

#' Create weighted estimates using the survey package
#' @noRd
#' @param fit_obj  An [R6][R6::R6Class] SurveyFit object
#' @param outcome The variable we are estimating
#' @param by The grouping variables
#' @return A table of size nlevels by 3 (level, estimate, sd).
create_wtd_ests <- function(fit_obj, outcome, by = NULL) {
  weights <- fit_obj$map()$samp_obj()$weights()
  if (is.null(weights)) {
    stop("Sample weights must be present", call. = FALSE)
  }
  design <- fit_obj$map()$samp_obj()$design()
  merged_data <- merge(fit_obj$map()$samp_obj()$mapped_data(),
                       fit_obj$map()$samp_obj()$survey_data()[c(outcome,".key")],
                       by = ".key")
  svy_dsn <- do.call(survey::svydesign, c(design, list(weights = weights, data = merged_data)))
  if (is.null(by)) {
    wtd_ests <- survey::svymean(stats::as.formula(paste0(c('~',outcome), collapse = "")), design = svy_dsn)
    wtd_ests <- data.frame(wtd_ests)
    rownames(wtd_ests) <- levels(merged_data[[outcome]])
    wtd_ests <- wtd_ests[seq(2, dim(wtd_ests)[1], 2),]
    colnames(wtd_ests) <- c("mean", "sd")
  } else {
    wtd_ests <- survey::svyby(
      formula = stats::as.formula(paste0("~",outcome)),
      by = stats::as.formula(paste0("~", by)),
      design = svy_dsn,
      survey::svymean
    )
    wtd_ests <- data.frame(wtd_ests)[, c(1,3,5)]
    colnames(wtd_ests) <- c(by, "mean", "sd")
  }
  wtd_ests
}


#' Generate approximate samples of predicted probabilities from a binomial
#' `glmerMod` object.
#'
#' @export
#' @param object A binomial model fit using `lme4::glmer()`.
#' @param newdata The data to use.
#' @param nsamples Number of samples to generate.
#' @return A matrix with `nrow(newdata)` rows and `nsamples` columns.
#'
sim_posterior_probs  <- function(object, newdata, nsamples = 4000) {
  require_suggested_package("merTools")
  if (!inherits(object, "glmerMod")) {
    stop("Object must have class 'glmerMod'.", call. = FALSE)
  }
  if (stats::family(object)$family != "binomial") {
    stop("Model family must be binomial.", call. = FALSE)
  }
  predict_nd <-
    merTools::predictInterval(
      merMod = object,
      newdata = newdata,
      level = 0.95,
      n.sims = nsamples,
      stat = "median",
      type = "probability",
      include.resid.var = TRUE,
      returnSims = TRUE,
      .parallel = FALSE
    )
  stats::plogis(attr(predict_nd, "sim.results"))
}


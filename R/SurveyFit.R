#' SurveyFit
#'
#' @name SurveyFit
#' @export
#' @description An [R6][R6::R6Class] `SurveyFit` object stores a fitted model
#'   object and provides methods for generating predicted probabilities for all
#'   poststrat cells, generating population and group estimates, and visualizing
#'   results.
#' @inherit SurveyMap examples
#'
SurveyFit <- R6::R6Class(
  classname = "SurveyFit",
  private = list(
    map_ = NULL,
    fit_ = NULL
  ),
  public = list(

    #' @description Create a new `SurveyFit` object. This method is called
    #'   internally by the `$fit()` method of the [`SurveyMap`] object and does
    #'   not need to be called directly by the user.
    #' @param fit A fitted model object.
    #' @param map A [`SurveyMap`] object.
    initialize = function(fit, map) {
      private$fit_ <- fit
      private$map_ <- map
      invisible(self)
    },

    #' @description Access the fitted model object
    fit = function() {
      private$fit_
    },

    #' @description Access the SurveyMap object
    map = function() {
      private$map_
    },

    #' @description Call the fitted model object's print method
    #' @param ... Optional arguments to pass the print method.
    print = function(...) {
      print(private$fit_, ...)
      invisible(self)
    },

    #' @description Use fitted model to add predicted probabilities to post-stratification dataset.
    #' @param fun The function to use to generate the predicted probabilities.
    #'   This should only be specified if you used a model fitting function
    #'   not natively supported by \pkg{mrpkit}.
    #'   For models fit using \pkg{rstanarm}, \pkg{brms}, or \pkg{lme4}, `fun`
    #'   is handled automatically. If `fun` is specified then:
    #'   * the first argument should be the fitted model object
    #'   * the second argument should be the poststratification data frame
    #'   * it can take an arbitrary number of other arguments
    #'   * the returned object should match the specifications in the 'Returns'
    #'    section below in order to be compatible with subsequent methods
    #' @param ... Arguments other than the fitted model object and
    #'   poststratification data frame to pass to `fun`.
    #' @return A matrix with rows corresponding to poststratification cells and
    #'   columns corresponding to posterior samples (or approximate ones
    #'   in the case of \pkg{lme4} models).
    population_predict = function(..., fun = NULL) {
      args <- list(...)
      if (!is.null(args$newdata) && is.null(fun)) {
        stop("The 'newdata' argument should not be specified.",
             call. = FALSE)
      }
      if (is.null(private$map_$poststrat_data())) {
        stop("Post-stratification data not found. ",
             "Please call the tabulate() method before fitting a model.",
             call. = FALSE)
      }
      poststrat <- private$map_$poststrat_data()

      if (is.null(fun)) {
        if ("stanreg" %in% class(private$fit_)) {
          require_suggested_package("rstanarm", "2.21.0")
          return(
            t(suppressMessages(rstanarm::posterior_linpred(
              object = private$fit_,
              newdata = poststrat,
              transform = TRUE,
              ...
            )))
          )
        } else if ("brmsfit" %in% class(private$fit_)) {
          require_suggested_package("brms")
          return(
            t(brms::posterior_epred(
              object = private$fit_,
              newdata = poststrat,
              dpar = "mu",
              allow_new_levels = TRUE,
              sample_new_levels =
                if (!is.null(args$sample_new_levels)) args$sample_new_levels
                else "gaussian",
              ...
            ))
          )
        } else if ("glmerMod" %in% class(private$fit_)) {
          require_suggested_package("lme4")
          return(
            sim_posterior_probs(
              object = private$fit_,
              newdata = poststrat,
              ...
            )
          )
        } else {
          stop("Custom population_predict method required. Please specifiy 'fun'.",
               call. = FALSE)
        }
      } else {
        fun <- match.fun(fun)
        fun(fitted_model, poststrat, ...)
      }
    },

    #' @description Aggregate estimates to the population level or by level of a grouping variable
    #' @param poststrat_estimates The object returned by `population_predict`.
    #' @param by Optionally a string specifying a grouping variable. If
    #'   specified the aggregation will happen by level of the named variable.
    #'   If not specified population-level estimates will be computed.
    #' @return A data frame. If `by` is not specified then the data frame will
    #'   have number of rows equal to the number of posterior draws. If `by` is
    #'   specified the data frame will have number of rows equal to the number
    #'   of posterior draws times the number of levels of the `by` variable,
    #'   and there will be an extra column indicating which level of the `by`
    #'   variable each row corresponds to.
    aggregate = function(poststrat_estimates, by = NULL) {
      poststrat_data <- private$map_$poststrat_data()
      if (!is.null(by)) {
        if (length(by) != 1) {
          stop("Currently only one variable can be named in 'by'.", call. = FALSE)
        }
        rotate_levels <- levels(private$map_$samp_obj()$mapped_data()[, by])
        out <- expand.grid(by = rotate_levels, draw = 1:ncol(poststrat_estimates), value = NA)
        colnames(out)[1] <- by
        for (focus_level in rotate_levels){
          level_loc <- poststrat_data[by] == focus_level
          out[out[by] == focus_level, "value"] <-
            apply(poststrat_estimates[level_loc, ], 2, function(x) sum(poststrat_data$N_j[level_loc]*x)/sum(poststrat_data$N_j[level_loc]))
        }
      } else {
        out <- data.frame(value = apply(poststrat_estimates, 2, function(x) sum(poststrat_data$N_j*x)/sum(poststrat_data$N_j)))
      }
      structure(out, class = c(class(out), "mrp_aggregate"))
    },

    plot = function(aggregated_estimates, weights = TRUE) {
      if (dim(aggregated_estimates)[2] > 2){
        focus_var <- colnames(aggregated_estimates)[1]
        which_q <- private$map_$item_map()[[focus_var]]$col_names()[1]
        svy_q <- private$map_$samp_obj()$questions()[[which_q]]
        gg <- ggplot2::ggplot(aggregated_estimates) +
          ggplot2::aes(x = .data[[focus_var]], y = .data[["value"]]) +
          ggplot2::geom_violin(fill = "darkblue", alpha = .3) +
          ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
          ggplot2::xlab(svy_q)
      } else {
        model_fit <- private$fit_
        lhs_var <- as.character(formula(model_fit))[[2]]
        svy_q <- private$map_$samp_obj()$questions()[[lhs_var]]
        gg <- ggplot2::ggplot(aggregated_estimates) +
          ggplot2::aes(x = .data[["value"]], y = ggplot2::after_stat(scaled)) +
          ggplot2::geom_density(fill = "darkblue", alpha = .3, ) +
          ggplot2::scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
          ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
          ggplot2::xlab(svy_q)
      }

      if (weights) {
        model_fit <- private$fit_
        lhs_var <- as.character(formula(model_fit))[[2]]
        if (dim(aggregated_estimates)[2] > 2) {
          by_var <- colnames(aggregated_estimates)[1]
          wtd_ests <- create_wtd_ests(self, lhs_var, by=by_var)
          gg <- gg +
            ggplot2::geom_point(data = wtd_ests, ggplot2::aes(x= .data[[by_var]], y = .data[["mean"]])) +
            ggplot2::geom_errorbar(
              data = wtd_ests,
              ggplot2::aes(x = .data[[by_var]],
                           ymin = .data[["mean"]] - 1.96*.data[["sd"]],
                           ymax = .data[["mean"]] + 1.96*.data[["sd"]]),
              inherit.aes = FALSE, alpha = .5)
        } else {
          wtd_ests <- create_wtd_ests(self, lhs_var)
          gg <- gg +
            ggplot2::geom_vline(data = wtd_ests, ggplot2::aes(xintercept = .data[["mean"]])) +
            ggplot2::annotate("rect",
                              xmin = wtd_ests$mean - 1.96*wtd_ests$sd, xmax = wtd_ests$mean + 1.96*wtd_ests$sd,
                              ymin = 0, ymax = 1,
                              alpha = .5, fill = "grey"
            )
        }
      }
      gg
    },

    #' @description Summarize MRP estimates.
    #' @details This method summarizes the MRP estimates and not the fitted
    #'   model. To instead summarize just the fitted model first extract
    #'   it using the `fit` method and then use `summary` to call the method
    #'   from the package that fit the model, i.e., `summary(SurveyFit$fit())`.
    #'
    #'   The efficiency of this method depends greatly on whether
    #'   `population_estimates` and `group_estimates` are specified by the user
    #'   or whether they are unspecified and thus must be computed internally.
    #'   If they are not specified then both the `population_predict` and
    #'   `aggregate` methods will be called to produce the necessary estimates
    #'   to summarize.
    #'
    #' @param population_estimates Optionally, population estimates returned by
    #'   the `aggregate` method. If not provided this is regenerated internally,
    #'   which may be slow for large models and data.
    #' @param group_estimates Optionally, group estimates returned by the
    #'   `aggregate` method. This can be a single data frame returned by
    #'   `aggregate` or a list of such data frame (to summarize by multiple
    #'   different variables). If not provided this is regenerated internally
    #'   (if `by` is specified), which may be slow for large models and data.
    #' @param ... Arguments passed to `print`, e.g. `digits`.
    #' @param by Character vector of variable names. If `group_estimates` is not
    #'   provided then `by` is used to specify which variables to summarize by.
    #' @param stats A named list of functions that each return a single value. The list
    #' names are used as column names in the summarized output. The list values
    #' must be names of functions or functions themselves. For example:
    #'
    #'     stats = list("mean", "std_dev" = sd, p70 = function(x) quantile(x, 0.7))
    #'
    #' will result in computing the mean (labelled "mean" in output, inferred because
    #' the function name was given as a string), the standard deviation (labelled "std_dev"
    #' in the output) and the 70th percentile (labelled "p70" in the output).
    #' The default is to compute just the mean and standard deviation.
    #'
    #' @return A named list containing the following components:
    #'   * `population`: A data frame with one row and columns `mean` and `sd`.
    #'   * `grouped`: A list with one data frame per grouping variable. Each
    #'   data frame has a column for the level of the grouping variable and also
    #'   columns `mean` and `sd`.
    #'
    #'   The list has an extra class `"mrp_summary"` with a custom `print` method.
    summary = function(population_estimates,
                       group_estimates,
                       ...,
                       by = NULL,
                       stats = list("mean", "sd")) {

      stats <- validate_stats(stats)
      stat_names <- names(stats)

      poststrat_estimates <- NULL
      if (missing(population_estimates)) {
        poststrat_estimates <- self$population_predict()
        population_estimates <- self$aggregate(poststrat_estimates)
      } else if (!inherits(population_estimates, "mrp_aggregate") ||
                 !is.data.frame(population_estimates) ||
                 !identical(colnames(population_estimates), "value")) {
        stop("If specified 'population_estimates' must be a data frame returned ",
             "by the aggregate method.", call. = FALSE)
      }

      if (!missing(group_estimates)) {
        if (is.data.frame(group_estimates)) {
          group_estimates <- list(group_estimates)
        }
        if (!all(sapply(group_estimates, inherits, "mrp_aggregate")) ||
            !all(sapply(group_estimates, is.data.frame))) {
          stop("If specified 'group_estimates' must be a data frame returned ",
               "by the aggregate method or a list of such data frames.", call. = FALSE)
        }
        if (!is.null(by)) {
          warning("'by' is ignored if 'group_estimates' is specified.", call. = FALSE)
          by <- NULL
        }
      } else if (!is.null(by)) {
        if (is.null(poststrat_estimates)) {
          poststrat_estimates <- self$population_predict()
        }
        group_estimates <- lapply(by, function(x) self$aggregate(poststrat_estimates, by = x))
      } else {
        group_estimates <- list()
      }

      population_summary <- population_estimates %>%
        dplyr::summarize (
          x = sapply(stats, function(f) f(.data$value)),
          name = stat_names
        ) %>%
        dplyr::mutate(id = "fake_id") %>%
        stats::reshape(direction = "wide", idvar = "id", timevar = "name")
      population_summary <- population_summary[, -1, drop = FALSE]
      colnames(population_summary) <- stat_names

      group_summary <- lapply(seq_along(group_estimates), function(j) {
        by_var <- colnames(group_estimates[[j]])[1]
        out <- group_estimates[[j]] %>%
          dplyr::group_by(.data[[by_var]]) %>%
          dplyr::summarize (
            x = sapply(stats, function(f) f(.data$value)),
            name = stat_names
          ) %>%
          as.data.frame() %>%
          stats::reshape(direction = "wide", idvar = by_var, timevar = "name")
        colnames(out)[-1] <- stat_names
        out
      })
      names(group_summary) <- sapply(group_summary, function(x) colnames(x)[1])
      out <- list(population = population_summary, grouped = group_summary)
      structure(out, class = c(class(out), "mrp_summary"), ...)
    }
  )
)

# Print method for objects generated by SurveyFit$summary()
#' @export
print.mrp_summary <- function(x, ...) {
  args <- list(...)
  args$digits <- args$digits %||% attributes(x)$digits
  args$row.names <- args$row.names %||% FALSE

  cat("\nPopulation estimate:\n")
  args$x <- x$population
  do.call(print, args)

  for (j in seq_along(x$grouped)) {
    cat("\nEstimates by ", names(x$grouped)[j], ":\n", sep = "")
    args$x <- x$grouped[[j]]
    do.call(print, args)
  }
  invisible(x)
}

validate_stats <- function(stats) {
  # take list names or name of function if specified as a string
  # otherwise error
  stat_names <- names(stats)
  if (is.null(stat_names)) {
    stat_names <- vector("character", length(stats))
  }
  for (j in seq_along(stat_names)) {
    if (!nzchar(stat_names[j])) {
      if (is.character(stats[[j]])) {
        stat_names[j] <- stats[[j]]
      } else {
        stop("Name must be provided for element ", j , " of 'stats'.", call. = FALSE)
      }
    }
  }
  stats <- lapply(stats, match.fun)
  stats <- setNames(stats, stat_names)
  for (j in seq_along(stats)) {
    if (length(stats[[j]](1:2)) != 1) {
      stop("Function '", stat_names[j], "' returns more than one value.", call. = FALSE)
    }
  }
  stats
}


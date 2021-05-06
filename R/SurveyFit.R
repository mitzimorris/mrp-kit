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
    #'   This should only be specified if using a custom model fitting function.
    #'   For models fit using \pkg{rstanarm}, \pkg{brms}, or \pkg{lme4}, `fun`
    #'   is handled automatically. If `fun` is a custom function then the first
    #'   argument should take in the fitted model object and the second argument
    #'   should take in the poststratification (`newdata`) data frame. The
    #'   function must return a matrix with rows corresponding to the columns of
    #'   the poststratification data and columns corresponding to simulations.
    #' @param ... Arguments other than the fitted model and `newdata` data frame
    #'   to pass to `fun`.
    #' @return A matrix with rows corresponding to poststrat cells and columns
    #'   corresponding to posterior samples.
    #'
    population_predict = function(fun = NULL, ...) {
      args <- list(...)
      if (!is.null(args$newdata)) {
        stop("The 'newdata' argument should not be specified.",
             call. = FALSE)
      }
      if (is.null(private$map_$poststrat_data())) {
        stop("Post-stratification data not found. ",
             "Please call the tabulate() method before fitting a model.",
             call. = FALSE)
      }
      poststrat <- private$map_$poststrat_data()

      if (is.null(args$fun)) {
        if ("stanreg" %in% class(private$fit_)){
          require_suggested_package("rstanarm", "2.21.0")
          return(
            t(suppressMessages(rstanarm::posterior_linpred(
              object = private$fit_,
              newdata = poststrat,
              transform = TRUE,
              ...
            )))
          )
        }
        if ("brmsfit" %in% class(private$fit_)){
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
        }
        if ("glmerMod" %in% class(private$fit_)) {
          require_suggested_package("lme4")
          return(
            sim_posterior_probs(
              object = private$fit_,
              newdata = poststrat,
              ...
            )
          )
        }
      } else {
        fun <- match.fun(fun)
        fun(fitted_model, poststrat, ...)
      }
    },
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
      out
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

    #' @description Summarize MRP estimates. To instead summarize the
    #'   fitted model itself first extract it using the `fit` method and then
    #'   use `summary` to call the method from the package that fit the model,
    #'   i.e., `summary(SurveyFit$fit())`.
    #' @param poststrat_estimates The object returned by the
    #'   `population_predict` method. If not provided this is regenerated
    #'   internally which may take some time for very large models and data. For
    #'   models that require a custom `population_predict` method (models not
    #'   fit using \pkg{rstanarm}, \pkg{brms}, or \pkg{lme4})
    #'   `poststrat_estimates` can't be automatically regenerated and so the
    #'   argument must be specified.
    #' @param by Optionally, a character vector of variable names. If not
    #'   specified only the population estimate will be summarized. If specified
    #'   then estimates for each level of the variable in `by` will be
    #'   summarized.
    #' @param ... Arguments passed to `print`, e.g. `digits`.
    #' @return A named list containing the following components:
    #'   * `population`: data frame with one row and columns `mean` and `sd`.
    #'   * `grouped`: data frame with a column for the level of the grouping
    #'   variable and also columns `mean` and `sd`.
    summary = function(poststrat_estimates, by = NULL, ...) {
      if (missing(poststrat_estimates)) {
        poststrat_estimates <- self$population_predict()
      }
      popn_ests <- self$aggregate(poststrat_estimates)
      popn_summary <- data.frame(
        mean = mean(popn_ests$value),
        sd = stats::sd(popn_ests$value)
      )
      if (!is.null(by)) {
        group_ests <- self$aggregate(poststrat_estimates, by = by)
        group_summary <-
          group_ests %>%
          dplyr::group_by(.data[[by]]) %>%
          dplyr::summarise(mean = mean(.data$value), sd = stats::sd(.data$value)) %>%
          as.data.frame()
      } else {
        group_summary <- NULL
      }
      structure(
        list(population = popn_summary, grouped = group_summary),
        class = c("mrp_summary", "list"),
        ...
      )
    }
  )
)

# Print method for objects generated by SurveyFit$summary()
#' @export
print.mrp_summary <- function(x, ...) {
  dots <- list(...)
  atts <- attributes(x)
  cat("\nPopulation estimate:\n")
  print(x$population, row.names = FALSE, digits = dots$digits %||% atts$digits, ...)

  if (!is.null(x$grouped)) {
    cat("\nEstimates by group:\n")
    print(x$grouped, row.names = FALSE, digits = dots$digits %||% atts$digits, ...)
  }
  invisible(x)
}

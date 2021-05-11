#' SurveyMap
#'
#' @name SurveyMap
#' @export
#'
#' @description An [R6][R6::R6Class] `SurveyMap` object holds the mapping
#' between a set of items in a survey and a population dataset.
#' The label is the item label in each dataset and the values
#' is a list of all possible values.  The values for the survey
#' and population must be aligned, i.e., the lists must have the
#' same number of elements and the values at index i in each list
#' are equivalent.  If there is a meaningful ordering over the values,
#' they should be listed in that order, either descending or ascending.
#'
#' @examples
#' library(dplyr)
#'
#' head(feline_survey)
#' feline_prefs <- SurveyData$new(
#'   data = feline_survey,
#'   questions = list(
#'     age1 = "Please identify your age group",
#'     gender = "Please select your gender",
#'     pet_own = "Which pet do you own?",
#'     y = "Response"
#'   ),
#'   responses = list(
#'     age1 = levels(feline_survey$age1),
#'     gender = levels(feline_survey$gender),
#'     pet_own = levels(feline_survey$pet_own),
#'     y = c("no","yes")
#'   ),
#'   weights = feline_survey$wt
#' )
#' feline_prefs$print()
#'
#' head(approx_popn)
#' popn_obj <- SurveyData$new(
#'   data = approx_popn,
#'   questions = list(
#'     age2 = "Which age group are you?",
#'     gender = "Gender?",
#'     pet_pref = "Which pet would you like to own?"
#'   ),
#'   # order doesn't matter (gender before age2 here) because
#'   # the list has the names of the variables
#'   responses = list(
#'     gender = levels(approx_popn$gender),
#'     age2 = levels(approx_popn$age2),
#'     pet_pref = levels(approx_popn$pet_pref)
#'   ),
#'   weights = approx_popn$wt
#' )
#' popn_obj$print()
#'
#' q1 <- SurveyQuestion$new(
#'   name = "age",
#'   col_names = c("age1","age2"),
#'   values_map = list(
#'     "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
#'     "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
#'   )
#' )
#' q2 <- SurveyQuestion$new(
#'   name = "pet",
#'   col_names = c("pet_own","pet_pref"),
#'   values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
#' )
#' q3 <- SurveyQuestion$new(
#'   name = "gender",
#'   col_names = c("gender","gender"),
#'   values_map = data.frame("male" = "m","female" = "f", "nonbinary" = "nb")
#' )
#'
#' # Create SurveyMap object
#' # can add all questions at once or incrementally
#' ex_map <- SurveyMap$new(samp_obj = feline_prefs, popn_obj = popn_obj, q1)
#' print(ex_map)
#' ex_map$validate()
#' ex_map$add(q3)
#' print(ex_map)
#' ex_map$delete(q3)
#' print(ex_map)
#' ex_map$add(q3)
#' ex_map$delete("gender")
#' print(ex_map)
#' ex_map$add(q2)
#' print(ex_map)
#' ex_map$replace(q1,q3)
#' print(ex_map)
#' ex_map$add(q1)
#' print(ex_map)
#' ex_map$validate()
#' ex_map$mapping()
#' ex_map$tabulate("age") # Just use age in the poststrat matrix
#' ex_map$tabulate() # Use all variables in the map
#'
#'
#' #' Example rstanarm usage
#' #' Returns a SurveyFit object
#' fit_1 <- ex_map$fit(
#'   fun = rstanarm::stan_glmer,
#'   formula = y ~ (1|age) + (1|gender),
#'   family = "binomial",
#'   refresh = 100,
#'   cores = 2
#' )
#'
#' \dontrun{
#' # Example brms usage
#' # Returns a SurveyFit object
#' # (not run because requires compilation)
#' fit_2 <- ex_map$fit(
#'   fun = brms::brm,
#'   formula = y ~ (1|age) + (1|gender),
#'   family = "bernoulli",
#'   refresh = 100,
#'   cores = 2
#' )
#' }
#'
#' # predicted probabilities
#' # returns matrix with rows for poststrat cells, cols for posterior draws
#' poststrat_estimates <- fit_1$population_predict()
#'
#' # estimates by age and gender
#' estimates_by_age <- fit_1$aggregate(poststrat_estimates, by = "age")
#' head(estimates_by_age)
#'
#' estimates_by_gender <- fit_1$aggregate(poststrat_estimates, by = "gender")
#' head(estimates_by_gender)
#'
#' # plot estimates by age and gender
#' fit_1$plot(estimates_by_age)
#' fit_1$plot(estimates_by_gender)
#'
#' # population estimate
#' estimates_popn <- fit_1$aggregate(poststrat_estimates)
#' mean(estimates_popn$value)
#'
#' # plot population estimate
#' fit_1$plot(estimates_popn)
#'
#' # summarize MRP estimates numerically
#' fit_1$summary()
#' fit_1$summary(by = c("age", "gender"), digits = 2)
#'
#' # summarize MRP estimates using precomputed aggregated estimates
#' # this is much faster for large models/data
#' fit_1$summary(estimates_popn, digits = 2)
#' fit_1$summary(estimates_popn, estimates_by_age, digits = 2)
#' fit_1$summary(estimates_popn, list(estimates_by_age, estimates_by_gender), digits = 2)
#'
#' # specify the summary stats and names
#' my_stats <- list("avg" = mean, "sdt_dev" = sd, "p70" = function(x) quantile(x, 0.7))
#' fit_1$summary(by = "age", stats = my_stats, digits = 2)
#' fit_1$summary(by = "age", stats = list("median", "mad"))
#'
#' # summarize fitted model itself (not the MRP estimates)
#' summary(fit_1$fit())
#'
SurveyMap <- R6::R6Class(
  classname = "SurveyMap",
  private = list(
    samp_obj_ = NULL,
    popn_obj_ = NULL,
    item_map_ = list(),
    poststrat_data_ = data.frame(NULL)
  ),
  public = list(

    #' @description Create a new SurveyMap object
    #' @param samp_obj The [SurveyData] object corresponding to the sample data.
    #' @param popn_obj The [SurveyData] object corresponding to the population data.
    #' @param ... [SurveyQuestion] objects.
    initialize = function(samp_obj, popn_obj, ...) {
      if (!inherits(samp_obj, "SurveyData")) {
        stop("samp_obj must be a SurveyData object.", call. = FALSE)
      }
      if (!inherits(popn_obj, "SurveyData")) {
        stop("popn_obj must be a SurveyData object.", call. = FALSE)
      }

      private$item_map_ <- list(...)
      for (i in seq_along(private$item_map_)) {
        names(private$item_map_)[i] <- private$item_map_[[i]]$name()
      }

      private$samp_obj_ <- samp_obj
      private$popn_obj_ <- popn_obj
      invisible(self)
    },

    #' @description Print a summary of the mapping.
    #' @param ... Currently ignored.
    print = function(...) {
      if (length(private$item_map_) > 0) {
        for (i in 1:length(private$item_map_)) {
          cat("==============",'\n')
          cat(private$item_map_[[i]]$col_names()[1], "=",
              private$item_map_[[i]]$col_names()[2], '\n')
          cat("--------------",'\n')
          for (j in 1:nrow(private$item_map_[[i]]$values())) {
            cat(as.character(private$item_map_[[i]]$values()[j,1]), "=",
                as.character(private$item_map_[[i]]$values()[j,2]), '\n')
          }
        }
      } else {
        cat("==============",'\n')
        cat("empty mapping",'\n')
      }
      invisible(self)
    },

    #' @description Add new [SurveyQuestion]s.
    #' @param ... The [SurveyQuestion]s to add.
    add = function(...) {
      dots <- list(...)
      for (i in 1:length(dots)) {
        ll_length <- length(private$item_map_)
        if (dots[[i]]$name() %in% names(private$item_map_)) {
          stop("Survey label '", dots[[i]]$name(), "' already defined.  ",
               "Use 'replace' method instead.", call. = FALSE)
        }
        private$item_map_[[ll_length + 1]] <- dots[[i]]
        names(private$item_map_)[ll_length + 1] <- private$item_map_[[ll_length + 1]]$name()
      }
      invisible(self)
    },

    #' @description Delete [SurveyQuestion]s.
    #' @param ... The [SurveyQuestion]s to delete.
    delete = function(...) {
      tmp_list <- list(...)
      for (i in length(tmp_list)) {
        if (inherits(tmp_list[[i]], "SurveyQuestion")) {
          loc_id <- names(private$item_map_) %in% tmp_list[[i]]$name()
          loc_name <- tmp_list[[i]]$name()
        } else {
          loc_id <- names(private$item_map_) %in% tmp_list[[i]]
          loc_name <- tmp_list[[i]]
        }
        if (sum(loc_id) == 0) {
          stop("Survey label '", loc_name, "' not defined.",
               call. = FALSE)
        } else {
          private$item_map_[[which(loc_id)]] <- NULL
        }
      }
      invisible(self)
    },

    #' @description Replace one [SurveyQuestion] with another.
    #' @param old_question The [SurveyQuestion] object to replace.
    #' @param new_question The [SurveyQuestion] object to use instead.
    #'
    replace = function(old_question, new_question) {
      self$delete(old_question)
      self$add(new_question)
      invisible(self)
    },

    #' @description Validate the mapping
    validate = function() {
      samp_dfnames <- colnames(private$samp_obj_$survey_data(key = FALSE))
      popn_dfnames <- colnames(private$popn_obj_$survey_data(key = FALSE))
      samp_mapnames <- character(length(private$item_map_))
      popn_mapnames <- character(length(private$item_map_))
      for (j in 1:length(private$item_map_)) {
        samp_mapnames[j] <- private$item_map_[[j]]$col_names()[1]
        popn_mapnames[j] <- private$item_map_[[j]]$col_names()[2]
        if (!is.factor(private$samp_obj_$survey_data()[, samp_mapnames[j]])) {
          private$samp_obj_$add_survey_data_column(
            name = samp_mapnames[j],
            value = as.factor(private$samp_obj_$survey_data()[, samp_mapnames[j]])
          )
          warning("Converting '", samp_mapnames[j], "' into a factor with levels ",
                  paste(levels(private$samp_obj_$survey_data()[, samp_mapnames[j]]), collapse = ", "),
                  call. = FALSE)
        }
        if (!is.factor(private$popn_obj_$survey_data()[, popn_mapnames[j]])) {
          private$popn_obj_$add_survey_data_column(
            name = popn_mapnames[j],
            value = as.factor(private$popn_obj_$survey_data()[, popn_mapnames[j]])
          )
          warning("Converting '", popn_mapnames[j], "' into a factor with levels ",
                  paste(levels(private$popn_obj_$survey_data()[, popn_mapnames[j]]), collapse = ", "),
                  call. = FALSE)
        }
        levels_map_samp <- levels(private$item_map_[[j]]$values()[, 1])
        levels_map_popn <- levels(private$item_map_[[j]]$values()[, 2])
        levels_data_samp <- levels(private$samp_obj_$survey_data()[, samp_mapnames[j]])
        levels_data_popn <- levels(private$popn_obj_$survey_data()[, popn_mapnames[j]])
        if (!samp_mapnames[j] %in% samp_dfnames) {
          stop("Variable '", samp_mapnames[j], "' not in sample", call. = FALSE)
        }
        if (!popn_mapnames[j] %in% popn_dfnames) {
          stop("Variable '", popn_mapnames[j], "' not in population",
               call. = FALSE)
        }
        if (sum(!levels_map_samp %in% levels_data_samp) > 0) {
          stop("Levels in '",samp_mapnames[j],"' ",
               paste(levels_map_samp[!levels_map_samp %in% levels_data_samp], collapse = ", "),
               " are included in the map but are not in the sample",
               call. = FALSE)
        }
        if (sum(!levels_data_samp %in% levels_map_samp) > 0) {
          stop("Levels in '",samp_mapnames[j], "' ",
               paste(levels_data_samp[!levels_data_samp %in% levels_map_samp], collapse = ", "),
               " are included in the sample but are not in the map",
               call. = FALSE)
        }
        if (sum(!levels_map_popn %in% levels_data_popn) > 0) {
          stop("Levels in '",popn_mapnames[j], "' ",
               paste(levels_map_popn[!levels_map_popn %in% levels_data_popn], collapse = ", "),
               " are included in the map but are not in the population data",
               call. = FALSE)
        }
        if (sum(!levels_data_popn %in% levels_map_popn)>0) {
          stop("Levels in '",popn_mapnames[j], "' ",
               paste(levels_data_popn[!levels_data_popn %in% levels_map_popn], collapse = ", "),
               " are included in the population data but are not in the map",
               call. = FALSE)
        }
      }
      if (sum(popn_mapnames %in% popn_dfnames) < length(popn_dfnames)) {
        warning("Variable(s) ", paste(sQuote(popn_dfnames[!popn_dfnames %in% popn_mapnames], q = "'"), collapse = ", "),
                " are available in the population but won't be used in the model ",
                call. = FALSE)
      }
      if (sum(!samp_dfnames %in% c(samp_mapnames,popn_dfnames)) == 0) {
        warning("At least one variable in the survey needs to be unknown in the population.",
                call. = FALSE)
      }
      invisible(self)
    },

    #' @description Prepare the mapped data
    mapping  = function() {
      for (j in 1:length(private$item_map_)) {
        samp_mapnames <- private$item_map_[[j]]$col_names()[1]
        popn_mapnames <- private$item_map_[[j]]$col_names()[2]
        levels_map_samp <- private$item_map_[[j]]$values()[, 1]
        levels_map_popn <- private$item_map_[[j]]$values()[, 2]
        new_varname <- private$item_map_[[j]]$name()
        new_levels_samp <- character(length(levels_map_samp))
        new_levels_popn <- character(length(levels_map_popn))
        multiple_collapse <- rep(x = FALSE, length(levels_map_samp))
        for (k in 1:length(levels_map_samp)) {
          is_samp_unique <- sum(levels_map_samp %in% levels_map_samp[k]) == 1
          is_popn_unique <- sum(levels_map_popn %in% levels_map_popn[k]) == 1
          new_levels_samp[k] <- as.character(levels_map_samp[k])
          new_levels_popn[k] <- as.character(levels_map_popn[k])
          if (is_samp_unique && !is_popn_unique) {
            names(new_levels_samp)[k] <- as.character(levels_map_popn[k])
            names(new_levels_popn)[k] <- as.character(levels_map_popn[k])
          } else if (is_samp_unique & is_popn_unique) {
            names(new_levels_samp)[k] <- as.character(levels_map_samp[k])
            names(new_levels_popn)[k] <- as.character(levels_map_samp[k])
          } else if (!is_samp_unique && is_popn_unique) {
            names(new_levels_samp)[k] <- as.character(levels_map_samp[k])
            names(new_levels_popn)[k] <- as.character(levels_map_samp[k])
          } else if (!is_samp_unique & !is_popn_unique) {
            multiple_collapse[k] = TRUE
          }
        }
        if(sum(multiple_collapse)>0){

          dup_labels_popn <- levels_map_popn[duplicated(levels_map_popn)]
          dup_labels_samp <- levels_map_samp[duplicated(levels_map_samp)]

          dup_labels_samp <- levels_map_samp %in% dup_labels_samp
          dup_labels_popn <- levels_map_popn %in% dup_labels_popn

          potential_indices <- dup_labels_samp + dup_labels_popn>0
          if (sum(multiple_collapse[potential_indices]>0)){
            collapsed_names <- paste0(unique(levels_map_samp[potential_indices]), collapse = " + ")
            names(new_levels_popn)[potential_indices] <- collapsed_names
            names(new_levels_samp)[potential_indices] <- collapsed_names
          } else{
            stop("There is an error with the mapping. Please investigate further.")
          }

        }
        private$samp_obj_$add_mapped_data_column(new_varname, forcats::fct_recode(private$samp_obj_$survey_data()[[samp_mapnames]], !!!new_levels_samp))
        private$popn_obj_$add_mapped_data_column(new_varname, forcats::fct_recode(private$popn_obj_$survey_data()[[popn_mapnames]], !!!new_levels_popn))
      }
      invisible(self)
    },

    #' @description Prepare the poststratification table
    #' @param ... The variables to include.
    tabulate  = function(...) {
      grouping_vars <- c(...)
      if (length(grouping_vars) == 0) {
        grouping_vars <- names(private$item_map_)
      }
      if (sum(!grouping_vars %in% names(private$item_map_)) > 0) {
        stop("At least one poststratification variable doesn't correspond to the map.", call. = FALSE)
      }
      private$poststrat_data_ <-
        private$popn_obj_$mapped_data() %>%
        dplyr::mutate(wts = private$popn_obj_$weights()) %>%
        dplyr::group_by_at(dplyr::all_of(grouping_vars)) %>%
        dplyr::summarize(N_j = sum(wts), .groups = 'drop')
      invisible(self)
    },

    #' @description Fit a model.
    #' @param fun The model fitting function to use. For example,
    #'   `fun=rstanarm::stan_glmer`, `fun=brms::brm`, `fun=lme4::glmer`.
    #'   If using a custom `fun` it must have a `data` argument that accepts
    #'   a data frame (like standard R modeling functions).
    #' @param formula The model formula. Can be either a string or a formula
    #'   object.
    #' @param ... Arguments other than `formula` and `data` to pass to `fun`.
    #'   The data argument will be automatically specified internally.
    #' @return A [SurveyFit] object.
    #'
    fit = function(fun, formula, ...) {
      fun <- match.fun(fun)
      args <- list(...)
      if (!is.null(args$data)) {
        stop("The 'data' argument should not be specified.",
             call. = FALSE)
      }
      if ("family" %in% names(formals(fun)) &&
          !family_is_binomial(args$family)) {
        stop("Currently only binomial and bernoulli models are supported.",
             call. = FALSE)
      }
      if (identical(colnames(private$samp_obj_$mapped_data()), ".key")) {
        stop("Mapped data not found. ",
             "Please call the mapping() method before fitting a model.",
             call. = FALSE)
      }
      if (nrow(private$poststrat_data_) == 0) {
        stop("Post-stratification data not found. ",
             "Please call the tabulate() method before fitting a model.",
             call. = FALSE)
      }

      admin_package <- as.character(getNamespaceName(environment(fun)))
      if (!any(admin_package %in% c("lme4","brms","rstanarm"))){
        warning("Only rstanarm, brms and lme4 are supported natively. ",
                "Other modeling tools will need a custom population_predict() method.",
                call. = FALSE)
      }

      formula <- as.formula(formula)
      mapped_data <- private$samp_obj_$mapped_data()
      rhs_vars <- all.vars(formula[-2])
      lhs_vars <- all.vars(update(formula, "~0"))

      if (sum(!rhs_vars %in% colnames(mapped_data))) {
        stop("Not all variables available in the data. ",
             paste("Missing vars: ", paste(rhs_vars[!rhs_vars %in% colnames(mapped_data)], collapse = ', ')),
             call. = FALSE)
      }
      if (sum(!lhs_vars %in% colnames(private$samp_obj_$survey_data()))) {
        stop("Outcome variable not present in data. ",
             call. = FALSE)
      }
      if (sum(!rhs_vars %in% colnames(private$poststrat_data_))) {
        stop("Predictor variables not known in population. ",
             "Please ensure all predictor variables are mapped from sample to population. ",
             paste("Missing vars:", paste(rhs_vars[!rhs_vars %in% colnames(private$poststrat_data_)], collapse = ', ')),
             call. = FALSE)
      }

      need_vars <- setdiff(all.vars(formula), colnames(mapped_data))
      y_and_x <- private$samp_obj_$survey_data()[, need_vars, drop = FALSE]
      args$formula <- formula
      args$data <- cbind(mapped_data, y_and_x)
      fit <- do.call(fun, args)
      SurveyFit$new(fit = fit, map = self)
    },

    #' @description Access the `item_map`
    item_map = function() private$item_map_,

    #' @description Access the [SurveyData] object containing the sample data
    samp_obj = function() private$samp_obj_,

    #' @description Access the [SurveyData] object containing the population data
    popn_obj = function() private$popn_obj_,

    #' @description Access the poststratification data frame
    poststrat_data = function() {
      if (is.null(private$poststrat_data_)) {
        stop("Please call the tabulate() method before accessing the poststrat data.",
             call. = FALSE)
      }
      private$poststrat_data_
    }
  )
)

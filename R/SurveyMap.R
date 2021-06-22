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
#' head(shape_survey)
#' box_prefs <- SurveyData$new(
#'   data = shape_survey,
#'   questions = list(
#'     age = "Please identify your age group",
#'     gender = "Please select your gender",
#'     vote_for = "Which party did you vote for in the 2018 election?",
#'     y = "If today is the election day, would you vote for the Box Party?"
#'   ),
#'   responses = list(
#'     age = levels(shape_survey$age),
#'     gender = levels(shape_survey$gender),
#'     vote_for = levels(shape_survey$vote_for),
#'     y = c("no","yes")
#'   ),
#'   weights = "wt",
#'   design = list(ids =~1)
#' )
#' box_prefs$print()
#' box_prefs$n_questions()
#'
#' head(approx_voters_popn)
#' popn_obj <- SurveyData$new(
#'   data = approx_voters_popn,
#'   questions = list(
#'     age_group = "Which age group are you?",
#'     gender = "Gender?",
#'     vote_pref = "Which party do you prefer to vote for?"
#'   ),
#'   # order doesn't matter (gender before age2 here) because
#'   # the list has the names of the variables
#'   responses = list(
#'     gender = levels(approx_voters_popn$gender),
#'     age_group = levels(approx_voters_popn$age_group),
#'     vote_pref = levels(approx_voters_popn$vote_pref)
#'   ),
#'   weights = "wt"
#' )
#' popn_obj$print()
#'
#' q_age <- QuestionMap$new(
#'   name = "age",
#'   col_names = c("age","age_group"),
#'   values_map = list(
#'     "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
#'     "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
#'   )
#' )
#' print(q_age)
#'
#' q_party_pref <- QuestionMap$new(
#'   name = "party_pref",
#'   col_names = c("vote_for","vote_pref"),
#'   values_map = list("Box Party" = "BP",  "BP" = "BP","Circle Party" = "CP", "CP" = "CP")
#' )
#' q_gender <- QuestionMap$new(
#'   name = "gender",
#'   col_names = c("gender", "gender"),
#'   values_map = list("male" = "m","female" = "f", "nonbinary" = "nb")
#' )
#'
#' # Create SurveyMap object adding all questions at once
#' ex_map <- SurveyMap$new(
#'   sample = box_prefs,
#'   population = popn_obj,
#'   q_age,
#'   q_party_pref,
#'   q_gender
#' )
#' print(ex_map) # or ex_map$print()
#'
#' # Or can add questions incrementally
#' ex_map <- SurveyMap$new(sample = box_prefs, population = popn_obj)
#' print(ex_map)
#'
#' ex_map$add(q_age, q_party_pref)
#' print(ex_map)
#'
#' ex_map$add(q_gender)
#' print(ex_map)
#'
#' # Create the mapping between sample and population
#' ex_map$mapping()
#'
#' # Create the poststratification data frame using all variables in the mapping
#' # (alternatively, can specify particular variables, e.g. tabulate("age"))
#' ex_map$tabulate()
#'
#' # Example rstanarm usage
#' # Returns a SurveyFit object
#' fit_1 <- ex_map$fit(
#'   fun = rstanarm::stan_glmer,
#'   formula = y ~ (1|age) + (1|gender),
#'   family = "binomial",
#'   seed = 1111,
#'   # just to keep the example fast and small
#'   chains = 1
#' )
#'
#' # Example lme4 usage
#' # fit_2 <- ex_map$fit(
#' #   fun = lme4::glmer,
#' #   formula = y ~ (1|age) + (1|gender),
#' #   family = "binomial"
#' # )
#' #
#' # Example brms usage
#' # fit_3 <- ex_map$fit(
#' #   fun = brms::brm,
#' #   formula = y ~ (1|age) + (1|gender),
#' #   family = "bernoulli",
#' #   seed = 1111
#' # )
#'
#'
#' # predicted probabilities
#' # returns matrix with rows for poststrat cells, cols for posterior draws
#' poststrat_estimates <- fit_1$population_predict()
#'
#' # estimates by age level
#' estimates_by_age <- fit_1$aggregate(poststrat_estimates, by = "age")
#' head(estimates_by_age)
#' estimates_by_age %>%
#'   group_by(age) %>%
#'   summarize(mean = mean(value), sd = sd(value))
#'
#' # plot estimates by age
#' fit_1$plot(estimates_by_age, weights = FALSE)
#' fit_1$plot(estimates_by_age, weights = TRUE)
#'
#' # population estimate
#' estimates_popn <- fit_1$aggregate(poststrat_estimates)
#' mean(estimates_popn$value)
#'
#' # plot population estimate
#' fit_1$plot(estimates_popn, weights = FALSE)
#' fit_1$plot(estimates_popn, weights = TRUE)
#'
SurveyMap <- R6::R6Class(
  classname = "SurveyMap",
  private = list(
    sample_ = NULL,
    population_ = NULL,
    item_map_ = list(),
    poststrat_data_ = data.frame(NULL),
    mapped_sample_data_ = NULL,
    mapped_population_data_ = NULL
  ),
  public = list(

    #' @description Create a new SurveyMap object
    #' @param sample The [SurveyData] object corresponding to the sample data.
    #' @param population The [SurveyData] object corresponding to the population data.
    #' @param ... [QuestionMap] objects.
    initialize = function(sample, population, ...) {
      if (!inherits(sample, "SurveyData")) {
        stop("'sample' must be a SurveyData object.", call. = FALSE)
      }
      if (!inherits(population, "SurveyData")) {
        stop("'population' must be a SurveyData object.", call. = FALSE)
      }

      private$item_map_ <- list(...)
      for (i in seq_along(private$item_map_)) {
        names(private$item_map_)[i] <- private$item_map_[[i]]$name()
      }

      private$sample_ <- sample
      private$population_ <- population
      private$mapped_sample_data_ <- data.frame(.key = 1:nrow(sample$survey_data()))
      private$mapped_population_data_ <- data.frame(.key = 1:nrow(population$survey_data()))
      self$validate()
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

    #' @description Add new [QuestionMap]s.
    #' @param ... The [QuestionMap]s to add.
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
      self$validate()
      invisible(self)
    },

    #' @description Delete [QuestionMap]s.
    #' @param ... The [QuestionMap]s to delete.
    delete = function(...) {
      tmp_list <- list(...)
      for (i in length(tmp_list)) {
        if (inherits(tmp_list[[i]], "QuestionMap")) {
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
      self$validate()
      invisible(self)
    },

    #' @description Replace one [QuestionMap] with another.
    #' @param old_question The [QuestionMap] object to replace.
    #' @param new_question The [QuestionMap] object to use instead.
    #'
    replace = function(old_question, new_question) {
      self$delete(old_question)
      self$add(new_question)
      invisible(self)
    },

    #' @description Validate the mapping
    validate = function() {
      if (length(private$item_map_) == 0) {
        return(invisible(self))
      }
      samp_dfnames <- colnames(private$sample_$survey_data(key = FALSE))
      popn_dfnames <- colnames(private$population_$survey_data(key = FALSE))
      samp_mapnames <- character(length(private$item_map_))
      popn_mapnames <- character(length(private$item_map_))
      for (j in 1:length(private$item_map_)) {
        samp_mapnames[j] <- private$item_map_[[j]]$col_names()[1]
        popn_mapnames[j] <- private$item_map_[[j]]$col_names()[2]
        if (!is.factor(private$sample_$survey_data()[, samp_mapnames[j]])) {
          private$sample_$add_survey_data_column(
            name = samp_mapnames[j],
            value = as.factor(private$sample_$survey_data()[, samp_mapnames[j]])
          )
          warning("Converting '", samp_mapnames[j], "' into a factor with levels ",
                  paste(levels(private$sample_$survey_data()[, samp_mapnames[j]]), collapse = ", "),
                  call. = FALSE)
        }
        if (!is.factor(private$population_$survey_data()[, popn_mapnames[j]])) {
          private$population_$add_survey_data_column(
            name = popn_mapnames[j],
            value = as.factor(private$population_$survey_data()[, popn_mapnames[j]])
          )
          warning("Converting '", popn_mapnames[j], "' into a factor with levels ",
                  paste(levels(private$population_$survey_data()[, popn_mapnames[j]]), collapse = ", "),
                  call. = FALSE)
        }
        levels_map_samp <- levels(private$item_map_[[j]]$values()[, 1])
        levels_map_popn <- levels(private$item_map_[[j]]$values()[, 2])
        levels_data_samp <- levels(private$sample_$survey_data()[, samp_mapnames[j]])
        levels_data_popn <- levels(private$population_$survey_data()[, popn_mapnames[j]])
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
      if (anyNA(private$sample_$survey_data()[, samp_mapnames]) ||
          anyNA(private$population_$survey_data()[, popn_mapnames])) {
        stop("NAs not allowed in variables mapped between sample and population.", call. = FALSE)
      }

      invisible(self)
    },

    #' @description The mapping method uses the given maps between questions to create new sample and population data
    #' that has unified variable names (i.e., if the underlying construct is called `age`, both sample and population will have
    #' an `age` column, even if in the the raw data both had different variable names). The method also unifies the levels of each
    #' variable in the sample and population so that the maximum set of consistent levels is created. Names of these new levels are
    #' made according the the sample level names. If multiple levels are combined, the new name will be existing levels seperated by a
    #' ` + `.
    #'
    mapping  = function() {
      for (j in 1:length(private$item_map_)) {
        samp_mapnames <- private$item_map_[[j]]$col_names()[1]
        popn_mapnames <- private$item_map_[[j]]$col_names()[2]
        levels_samp <- levels(private$item_map_[[j]]$values()[, 1])
        levels_popn <- levels(private$item_map_[[j]]$values()[, 2])
        levels_map_samp <- private$item_map_[[j]]$values()[, 1]
        levels_map_popn <- private$item_map_[[j]]$values()[, 2]
        new_varname <- private$item_map_[[j]]$name()
        new_levels_samp <- character(length(levels_map_samp))
        new_levels_popn <- character(length(levels_map_popn))

        # Find the naming clusters of levels  (loosely inspired by complete linkage clustering algorithms.
        # Three major steps:
        # Step 1:  Create a matrix that has the factor level labels for the sample as row names, and factor level labels for the
        # population as the column names. If the levels are linked, record 1 in the corresponding entry, otherwise record 0.
        mapped_levels <- matrix(0,nrow=length(unique(levels_map_samp)),ncol = length(unique(levels_map_popn)))
        colnames(mapped_levels)<-levels_popn
        row.names(mapped_levels)<-levels_samp
        for (unique_samp_levels in 1:length(mapped_levels[,1])) {
          which_samp_levels <- levels_map_samp == levels_samp[unique_samp_levels]
          for(corres_popn in levels_map_popn[which_samp_levels]){
            mapped_levels[unique_samp_levels,unique(levels_popn) == corres_popn] <-1
          }
        }
        if(sum(rowSums(mapped_levels, na.rm=TRUE)!=0)<nrow(mapped_levels)){
          stop("Levels: ",paste(row.names(sum(rowSums(mapped_levels, na.rm=TRUE)!=0)),collapse = " ")," do not have a match in the population.")
        }
        if(sum(colSums(mapped_levels, na.rm=TRUE)!=0)<ncol(mapped_levels)){
          stop("Levels: ",paste(colnames(sum(colSums(mapped_levels, na.rm=TRUE)!=0)),collapse = " ")," do not have a match in the sample.")
        }

        # Moving down the rows (i.e., for each factor level in the sample), identify if there are multiple corresponding population
        # values. If there are, merge the columns for the values together and combine the col names together, seperated by a ' + ' symbol.
        tmp_mapped_levels <- mapped_levels
        mapped_levels_new <- mapped_levels
        for (unique_samp_levels in 1:length(mapped_levels[,1])) {
          if(sum(mapped_levels[unique_samp_levels,],na.rm = TRUE)>1){
            mapped_levels_new <- tmp_mapped_levels[,tmp_mapped_levels[unique_samp_levels,]!=1, drop = FALSE]
            ll = ncol(mapped_levels_new)
            mapped_levels_new <- cbind(mapped_levels_new,ifelse(rowSums(tmp_mapped_levels[,tmp_mapped_levels[unique_samp_levels,]==1])>=1,1,0))
            colnames(mapped_levels_new)[ll+1] <- paste0(colnames(tmp_mapped_levels[,tmp_mapped_levels[unique_samp_levels,,drop = FALSE]==1,drop = FALSE]), collapse = " + ")
            tmp_mapped_levels<- mapped_levels_new        }
        }

        # Moving across the now abridged columns (i.e., for every new facter level in the popn), identify if there are multiple
        # corresponding sample values. If there are, merge the rows for the values together and combine the row names (sample factor levels)
        # together, seperated by a " + " symbol

        tmp_mapped_levels <- mapped_levels_new
        mapped_levels_fin <- mapped_levels_new
        for (unique_popn_levels in 1:length(mapped_levels_fin[1,])) {
          if(sum(mapped_levels_new[,unique_popn_levels],na.rm = TRUE)>1){
            mapped_levels_new <- tmp_mapped_levels[tmp_mapped_levels[,unique_popn_levels]!=1,, drop = FALSE]
            ll = nrow(mapped_levels_new)
            mapped_levels_new <- rbind(mapped_levels_new,ifelse(colSums(mapped_levels_fin[mapped_levels_fin[,unique_popn_levels]==1,])>=1,1,0))
            row.names(mapped_levels_new)[ll+1] <- paste0(row.names(tmp_mapped_levels[tmp_mapped_levels[,unique_popn_levels, drop = FALSE]==1,,drop = FALSE]), collapse = " + ")
            tmp_mapped_levels<- mapped_levels_new
          }
        }
        mapped_levels_fin <- mapped_levels_new

        # Create the name remapping, with the old names as the previous factor levels, and the names of the vector as the new levels
        # For the sample
        new_levels_samp_names <- length(new_levels_samp)
        for(samp_level in 1:length(levels_samp)){
          if(length(grep(levels_samp[samp_level],row.names(mapped_levels_fin)))>1){
            #specific match
            new_levels_samp_names[samp_level] <- row.names(mapped_levels_fin)[grep(paste0("^",levels_samp[samp_level],"$"),row.names(mapped_levels_fin))]
          }else{
            #partial match
            new_levels_samp_names[samp_level] <- row.names(mapped_levels_fin)[grep(levels_samp[samp_level],row.names(mapped_levels_fin))]
          }
        }
        new_levels_samp <- levels_samp
        names(new_levels_samp) <- new_levels_samp_names
        # For the population
        new_levels_popn_names <- length(new_levels_popn)
        for(popn_level in 1:length(levels_popn)){
          if(length(grep(paste0(levels_popn[popn_level]),colnames(mapped_levels_fin)))>1){
            popn_level_loc <- grep(paste0("^",levels_popn[popn_level],"$"),colnames(mapped_levels_fin))
          } else{
            popn_level_loc <- grep(levels_popn[popn_level],colnames(mapped_levels_fin))
          }
          #name the population levels according to the sample data
          new_levels_popn_names[popn_level] <- row.names(mapped_levels_fin)[mapped_levels_fin[,popn_level_loc]==1]
        }
        new_levels_popn <- levels_popn
        names(new_levels_popn) <- new_levels_popn_names

        # Finally rename the data
        private$mapped_sample_data_[[new_varname]] <- forcats::fct_recode(private$sample_$survey_data()[[samp_mapnames]], !!!new_levels_samp)
        private$mapped_population_data_[[new_varname]] <- forcats::fct_recode(private$population_$survey_data()[[popn_mapnames]], !!!new_levels_popn)
      }
    invisible(self)
    },

    #' @description Prepare the poststratification table
    #' @param ... The variables to include.
    tabulate = function(...) {
      if (ncol(self$mapped_population_data(key = FALSE)) == 0) {
        stop("Please call the mapping() method before tabulate()", call. = FALSE)
      }
      grouping_vars <- c(...)
      if (length(grouping_vars) == 0) {
        grouping_vars <- names(private$item_map_)
      }
      if (sum(!grouping_vars %in% names(private$item_map_)) > 0) {
        stop("At least one poststratification variable doesn't correspond to the map.", call. = FALSE)
      }
      private$poststrat_data_ <-
        private$mapped_population_data_ %>%
        dplyr::mutate(wts = private$population_$weights()) %>%
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
      if (ncol(self$mapped_population_data(key = FALSE)) == 0) {
        stop("Mapped data not found. ",
             "Please call the mapping() method before fitting a model.",
             call. = FALSE)
      }
      if (nrow(private$poststrat_data_) == 0) {
        stop("Post-stratification data not found. ",
             "Please call the tabulate() method before fitting a model.",
             call. = FALSE)
      }

      formula <- as.formula(formula)
      mapped_data <- private$mapped_sample_data_
      rhs_vars <- all.vars(formula[-2])
      lhs_vars <- all.vars(update(formula, "~0"))

      if (sum(!rhs_vars %in% colnames(mapped_data))) {
        stop("Not all variables available in the data. ",
             paste("Missing vars: ", paste(rhs_vars[!rhs_vars %in% colnames(mapped_data)], collapse = ', ')),
             call. = FALSE)
      }
      if (sum(!lhs_vars %in% colnames(private$sample_$survey_data()))) {
        stop("Outcome variable not present in data. ",
             call. = FALSE)
      }
      if (sum(!rhs_vars %in% colnames(private$poststrat_data_))) {
        stop("Predictor variables not known in population. ",
             "Please ensure all predictor variables are mapped from sample to population. ",
             paste("Missing vars:", paste(rhs_vars[!rhs_vars %in% colnames(private$poststrat_data_)], collapse = ', ')),
             call. = FALSE)
      }
      if (anyNA(private$sample_$survey_data()[[lhs_vars]])) {
        warning("Outcome variable has missing values that may be dropped ",
                "by the model fitting package.",
                call. = FALSE
        )
      }
      if (!any(getNamespaceName(environment(fun)) %in% c("lme4","brms","rstanarm"))) {
        warning("Only rstanarm, brms and lme4 are supported natively. ",
                "Other modeling tools will need a custom population_predict() method.",
                call. = FALSE)
      }

      need_vars <- setdiff(all.vars(formula), colnames(mapped_data))
      y_and_x <- private$sample_$survey_data()[, need_vars, drop = FALSE]
      args$formula <- formula
      args$data <- cbind(mapped_data, y_and_x)
      fit <- do.call(fun, args)
      SurveyFit$new(fit = fit, map = self)
    },

    #' @description Access the `item_map`
    item_map = function() private$item_map_,

    #' @description Access the [SurveyData] object containing the sample data
    sample = function() private$sample_,

    #' @description Access the [SurveyData] object containing the population data
    population = function() private$population_,

    #' @description Access the poststratification data frame created by the `tabulate` method
    poststrat_data = function() {
      if (is.null(private$poststrat_data_)) {
        stop("Please call the tabulate() method before accessing the poststrat data.",
             call. = FALSE)
      }
      private$poststrat_data_
    },

    #' @description Access the data frame containing the mapped sample data
    #'   created by the `mapping` method
    #' @param key Should the `.key` column be included? This column just
    #'   indicates the original order of the rows and is primarily intended
    #'   for internal use.
    mapped_sample_data = function(key = TRUE) {
      if (key) {
        private$mapped_sample_data_
      } else {
        private$mapped_sample_data_[, colnames(private$mapped_sample_data_) != ".key", drop = FALSE]
      }
    },

    #' @description Access the data frame containing the mapped population data
    #'   created by the `mapping` method
    #' @param key Should the `.key` column be included? This column just
    #'   indicates the original order of the rows and is primarily intended
    #'   for internal use.
    mapped_population_data = function(key = TRUE) {
      if (key) {
        private$mapped_population_data_
      } else {
        private$mapped_population_data_[, colnames(private$mapped_population_data_) != ".key", drop = FALSE]
      }
    }
  )
)

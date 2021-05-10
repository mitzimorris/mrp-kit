#' SurveyData objects
#'
#' @name SurveyData
#' @export
#' @description
#' A `SurveyData` object represents a survey and its metadata. The survey itself
#' is a data frame. The survey metadata consists of the text of the survey
#' questions, the allowed response values, and optionally survey weights and a
#' survey design specification.
#'
#' @examples
#'
#' head(nlp_survey)
#' nlp_prefs <- SurveyData$new(
#'   data = nlp_prefs,
#'   questions = list(
#'     age = "Please identify your age group",
#'     gender = "Please select your gender",
#'     vote_for = "Which party did you vote for in the 2018 election?",
#'     highest_educ = "Please identify your completed highest education",
#'     state = "Please identify the state do you live in"
#'     y = "Will you vote for the NLP if today is the election day?"
#'   ),
#'   responses = list(
#'     age = levels(nlp_survey$age),
#'     gender = levels(nlp_survey$gender),
#'     vote_for = levels(nlp_survey$vote_for),
#'     highest_educ = levels(nlp_survey$highest_educ),
#'     state = levels(nlp_survey$state),
#'     y = c("no","yes")
#'   ),
#'   weights = nlp_survey$wt,
#'   design = list(ids =~1)
#' )
#' nlp_prefs$print()
#'
#' head(approx_voters_popn)
#' popn_obj <- SurveyData$new(
#'   data = approx_voters_popn,
#'   questions = list(
#'     age_group = "Which age group are you?",
#'     gender = "Gender?",
#'     vote_pref = "Which party do you prefer to vote for?",
#'     education = "What is your highest grade completed?",
#'     state = "Which state do you live in?"
#'   ),
#'   # order doesn't matter (gender before age2 here) because
#'   # the list has the names of the variables
#'   responses = list(
#'     gender = levels(approx_voters_popn$gender),
#'     age_group = levels(approx_voters_popn$age_group),
#'     vote_pref = levels(approx_voters_popn$vote_pref),
#'     education = levels(approx_voters_popn$education),
#'     state = levels(approx_voters_popn$state)
#'   ),
#'   weights = approx_voters_popn$wt
#' )
#' popn_obj$print()
#'
SurveyData <- R6::R6Class(
    classname = "SurveyData",
    private = list(
        survey_data_ = data.frame(NULL),
        mapped_data_ = data.frame(NULL),
        questions_ = list(),
        responses_ = list(),
        weights_ = numeric(),
        design_ = list()
    ),
    public = list(
        #' @description Create a new SurveyData object
        #' @param data A data frame containing the survey data.
        #' @param questions,responses Named lists containing the text of the
        #'   survey questions and the allowed responses, respectively. The names
        #'   must correspond to the names of variables in `data`. See
        #'   **Examples**.
        #' @param weights Optionally, a vector of survey weights.
        #' @param design Optionally, a named list of arguments to pass to `survey::svydesign()`.
        initialize = function(data,
                              questions = list(),
                              responses = list(),
                              weights = numeric(),
                              design = list(ids =~1)) {
            if (ncol(data) == 0 || nrow(data) == 0) {
                stop("'data' cannot be empty.", call. = FALSE)
            }
            # allow no question/response, else require info for all columns
            if (length(questions) != 0 || length(responses) != 0) {
                if (ncol(data) != length(questions) &
                    length(questions) == sum(stats::complete.cases(data))) {
                    stop("Mismatch between number of data columns and questions.",
                         call. = FALSE)
                }
                if (length(responses) != length(questions)) {
                    stop("Mismatch between number of survey questions and answers.",
                         call. = FALSE)
                }
            }
            # allow no weights, else require weights for all rows
            if (length(weights) != 0 && nrow(data) != length(weights)) {
                stop("Mismatch between number of data columns and weights.",
                     call. = FALSE)
            }

            nms_q <- sort(names(questions))
            nms_r <- sort(names(responses))
            if (is.null(nms_q) || sum(nzchar(nms_q)) != length(nms_q)) {
                stop("All elements of 'questions' and 'responses' list must have names.", call. = FALSE)
            }
            if (!identical(nms_q, nms_r)) {
                stop("Names in 'questions' and 'responses' lists must be the same.")
            }
            if (!all(nms_q %in% colnames(data))) {
                stop("Names of 'questions' much match column names in 'data'.", call. = FALSE)
            }
            questions <- questions[nms_q]
            responses <- responses[nms_q]

            private$questions_ <- questions
            private$responses_ <- responses
            private$weights_ <- weights
            private$design_ <- design
            private$survey_data_ <- data.frame(.key = 1:nrow(data), data)
            private$mapped_data_ <- data.frame(.key = 1:nrow(data))
            invisible(self)
        },

        #' @description Number of observations in the survey data
        n_obs = function() nrow(private$survey_data_),

        #' @description Number of survey questions
        n_questions = function() length(private$questions_),

        #' @description Print a summary of the survey data
        #' @param ... Currently ignored.
        print = function(...) {
            cat("Survey with",
                self$n_obs(), "observations,",
                self$n_questions(), "questions",
                "\n")

            print_survey_design(private$design_, private$weights_, private$survey_data_)
            print_questions_and_responses(private$questions_, private$responses_)
            invisible(self)
        },

        #' @description Add a column to the sample data. This is primarily
        #'   intended for internal use but may occasionally be useful.
        #' @param name,value The name of the new variable (a string) and the
        #' vector of values to add to the data frame.
        add_survey_data_column = function(name, value) {
            private$survey_data_[[name]] <- value
            invisible(self)
        },
        #' @description Add a column to the mapped data. This is primarily
        #'   intended for internal use but may occasionally be useful.
        #' @param name,value The name of the new variable (a string) and the
        #' vector of values to add to the data frame.
        add_mapped_data_column = function(name, value) {
            if (ncol(private$mapped_data_)  == 0) {
                private$mapped_data_ <- data.frame(value)
                colnames(private$mapped_data_) <- name
            } else {
                private$mapped_data_[[name]] <- value
            }
            invisible(self)
        },

        #' @description Access the data frame containing the sample data.
        #' @param key Should the `.key` column be included? This column just
        #'   indicates the original order of the rows and is primarily intended
        #'   for internal use.
        survey_data = function(key = TRUE) {
            if (key) {
                private$survey_data_
            } else {
                private$survey_data_[, colnames(private$survey_data_) != ".key", drop = FALSE]
            }
        },

        #' @description Access the data frame containing the mapped data.
        #' @param key Should the `.key` column be included? This column just
        #'   indicates the original order of the rows and is primarily intended
        #'   for internal use.
        mapped_data = function(key = TRUE) {
            if (key) {
                private$mapped_data_
            } else {
                private$mapped_data_[, colnames(private$mapped_data_) != ".key", drop = FALSE]
            }
        },

        #' @description Access the list of survey questions
        questions = function() private$questions_,
        #' @description Access the list of allowed survey responses
        responses = function() private$responses_,
        #' @description Access the survey weights
        weights = function() private$weights_,
        #' @description Access the survey design
        design = function() private$design_
    )
)


# internal ----------------------------------------------------------------

# print 1-line summary of survey design
print_survey_design <- function(design, weights, data) {
  svy_design <- do.call(survey::svydesign, c(design, list(weights = weights, data = data)))
  svy_design$call <- NULL
  cat(utils::capture.output(print(svy_design))[1], "\n")
}

print_questions_and_responses <- function(questions, responses) {
  for (i in seq_along(questions)) {
    cat("\nColumn label:", names(questions)[i], "\n")
    cat("Question:", questions[[i]], "\n")
    cat("Allowed answers:", paste(responses[[i]], collapse = ", "), "\n")
  }
}

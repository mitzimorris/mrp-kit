#' SurveyQuestion
#'
#' @name SurveyQuestion
#' @export
#'
#' @description An [R6][R6::R6Class] `SurveyQuestion` object holds the mapping
#'   for one question or demographic between the survey and population dataset.
#'
#' @examples
#' q1 <- SurveyQuestion$new(
#'   name = "age",
#'   col_names = c("age","age_group"),
#'   values_map = list(
#'     "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
#'     "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
#'   )
#' )
#' print(q1)
#'
#' q2 <- SurveyQuestion$new(
#'   name = "party_pref",
#'   col_names = c("vote_for","vote_pref"),
#'   values_map = list("Neverland Labor Party" = "NLP", "NLP" = "NLP","Neverland Democrats" = "The Democrats","The Democrats" = "The Democrats")
#' )
#' q3 <- SurveyQuestion$new(
#'   name = "gender",
#'   col_names = c("gender", "gender"),
#'   values_map = list("male" = "m","female" = "f", "nonbinary" = "nb")
#' )
#' q4 <- SurveyQuestion$new(
#'   name = "highest_education",
#'   col_names = c("highest_educ", "education"),
#'   values_map = list("no high school" = "no high school",
#'                     "high school" = "high school",
#'                     "some college" = "some college",
#'                     "4-year college" = "4-year college",
#'                     "post-graduate" = "post-grad")
#' )
#' q5 <- SurveyQuestion$new(
#'   name = "state",
#'   col_names = c("state", "state"),
#'   values_map = list("State A" = "A", "State B" = "B", "State C" = "C",
#'                     "State D" = "D", "State E" = "E")
#' )
SurveyQuestion <- R6::R6Class(
  classname = "SurveyQuestion",
  private = list(
    name_ = character(),
    col_names_ = character(),
    values_ = list()
  ),
  public = list(
    #' @description Create a new `SurveyQuestion`
    #' @param name The name of the underlying construct. For example if the
    #'   sample data uses the name `age1` and the population data uses
    #'   the name `age2` then `name` could be `age`.
    #' @param col_names A character vector of length two giving the column names
    #'   of the question in the sample data (first element) and the population
    #'   data (second element).
    #' @param values_map A named list where the names correspond to the
    #'   responses to the question in the sample data and the values correspond
    #'   to the responses in the population data. If there is a meaningful
    #'   ordering over the values, they should be listed in that order, either
    #'   descending or ascending. See **Examples**.
    #' @return A `SurveyQuestion` object that can be added to a [`SurveyMap`]
    #'   object.
    #'
    initialize = function(name, col_names, values_map) {
      if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single string.", call. = FALSE)
      }
      if (is.na(name)) {
        stop("'name' cannot be NA.", call. = FALSE)
      }
      if (!is.character(col_names) || length(col_names) != 2) {
        stop("'col_names' must be a character vector of length 2.", call. = FALSE)
      }
      if (anyNA(col_names)) {
        stop("NAs not allowed in 'col_names'.", call. = FALSE)
      }
      if (!is.list(values_map)) {
        stop("'values_map' must be a list", call. = FALSE)
      }
      if (!all(nzchar(names(values_map)))) {
        stop("All elements of 'values_map' must have names.", call. = FALSE)
      }
      if (anyNA(values_map)) {
        stop("NAs not allowed in 'values_map'", call. = FALSE)
      }

      private$name_ <- name
      private$col_names_ <- col_names
      private$values_ <- data.frame(names(values_map), as.character(unlist(values_map)),
                                    stringsAsFactors = TRUE)
      if (sum(duplicated(private$values_)) > 0) {
        warning("Duplicated values in map, removing duplicates", call. = FALSE)
        private$values_ <- private$values_[!duplicated(private$values_), , drop = FALSE]
      }
      names(private$values_) <- col_names
      invisible(self)
    },

    #' @description Print a `SurveyQuestion`.
    #' @param ... Currently ignored.
    #' @return The `SurveyQuestion` object, invisibly.
    print = function(...) {
      cat("--------------",'\n')
      cat(private$name_,'\n')
      cat(private$col_names_[1], "=", private$col_names_[2], '\n')
      cat("--------------",'\n')
      for (j in 1:nrow(private$values_)) {
        cat(as.character(private$values_[j, 1]), "=",
            as.character(private$values_[j, 2]), '\n')
      }
      invisible(self)
    },

    #' @description Access `name`
    name = function() private$name_,

    #' @description Access `col_names`
    col_names = function() private$col_names_,

    #' @description Access `values`
    values = function() private$values_
  )
)

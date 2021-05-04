#' Example survey for package testing
#' @source Created in CreateTestData.R
#' @format A data frame with columns:
#' \describe{
#'  \item{age1}{Age group}
#'  \item{gender}{Gender}
#'  \item{pet_own}{Pet ownership}
#'  \item{y}{Example binary outcome variable}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#' str(feline_survey)
#'
"feline_survey"

#' Example approximate population (modeled after ACS)
#' @source Created in CreateTestData.R
#' @format A data frame with columns:
#' \describe{
#'  \item{age2}{Age group}
#'  \item{gender}{Gender}
#'  \item{pet_pref}{Pet preference}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#'
#'  approx_popn
#'
"approx_popn"


#' A Simulated Data of NLP Voters Survey
#'
#' This is a simulated survey of voter preferences for the Neverland Labor Party
#' (NLP) that also includes several demographic variables. We use this data
#' throughout the documentation to demonstrate the usage of the package.
#'
#' @source Created in NeverlandVoters.R
#'
#' @format A data frame with columns:
#'
#' \describe{
#'  \item{age}{Age group}
#'  \item{gender}{Gender}
#'  \item{vote_for}{Party voted for in the 2016 election}
#'  \item{highest_educ}{Highest grade of education completed}
#'  \item{state}{State}
#'  \item{y}{Preference for NLP in 2020 election}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#'
#'  nlp_survey
#'
"nlp_survey"

#' Simulated approximate population data for the NLP Voters survey (modeled after ACS)
#'
#' @source Created in NeverlandVoters.R
#'
#' @format A data frame with columns:
#'
#' \describe{
#'  \item{age_group}{Age group}
#'  \item{gender}{Gender}
#'  \item{vote_for}{Party voted for in the 2016 election}
#'  \item{education}{Highest grade of education completed}
#'  \item{state}{State}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#'
#'  approx_voters_popn
#'
"approx_voters_popn"


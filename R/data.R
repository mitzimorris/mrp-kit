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

#' Example approximate population (modelled after ACS)
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
#' This is a generated data shows the survey result of the preference
#' to Neverland Labor Party (NLP)
#' along with several demographic variables to demonstrate the usage of the package.
#'
#' @source Created in NeverlandVoters.R
#'
#' @format A data frame with columns:
#'
#' \describe{
#'  \item{age}{Age group}
#'  \item{gender}{Gender}
#'  \item{vote_for}{The party that the respondent vote for in the 2016 election}
#'  \item{highest_educ}{Highest grade of education that is completed}
#'  \item{state}{State}
#'  \item{y}{Preference to NLP in 2020 election}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#'
#'  nlp_survey
#'
"nlp_survey"

#' A Simulated Data of NLP Voters Approximate Population (Modelled after ACS)
#'
#' This is a generated data shows the approximate population of the
#' Neverland Labor Party (NLP) Survey to demonstrate the usage of the package.
#'
#' @source Created in NeverlandVoters.R
#'
#' @format A data frame with columns:
#'
#' \describe{
#'  \item{age_group}{Age group}
#'  \item{gender}{Gender}
#'  \item{vote_for}{The party that being voted in the 2016 election}
#'  \item{education}{Highest grade of education that is completed}
#'  \item{state}{State}
#'  \item{wt}{Raked and trimmed survey weight}
#' }
#' @examples
#'
#'  approx_voters_popn
#'
"approx_voters_popn"


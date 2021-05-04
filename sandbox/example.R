
# Read in library
library(mrpkit)


#### Set up workspace ####
library(tidyverse)

# The user needs to have written a CSV that aligns the two survey and poststat
# Here we just read them in
aligned_levels <- "aligned.csv"
# If we have misaligned - 18-38 then need to go to make the one that is smallest.
# TODO: This could be a helper function? For age-groups: put in one set of break-points and another, then it returns the best options. For race - could map?
aligned_levels <- readr::read_csv(aligned_levels)

# Make simulated ages and gender responses
survey_data <- tibble(age =
                        sample(c("<18", "18-24", "25-34", "35-44", "45-55", "55-64", "65-74", "75-84", "85-94", "94+"),
                               100,
                               replace = TRUE),
                      gender =
                        sample(c("male", "female", "other", "baby"),
                               100,
                               replace = TRUE),
                      dog_owner =
                        sample(c("yes", "no", "cat"),
                               100,
                               replace = TRUE)
                      )
# Then save it as a CSV
survey_data$votes <- sample(c("R", "D"), 100, replace = TRUE)
write_csv(survey_data, "survey_data.csv")




#### Actual example ####
# Survey
# Use read_and_prepare_survey(). Just look at response and two columns.
my_survey_data <-
  mrpkit::read_and_prepare_survey("survey_data.csv", c("votes", "age", "gender"), c("character", "character", "character"))

# Pretend that we got the class wrong
mrpkit::read_and_prepare_survey("survey_data.csv", c("votes", "age", "gender"), c("integer", "character", "character"))

# TODO: Replace the function with a name prepare_survey, where the CSV input is a bit more general - can take a CSV or a dataframe.

# Poststat
# ATM this is just the same
# TODO: Maybe this is prepare_population? Toward our 'grammar'. If the ACS, then it's weights.


# Mash them
# Survey
# Use read_and_prepare_survey(). Just look at response and two columns.
mrpkit::bring_together_survey_and_poststrat(survey_data, aligned_levels, c("age", "gender"))
# TODO: Rename this to bring_together_survey_and_population. Then the next function is make_post_strat()



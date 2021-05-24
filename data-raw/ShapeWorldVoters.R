library(DeclareDesign)
library(brms)
library(survey)
library(dplyr)
library(forcats)


#declare the design of the population

design <- declare_population(N = 50000,
                             age_group = sample.int(7, size=N, replace = TRUE, prob = c(1,1,2,2,3,2,1)),
                             gender = sample.int(3, size=N, replace = TRUE, prob = c(.48,.49,.03)),
                             education = sample.int(6, size=N, replace = TRUE, prob = c(.05,.25,.27,.1,.3,.03)),
                             state = sample.int(5, size=N, replace = TRUE, prob = c(.28,.4,.12,.15,.05)),
                             vote_for = sample.int(4, size=N, replace = TRUE, prob = abs(rnorm(4,0,1))),
                             y = rbinom(N,1,inv_logit_scaled(c(1,-.5,.5)[gender] +
                                                               c(-.4,.4,.5,1.0)[vote_for] +
                                                               c(-.2,-.5,.5,.5,1.0,.4,.7)[age_group] +
                                                               c(-.2,.7,-.5,-.5,.4,.7)[education] +
                                                               c(-.3,-.6,.6,.3,.4)[state])))

### get population
popn <- draw_data(design +
                    declare_sampling(n=50000))

### function to create sample
get_sample <- function(data, n) {
  p = inv_logit_scaled(.5 + c(-1,-2,-1,-.5,1,.5,.1)[data$age_group]+
                         c(2,-.5,0)[data$gender] +
                         c(1,-.5,-1,.5)[data$vote_for] +
                         c(-1,1,-2,-.5,1,.5)[data$education] +
                         c(-.2,-.7,.6,.5,.4)[data$state])
  data[sample.int(nrow(data),n,replace=TRUE,prob = p), ]
}

### function to create approximate population
get_pop_sample <- function(data, n) {
  p = inv_logit_scaled(.5 + c(-.1,-.2,-.1,-.5,.1,.5,.01)[data$age_group]+
                         c(.2,-.5,0)[data$gender] +
                         c(.1,-.5,-.1,.5)[data$vote_for] +
                         c(-1,1,-2,-.5,1,.3)[data$education] +
                         c(-.2,-.5,.6,.5,.3)[data$state])
  data[sample.int(nrow(data),n,replace=TRUE,prob = p), ]
}

### create the survey data
bp_survey <- get_sample(popn,500)

### create the approx.population data
approx_voters_popn <- get_pop_sample(popn, 5000) %>%
  select(-y)

### create margin for all of the survey questions
age_margin <- xtabs(~age_group, data = popn)
gender_margin <- xtabs(~gender, data = popn)
vote_margin <- xtabs(~vote_for, data = popn)
education_margin <- xtabs(~education, data = popn)
state_margin <- xtabs(~state, data = popn)

### add weight to approx.population data
approx_voters_popn$wt = 1

svy_design_approxpopn <- svydesign(ids=~1,weights = ~wt, data= approx_voters_popn)

rkd_obj_approx_popn <- rake(design = svy_design_approxpopn, sample.margins = list(
  ~ gender, ~ age_group, ~ vote_for, ~ education, ~ state), population.margins = list(gender_margin,
                                                               age_margin, vote_margin,
                                                               education_margin, state_margin))
wts_trim_approx_popn <- trimWeights(rkd_obj_approx_popn,
                                    upper = quantile(weights(rkd_obj_approx_popn),.975))
approx_voters_popn$wt <- weights(wts_trim_approx_popn)

#### Voters Survey ####

# Create weights #
bp_survey$wt = 1

svy_design_bp_survey <- svydesign(ids=~1,weights = ~wt, data= bp_survey)

rkd_obj_bp_survey <- rake(design = svy_design_bp_survey, sample.margins = list(
  ~ gender, ~ age_group, ~ vote_for, ~ education, ~ state),
  population.margins = list(gender_margin,age_margin, vote_margin,
                            education_margin, state_margin))
wts_trim_bp_survey <- trimWeights(rkd_obj_bp_survey,
                                   upper = quantile(weights(rkd_obj_bp_survey),.975))
bp_survey$wt <- weights(wts_trim_bp_survey)


# Adjust measurement  #

bp_survey <- bp_survey %>%
  select(-S_inclusion_prob) %>%
  mutate_at(c("age_group", "gender", "vote_for", "education", "state" ,"y"), as.factor) %>%
  mutate(age =fct_recode(age_group, `18-25` = "1", `26-35` = "2",
                         `36-45` = "3",`46-55` = "4",
                         `56-65` = "5", `66-75` = "6", `76-90` = "7"),
         gender = fct_recode(gender, "male" = "1",
                             "female" = "2",
                             "nonbinary" = "3"),
         vote_for = fct_recode(vote_for,
                               "Box Party" = "1",
                               "BP" = "2",
                               "Circle Party" = "3",
                               "CP" = "4"),
         highest_educ = fct_recode(education,
                               "no high school" = "1",
                               "high school" = "2",
                               "some college" = "3",
                               "associates" = "4",
                               "4-year college" = "5",
                               "post-graduate" = "6"),
         state = fct_recode(state,
                               "State A" = "1",
                               "State B" = "2",
                               "State C" = "3",
                               "State D" = "4",
                               "State E" = "5"),
         y = fct_recode(y,"no" = "0", "yes" = "1")) %>%
  select(c(age, gender, vote_for, highest_educ, state, y, wt))


approx_voters_popn <- approx_voters_popn %>%
  select(-S_inclusion_prob) %>%
  mutate_at(c("age_group","gender","vote_for", "education", "state"), as.factor) %>%
  mutate(age_group =fct_recode(age_group, `18-35` = "1", `18-35` = "2",
                               `36-55` = "3",`36-55` = "4",
                               `56-65` = "5", `66+` = "6", `66+` = "7"),
         gender = fct_recode(gender, "m" = "1", "f" = "2", "nb" = "3"),
         vote_pref = fct_recode(vote_for, "BP" = "1",
                                "BP" = "2",
                                "CP" = "3",
                                "CP" = "4"),
         education = fct_recode(education,
                                "no high school" = "1",
                                "high school" = "2",
                                "some college" = "3",
                                "some college" = "4",
                                "4-years college" = "5",
                                "post-grad" = "6"),
         state = fct_recode(state,
                            "A" = "1",
                            "B" = "2",
                            "C" = "3",
                            "D" = "4",
                            "E" = "5")) %>%
  select(c(age_group, gender, vote_pref, wt, education, state))

usethis::use_data(bp_survey, approx_voters_popn, overwrite = TRUE)

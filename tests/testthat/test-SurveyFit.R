samp_obj <- SurveyData$new(
  data = feline_survey,
  questions = list(
    age1 = "Please identify your age group",
    gender = "Please select your gender",
    pet_own = "Which pet do you own?",
    y = "Response"
  ),
  responses = list(
    age1 = levels(feline_survey$age1),
    gender = levels(feline_survey$gender),
    pet_own = levels(feline_survey$pet_own),
    y = c("no","yes")
  ),
  weights = feline_survey$wt,
  design = formula("~.")
)

popn_obj <- SurveyData$new(
  data = approx_popn,
  questions = list(
    age2 = "Which age group are you?",
    gender = "Gender?",
    pet_pref = "Which pet would you like to own?"
  ),
  responses = list(
    gender = levels(approx_popn$gender),
    age2 = levels(approx_popn$age2),
    pet_pref = levels(approx_popn$pet_pref)
  ),
  weights = approx_popn$wt,
  design = formula("~.")
)

q1 <- SurveyQuestion$new(
  name = "age",
  col_names = c("age1","age2"),
  values_map = list(
    "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
    "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
  )
)
q2 <- SurveyQuestion$new(
  name = "pet",
  col_names = c("pet_own","pet_pref"),
  values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
)
q3 <- SurveyQuestion$new(
  name = "gender",
  col_names = c("gender","gender"),
  values_map = data.frame("male" = "m","female" = "f", "nonbinary" = "nb")
)
ex_map <- SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1,q2,q3)
ex_map$mapping()
ex_map$tabulate()

suppressWarnings(
  fit_glm <- ex_map$fit(
    fun = stats::glm,
    formula = y ~ age + gender,
    family = "binomial"
  )
)

if (requireNamespace("rstanarm", quietly = TRUE)) {
  suppressWarnings(
    fit_stan_glmer <- ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = "binomial",
      iter = 10,
      chains = 1,
      refresh = 0
    )
  )
  suppressWarnings(
    fit_stan_glm <- ex_map$fit(
      fun = rstanarm::stan_glm,
      formula = y ~ age + gender,
      family = "binomial",
      iter = 10,
      chains = 1,
      refresh = 0
    )
  )
}

if (requireNamespace("brms", quietly = TRUE)) {
  suppressWarnings(
    fit_brms <- ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = "bernoulli",
      iter = 10,
      chains = 1,
      refresh = 0
    )
  )
}

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_glmer <- ex_map$fit(
    fun = lme4::glmer,
    formula = y ~ (1|age) + (1|gender),
    family = "binomial"
  )
}

test_that("population_predict runs without errors",{
  skip_if_not_installed("rstanarm")
  expect_error(fit_stan_glmer$population_predict(), regexp = NA)
  expect_error(fit_stan_glm$population_predict(), regexp = NA)

  skip_if_not_installed("lme4")
  expect_error(fit_glmer$population_predict(), regexp = NA)

  skip_if_not_installed("brms")
  expect_error(fit_brms$population_predict(), regexp = NA)
})

test_that("population_predict errors if custom fun required but not specified", {
  expect_error(
    fit_glm$population_predict(),
    "Custom population_predict method required"
  )
})

test_that("Aggregate runs without errors",{
  skip_if_not_installed("rstanarm")
  expect_error(fit_stan_glmer$aggregate(fit_stan_glmer$population_predict()),
               regexp = NA)
  expect_error(fit_stan_glm$aggregate(fit_stan_glm$population_predict()),
               regexp = NA)

  skip_if_not_installed("lme4")
  expect_error(fit_glmer$aggregate(fit_glmer$population_predict()),
               regexp = NA)

  skip_if_not_installed("brms")
  expect_error(fit_brms$aggregate(fit_brms$population_predict()),
               regexp = NA)
})

# the objects in this test still need to be created
# test_that("Populations are within acceptable tolerance of previous runs (+/- 2% points)",{
#   skip_if_not_installed("rstanarm")
#   expect_lt(mean(popn_ests_rstanarm$value), .72 + .02)
#   expect_gt(mean(popn_ests_rstanarm$value), .72 - .02)
#
#   skip_if_not_installed("lme4")
#   #Benchmark to a run from previous, so different benchmark values
#   expect_lt(mean(popn_ests_lme4$value), .68 + .02)
#   expect_gt(mean(popn_ests_lme4$value), .68 - .02)
#
#   skip_if_not_installed("brms")
#   #Benchmark to a run from previous, so different benchmark values
#   expect_lt(mean(popn_ests_brms$value), .85 + .02)
#   expect_gt(mean(popn_ests_brms$value), .85 - .02)
# })


test_that("summary method throws correct errors", {
  expect_error(
    fit_stan_glmer$summary(stats = list(f = function(x) c(x, x))),
    "Function 'f' returns more than one value"
  )
  expect_error(
    fit_stan_glmer$summary(stats = list("mean", function(x) c)),
    "Name must be provided for element 2 of 'stats'"
  )
})


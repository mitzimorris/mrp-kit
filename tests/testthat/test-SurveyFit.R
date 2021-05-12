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

test_that("population_predict returns correct objects",{
  # expect 5 draws (because iter = 10 above) for columns
  # expect same number of rows as poststrat data
  expected_dims <- c(nrow(ex_map$poststrat_data()), 5)

  skip_if_not_installed("rstanarm")
  expect_equal(dim(fit_stan_glmer$population_predict()), expected_dims)
  expect_equal(dim(fit_stan_glm$population_predict()), expected_dims)

  skip_if_not_installed("brms")
  expect_equal(dim(fit_brms$population_predict()), expected_dims)

  skip_if_not_installed("lme4")
  expect_equal(dim(fit_glmer$population_predict(nsamples = 5)), expected_dims)
})

test_that("population_predict errors if custom fun required but not specified", {
  expect_error(
    fit_glm$population_predict(),
    "Custom population_predict method required"
  )
})

test_that("Aggregate (to population) returns correct objects", {
  expected_dims <- c(5, 1) # 5 = 10 iter / 2

  skip_if_not_installed("rstanarm")
  x <- fit_stan_glmer$aggregate(fit_stan_glmer$population_predict())
  expect_s3_class(x, "data.frame")
  expect_named(x, "value")
  expect_equal(dim(x), expected_dims)

  x <- fit_stan_glmer$aggregate(fit_stan_glm$population_predict())
  expect_s3_class(x, "data.frame")
  expect_named(x, "value")
  expect_equal(dim(x), expected_dims)

  skip_if_not_installed("brms")
  x <- fit_brms$aggregate(fit_brms$population_predict())
  expect_s3_class(x, "data.frame")
  expect_named(x, "value")
  expect_equal(dim(x), expected_dims)

  skip_if_not_installed("lme4")
  x <- fit_glmer$aggregate(fit_glmer$population_predict(nsamples = 5))
  expect_s3_class(x, "data.frame")
  expect_named(x, "value")
  expect_equal(dim(x), expected_dims)
})

test_that("Aggregate (by variable level) returns correct objects", {
  expected_dims <- c(5 * nlevels(popn_obj$survey_data()$age), 3) # 5 = 10 iter / 2
  expected_names <- c("age", "draw", "value")

  skip_if_not_installed("rstanarm")
  x <- fit_stan_glmer$aggregate(fit_stan_glmer$population_predict(), by = "age")
  expect_s3_class(x, "data.frame")
  expect_named(x, expected_names)
  expect_equal(dim(x), expected_dims)

  x <- fit_stan_glm$aggregate(fit_stan_glm$population_predict(), by = "age")
  expect_s3_class(x, "data.frame")
  expect_named(x, expected_names)
  expect_equal(dim(x), expected_dims)

  skip_if_not_installed("brms")
  x <- fit_brms$aggregate(fit_brms$population_predict(), by = "age")
  expect_s3_class(x, "data.frame")
  expect_named(x, expected_names)
  expect_equal(dim(x), expected_dims)

  skip_if_not_installed("lme4")
  x <- fit_glmer$aggregate(fit_glmer$population_predict(nsamples = 5), by = "age")
  expect_s3_class(x, "data.frame")
  expect_named(x, expected_names)
  expect_equal(dim(x), expected_dims)
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


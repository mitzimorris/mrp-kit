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
ex_map <- SurveyMap$new(samp_obj, popn_obj, q1,q2,q3)
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

test_that("fit objects have correct R6 class", {
  expect_r6_class(fit_stan_glmer, "SurveyFit")
  expect_r6_class(fit_stan_glm, "SurveyFit")
  expect_r6_class(fit_brms, "SurveyFit")
  expect_r6_class(fit_glmer, "SurveyFit")
  expect_r6_class(fit_glm, "SurveyFit")
})

test_that("map returns SurveyMap", {
  expect_r6_class(fit_stan_glmer$map(), "SurveyMap")
  expect_r6_class(fit_stan_glm$map(), "SurveyMap")
  expect_r6_class(fit_brms$map(), "SurveyMap")
  expect_r6_class(fit_glmer$map(), "SurveyMap")
  expect_r6_class(fit_glm$map(), "SurveyMap")
})

test_that("fit returns fitted model object", {
  expect_s3_class(fit_stan_glmer$fit(), "stanreg")
  expect_s3_class(fit_stan_glm$fit(), "stanreg")
  expect_s3_class(fit_brms$fit(), "brmsfit")
  expect_s4_class(fit_glmer$fit(), "glmerMod")
  expect_s3_class(fit_glm$fit(), "glm")
})

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

test_that("plot returns ggplot object", {
  popn <- fit_stan_glmer$aggregate(fit_stan_glmer$population_predict())
  expect_s3_class(fit_stan_glmer$plot(popn, weights = TRUE), "ggplot")
  expect_s3_class(fit_stan_glmer$plot(popn, weights = FALSE), "ggplot")

  by_age <- fit_stan_glmer$aggregate(fit_stan_glmer$population_predict(), by = "age")
  expect_s3_class(fit_stan_glmer$plot(by_age, weights = TRUE), "ggplot")
  expect_s3_class(fit_stan_glmer$plot(by_age, weights = FALSE), "ggplot")
})

test_that("plot appearance hasn't changed", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")

  # need to use deterministic inputs to the plots to test appearance changes
  popn <- data.frame(value = c(0.60, 0.65, 0.70, 0.75, 0.80))
  by_age <- data.frame(age = factor(rep(1:4, 5), labels = popn_obj$responses()$age),
                       draw = rep(1:5, each = 4),
                       value = c(0.299514804640785, 0.217973976861686, 0.258871184848249, 0.271914651058614,
                                 0.31649006572552, 0.697755558183417, 0.209188556019217, 0.431414544349536,
                                 0.74861360588111, 0.209171398542821, 0.385790600813925, 0.45710161360912,
                                 0.537278639385477, 0.752459280192852, 0.510153451515362, 0.644023981876671,
                                 0.570707685500383, 0.323449705773965, 0.505213424470276, 0.236637408705428)
                       )

  vdiffr::expect_doppelganger(
    "plot-population",
    fit_glmer$plot(popn, weights = FALSE),
    path = "SurveyFit"
  )
  vdiffr::expect_doppelganger(
    "plot-population-weights",
    fit_glmer$plot(popn, weights = TRUE),
    path = "SurveyFit"
  )

  vdiffr::expect_doppelganger(
    "plot-group",
    fit_glmer$plot(by_age, weights = FALSE),
    path = "SurveyFit"
  )

  vdiffr::expect_doppelganger(
    "plot-group-weights",
    fit_glmer$plot(by_age, weights = TRUE),
    path = "SurveyFit"
  )
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


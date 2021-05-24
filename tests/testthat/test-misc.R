test_that("sim_posterior_probs throws error if not given glmerMod object", {
  skip_if_not_installed("merTools")
  skip_if_not_installed("rstanarm")
  skip_if_not_installed("lme4")

  example_newdata <- feline_survey[sample.int(nrow(feline_survey),150),]
  rstanarm_fit <- suppressWarnings(
    rstanarm::stan_glmer(y ~ (1|age1) + (1|gender),
                         data = feline_survey[1:300,],
                         family = binomial(link = "logit"),
                         refresh = 0, iter = 10, chains = 1)
  )

  glmer_fit <- lme4::glmer(y ~ (1|age1) + (1|gender), data = feline_survey,
                           family = binomial(link = "logit"))

  glm_fit <- stats::glm(y ~ age1 + gender, data = feline_survey,
                        family = binomial(link = "logit"))
  lm_fit <- stats::lm(as.numeric(y) ~ age1 + gender, data = feline_survey)

  expect_error(
    suppressWarnings(sim_posterior_probs(rstanarm_fit,newdata = example_newdata)),
    "Object must have class 'glmerMod'."
  )
  expect_error(
    suppressWarnings(sim_posterior_probs(glm_fit,newdata = example_newdata)),
    "Object must have class 'glmerMod'."
  )
  expect_error(
    suppressWarnings(sim_posterior_probs(lm_fit,newdata = example_newdata)),
    "Object must have class 'glmerMod'."
  )
  expect_equal(
    dim(sim_posterior_probs(glmer_fit,newdata = example_newdata, nsamples = 10)),
    c(150,10)
  )
  expect_equal(
    sum(is.na(dim(sim_posterior_probs(glmer_fit,newdata = example_newdata, nsamples = 10)))),
    0)
})

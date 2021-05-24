suppressWarnings({
  samp <- SurveyData$new(data = feline_survey, weights = feline_survey$wt)
  popn <- SurveyData$new(data = approx_popn, weights = approx_popn$wt)
})
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

test_that("object has correct R6 class", {
  expect_r6_class(SurveyMap$new(samp, popn), "SurveyMap")
})


test_that("error thrown if inputs are not SurveyData objects", {
  expect_error(
    SurveyMap$new("ABC", popn),
    "'sample' must be a SurveyData object"
  )
  expect_error(
    SurveyMap$new(samp, "ABC"),
    "'population' must be a SurveyData object"
  )
})

test_that("initializing with 0 questions doesn't error", {
  expect_silent(SurveyMap$new(samp, popn))
})

test_that("initializing with >0 questions doesn't error", {
  expect_silent(SurveyMap$new(samp, popn, q1, q2))
})

test_that("add() errors if name already exists", {
  ex_map <- SurveyMap$new(samp, popn, q1)
  expect_error(
    ex_map$add(q1),
    "Survey label 'age' already defined"
  )
})

test_that("validate creates correct levels (example1)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26-45", "46+"), 10)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26-45","46+"),
      y = c("no","yes")
    ),
    weights = rnorm(30 ,0, 1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45", "46+"), 50))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45", "46+")
    ),
    weights = rnorm(100, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-45", "26-45" = "18-45", "46+" = "46+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age,
                  c('18-45', "46+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c('18-45', "46+"))
})

test_that("validate creates correct levels (example2)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26-45","46+"), 10)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26-45","46+"),
      y = c("no","yes")
    ),
    weights = rnorm(30,0,1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45", "46+"), 50))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45","46+")
    ),
    weights = rnorm(100, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-45", "26-45" = "18-45", "46+" = "46+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age[1:10],
                  c("18-45", "18-45", "46+", "18-45", "18-45", "46+", "18-45", "18-45", "46+", "18-45"))
  expect_setequal(ex_map$mapped_population_data()$age[1:10],
                  c("18-45", "46+", "18-45", "46+", "18-45", "46+", "18-45", "46+", "18-45", "46+"))
})


test_that("validate creates correct levels (example3)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26+"), 15)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26+"),
      y = c("no", "yes")
    ),
    weights = rnorm(30, 0, 1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-34", "35+"), 30))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-34", "35+")
    ),
    weights = rnorm(90, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c('18-25', "26+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c('18-25', "26+"))
})

test_that("validate creates correct levels (example4)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26+"), 15)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26+"),
      y = c("no", "yes")
    ),
    weights = rnorm(30, 0, 1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-34", "35+"), 30))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-34", "35+")
    ),
    weights = rnorm(90, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list("18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age[7:16],
                  c("18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+"))
  expect_setequal(ex_map$mapped_population_data()$age[7:16],
                  c("18-25", "26+", "26+", "18-25", "26+", "26+", "18-25", "26+", "26+", "18-25"))
})

test_that("validate creates correct levels (example5)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26-30", "31-40", "41-55", "56+"), 20)),
                      y = factor(rbinom(100, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26-30", "31-40", "41-55", "56+"),
      y = c("no","yes")
    ),
    weights = rnorm(100, 0, 1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-35", "36-45", "46-55", "56+"), 40))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-35", "36-45", "46-55", "56+")
    ),
    weights = rnorm(200, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-25", "26-30" = "26-35", "31-40" = "26-35", "31-40" = "36-45",
      "41-55" = "36-45", "41-55"="46-55", "56+"="56+"
    )
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
})

test_that("validate creates correct levels (example6)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("A","A","B","C","D","D","E"), 20)),
                      y = factor(rbinom(140, 1, .5), levels = c("no", "yes"))),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("A", "B", "C", "D", "E"),
      y = c("no", "yes")
    ),
    weights = rnorm(140, 0, 1),
    design = formula("~.")
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep( c("Z","Y","Y","C","X","Q","Q"), 40))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 =  c("Z","Y","C","X","Q")
    ),
    weights = rnorm(280, 0, 1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "A" = "Z", "A" = "Y","B" = "Y","C" = "C",
      "D" = "X","D"="Q","E"="Q")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q1))
  ex_map$validate()
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c("A + B","C","D + E"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c("A + B","C","D + E"))
})


# This test is currently failing. Should it be?
# Expected match: "Predictor variables not known in population."
# Actual message: "Not all variables available in the data. Missing vars:  gender"
test_that("Error if predictor vars not included in poststrat matrix",{
  ex_map <- SurveyMap$new(samp, popn, q1, q2)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map <-  SurveyMap$new(samp, popn, q1, q2)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map <-  SurveyMap$new(samp, popn, q1, q2)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
})

test_that("Error if vars not included in data",{
  ex_map <-  SurveyMap$new(samp, popn, q1,q2,q3)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender) + pineapple,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender) +pineapple ,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data. ", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender) +pineapple,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data.", fixed = TRUE
  )

  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = pineapple ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = pineapple ~ (1|age) + (1|gender) ,
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = pineapple ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
})

test_that("Error if not fitting a bernoulli/binomial model", {
  ex_map <- SurveyMap$new(samp, popn, q1, q2, q3)
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )

})

test_that("Error if data hasn't been mapped yet",{
  ex_map <-  SurveyMap$new(samp, popn, q1, q2, q3)
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
})

test_that("Error if poststrat matrix hasn't been created yet",{
  ex_map <- SurveyMap$new(samp, popn, q1, q2, q3)
  ex_map$mapping()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
})

test_that("Error if data is given as input",{
  ex_map <- SurveyMap$new(samp, popn, q1, q2, q3)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
})

test_that("Warning if fitting using packages that are not lme4, brms, rstanarm ", {
  ex_map <- SurveyMap$new(samp, popn, q1, q2, q3)
  ex_map$mapping()
  ex_map$tabulate()
  expect_warning(
    ex_map$fit(
      fun = stats::glm,
      formula = y ~ age + gender,
      family = "binomial"
    ),
    "Only rstanarm, brms and lme4 are supported natively.", fixed = TRUE
  )
})

test_that("Model fits do not cause errors if specified correctly",{
  ex_map <- SurveyMap$new(samp, popn, q1, q2, q3)
  ex_map$mapping()
  ex_map$tabulate()
  skip_if_not_installed("rstanarm")
  expect_error(suppressWarnings(ex_map$fit(
    fun = rstanarm::stan_glmer,
    formula = y ~ (1|age) + (1|gender),
    refresh = 0,
    cores = 1,
    chains = 1,
    iter = 10,
    family = binomial(link = "logit")
  )), regexp = NA)
  expect_error(suppressWarnings(ex_map$fit(
    fun = rstanarm::stan_glm,
    formula = y ~ age + gender,
    refresh = 0,
    chains = 1,
    cores = 1,
    iter = 10,
    family = "binomial"
  )), regexp = NA)

  skip_if_not_installed("brms")
  expect_error(suppressWarnings(ex_map$fit(
    fun = brms::brm,
    formula = y ~ (1|age) + (1|gender),
    refresh = 0,
    iter = 10,
    chains = 1,
    cores = 1,
    family = "bernoulli"
  )), regexp = NA)

  skip_if_not_installed("lme4")
  expect_error(ex_map$fit(
    fun = lme4::glmer,
    formula = y ~ (1|age) + (1|gender),
    family = "binomial"
  ), regexp = NA)
})

test_that("mapped_sample_data and mapped_population_data work", {
  ex_map <-  SurveyMap$new(samp, popn, q1, q2, q3)
  empty_samp <- ex_map$mapped_sample_data(key = FALSE)
  empty_pop <- ex_map$mapped_population_data(key = FALSE)
  expect_equal(dim(empty_samp), c(nrow(samp$survey_data()), 0))
  expect_equal(dim(empty_pop), c(nrow(popn$survey_data()), 0))

  just_key_samp <- ex_map$mapped_sample_data(key = TRUE)
  just_key_pop <- ex_map$mapped_population_data(key = TRUE)
  expect_named(just_key_samp, ".key")
  expect_named(just_key_pop, ".key")
  expect_equal(dim(just_key_samp), c(nrow(samp$survey_data()), 1))
  expect_equal(dim(just_key_pop), c(nrow(popn$survey_data()), 1))

  ex_map$mapping()
  expect_named(ex_map$mapped_sample_data(key = FALSE), c("age", "pet", "gender"))
  expect_named(ex_map$mapped_population_data(key = FALSE), c("age", "pet", "gender"))
})

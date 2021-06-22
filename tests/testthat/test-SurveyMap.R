suppressWarnings({
  samp <- SurveyData$new(data = feline_survey, weights = "wt")
  popn <- SurveyData$new(data = approx_popn, weights = "wt")
})
q_age <- QuestionMap$new(
  name = "age",
  col_names = c("age1","age2"),
  values_map = list(
    "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
    "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
  )
)
q_pet <- QuestionMap$new(
  name = "pet",
  col_names = c("pet_own","pet_pref"),
  values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
)
q_gender <- QuestionMap$new(
  name = "gender",
  col_names = c("gender","gender"),
  values_map = data.frame("male" = "m","female" = "f", "nonbinary" = "nb")
)

test_that("object has correct R6 class", {
  expect_r6_class(SurveyMap$new(samp, popn), "SurveyMap")
})

test_that("print output hasn't changed", {
  expect_known_output(
    print(SurveyMap$new(samp, popn)),
    file = test_path("answers/SurveyMap-print-empty")
  )

  expect_known_output(
    print(suppressWarnings(SurveyMap$new(samp, popn, q_age))),
    file = test_path("answers/SurveyMap-print-1-question")
  )
  expect_known_output(
    print(suppressWarnings(SurveyMap$new(samp, popn, q_age))),
    file = test_path("answers/SurveyMap-print-1-question")
  )
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
  ex_map <- SurveyMap$new(samp, popn)
  expect_length(ex_map$item_map(), 0)
})

test_that("initializing with >0 questions doesn't error", {
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
  expect_length(ex_map$item_map(), 3)
})

test_that("add() errors if name already exists", {
  ex_map <- suppressWarnings(SurveyMap$new(samp, popn, q_age))
  expect_error(
    ex_map$add(q_age),
    "Survey label 'age' already defined"
  )
})

test_that("validate creates correct levels (example1)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26-45", "46+"), 10)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(30 ,0, 1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26-45","46+"),
      y = c("no","yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45", "46+"), 50)), wt = rnorm(100, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45", "46+")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-45", "26-45" = "18-45", "46+" = "46+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age,
                  c('18-25 + 26-45', "46+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c('18-25 + 26-45', "46+"))
})

test_that("validate creates correct levels (example2)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26-45","46+"), 10)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(30,0,1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26-45","46+"),
      y = c("no","yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45", "46+"), 50)),
                      wt = rnorm(100, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45","46+")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-45", "26-45" = "18-45", "46+" = "46+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age[1:10],
                  c("18-25 + 26-45", "18-25 + 26-45", "46+","18-25 + 26-45", "18-25 + 26-45", "46+", "18-25 + 26-45",
                    "18-25 + 26-45", "46+", "18-25 + 26-45"))
  expect_setequal(ex_map$mapped_population_data()$age[1:10],
                  c("18-25 + 26-45", "46+", "18-25 + 26-45", "46+", "18-25 + 26-45", "46+", "18-25 + 26-45", "46+", "18-25 + 26-45", "46+"))
})


test_that("validate creates correct levels (example3)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26+"), 15)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(30, 0, 1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26+"),
      y = c("no", "yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-34", "35+"), 30)),
                      wt = rnorm(90, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-34", "35+")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1", "age2"),
    values_map = list("18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c('18-25', "26+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c('18-25', "26+"))
})

test_that("validate creates correct levels (example4)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26+"), 15)),
                      y = factor(rbinom(30, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(30, 0, 1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26+"),
      y = c("no", "yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-34", "35+"), 30)),
                      wt = rnorm(90, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-34", "35+")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list("18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_setequal(ex_map$mapped_sample_data()$age[7:16],
                  c("18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+"))
  expect_setequal(ex_map$mapped_population_data()$age[7:16],
                  c("18-25", "26+", "26+", "18-25", "26+", "26+", "18-25", "26+", "26+", "18-25"))
})

test_that("validate creates correct levels (example5)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25", "26-30", "31-40", "41-55", "56+"), 20)),
                      y = factor(rbinom(100, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(100, 0, 1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25", "26-30", "31-40", "41-55", "56+"),
      y = c("no","yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25", "26-35", "36-45", "46-55", "56+"), 40)),
                      wt = rnorm(200, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25", "26-35", "36-45", "46-55", "56+")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-25", "26-30" = "26-35", "31-40" = "26-35", "31-40" = "36-45",
      "41-55" = "36-45", "41-55"="46-55", "56+"="56+"
    )
  )
  expect_silent(ex_map <- SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
  expect_setequal(ex_map$mapped_population_data()$age,
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
})

test_that("validate creates correct levels (example6)", {
  samp <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("A","A","B","C","D","D","E"), 20)),
                      y = factor(rbinom(140, 1, .5), levels = c("no", "yes")),
                      wt = rnorm(140, 0, 1)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("A", "B", "C", "D", "E"),
      y = c("no", "yes")
    ),
    weights = "wt"
  )

  popn <- SurveyData$new(
    data = data.frame(age2 = factor(rep( c("Z","Y","Y","C","X","Q","Q"), 40)),
                      wt = rnorm(280, 0, 1)),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 =  c("Z","Y","C","X","Q")
    ),
    weights = "wt"
  )

  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "A" = "Z", "A" = "Y","B" = "Y","C" = "C",
      "D" = "X","D"="Q","E"="Q")
  )
  ex_map <- SurveyMap$new(samp, popn, q_age)
  ex_map$mapping()
  expect_setequal(levels(ex_map$mapped_sample_data()$age),
                  c("A + B","C","D + E"))
  expect_setequal(levels(ex_map$mapped_population_data()$age),
                  c("A + B","C","D + E"))
})

test_that("validate errors if NAs in adjustment variables in sample data", {
  d <- feline_survey
  d$age1[3] <- NA
  suppressWarnings(samp_NA <- SurveyData$new(data = d))
  expect_error(
    SurveyMap$new(samp_NA, popn, q_age, q_pet, q_gender),
    "NAs not allowed in variables mapped between sample and population"
  )
})

test_that("validate errors if NAs in population data", {
  d <- approx_popn
  d$age2[3] <- NA
  suppressWarnings(popn_NA <- SurveyData$new(data = d, weights = "wt"))
  expect_error(
    SurveyMap$new(samp, popn_NA, q_age, q_pet, q_gender),
    "NAs not allowed in variables mapped between sample and population"
  )
})

test_that("Error if predictor vars not included in poststrat matrix",{
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet,q_gender)
  ex_map$mapping()
  ex_map$tabulate()

  # change the internal poststrat to trigger error
  # This is a very unlikely thing to happen, but just in case!
  ex_map$.__enclos_env__$private$poststrat_data_<-ex_map$.__enclos_env__$private$poststrat_data_[c("age","pet","N_j")]
  skip_if_not_installed("lme4")
  expect_error(
    ex_map$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map <-  SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
  ex_map$mapping()
  ex_map$tabulate()
  ex_map$.__enclos_env__$private$poststrat_data_<-ex_map$.__enclos_env__$private$poststrat_data_[c("age","pet","N_j")]
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map <-  SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
  ex_map$mapping()
  ex_map$tabulate()
  ex_map$.__enclos_env__$private$poststrat_data_<-ex_map$.__enclos_env__$private$poststrat_data_[c("age","pet","N_j")]
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

test_that("Error if vars not included in data", {
  ex_map <-  SurveyMap$new(samp, popn, q_age,q_pet,q_gender)
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
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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
  ex_map <-  SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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

test_that("Error if calling tabulate before mapping", {
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
  expect_error(
    ex_map$tabulate(),
    "Please call the mapping() method",
    fixed = TRUE
  )
})

test_that("Error if poststrat matrix hasn't been created yet",{
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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

test_that("Warning if missingness in outcome", {
  d <- feline_survey
  d$y[3] <- NA
  suppressWarnings(
    samp2 <- SurveyData$new(
      data = d,
      weights = "wt"
    )
  )
  ex_map <- SurveyMap$new(samp2, popn, q_age, q_pet, q_gender)
  ex_map$mapping()
  ex_map$tabulate()
  expect_warning(
    fit <- ex_map$fit(lme4::glmer, formula = y ~ (1|age), family = "binomial"),
    "Outcome variable has missing values that may be dropped by the model fitting package."
  )
  expect_r6_class(fit, "SurveyFit")
})

test_that("Model fits do not cause errors if specified correctly",{
  ex_map <- SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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
  ex_map <-  SurveyMap$new(samp, popn, q_age, q_pet, q_gender)
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

test_that("tabulate doesn't error if no weights were specified", {
  suppressWarnings({
    samp <- SurveyData$new(feline_survey)
    popn <- SurveyData$new(approx_popn)
  })
  q_age <- QuestionMap$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-35", "26-35" = "18-35","36-45" = "36-55",
      "46-55" = "36-55", "56-65" = "56-65", "66-75" = "66+", "76-90" = "66+"
    )
  )
  ex_map <- suppressWarnings(SurveyMap$new(samp, popn, q_age))
  ex_map$mapping()
  expect_equal(nrow(ex_map$poststrat_data()), 0)
  expect_silent(ex_map$tabulate())
  expect_equal(dim(ex_map$poststrat_data()), c(4, 2))
})


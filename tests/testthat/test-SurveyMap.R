# objects to use in the tests
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


test_that("error thrown if inputs are not SurveyData objects", {
  expect_error(
    SurveyMap$new("ABC", popn_obj),
    "samp_obj must be a SurveyData object"
  )
  expect_error(
    SurveyMap$new(samp_obj, "ABC"),
    "popn_obj must be a SurveyData object"
  )
})

test_that("initializing with 0 questions doesn't error", {
  expect_silent(SurveyMap$new(samp_obj, popn_obj))
})

test_that("initializing with >0 questions doesn't error", {
  expect_silent(SurveyMap$new(samp_obj, popn_obj, q1, q2))
})

test_that("add() errors if name already exists", {
  x <- SurveyMap$new(samp_obj, popn_obj, q1)
  expect_error(
    x$add(q1),
    "Survey label 'age' already defined"
  )
})

test_that("validate creates correct levels (example1)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26-45","46+"),10)),
                      y = rbinom(30,1,.5)),
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

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45","46+"),50))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45","46+")
    ),
    weights = rnorm(100,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-45", "26-45" = "18-45","46+" = "46+")
  )
  expect_silent(ex_mapping1 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping1$validate()
  ex_mapping1$mapping()
  expect_setequal(ex_mapping1$.__enclos_env__$private$samp_obj_$mapped_data()$age,
                  c('18-45', "46+"))
  expect_setequal(ex_mapping1$.__enclos_env__$private$popn_obj_$mapped_data()$age,
                  c('18-45', "46+"))
})

test_that("validate creates revels the data correctly (example1)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26-45","46+"),10)),
                      y = rbinom(30,1,.5)),
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

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-45","46+"),50))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-45","46+")
    ),
    weights = rnorm(100,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-45", "26-45" = "18-45","46+" = "46+")
  )
  expect_silent(ex_mapping1 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping1$validate()
  ex_mapping1$mapping()
  expect_setequal(ex_mapping1$.__enclos_env__$private$samp_obj_$mapped_data()$age[1:10],
                  c("18-45", "18-45", "46+", "18-45", "18-45", "46+", "18-45", "18-45", "46+", "18-45"))
  expect_setequal(ex_mapping1$.__enclos_env__$private$popn_obj_$mapped_data()$age[1:10],
                  c("18-45", "46+", "18-45", "46+", "18-45", "46+", "18-45", "46+", "18-45", "46+"))
})


test_that("validate creates correct levels (example2)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26+"),15)),
                      y = rbinom(30,1,.5)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26+"),
      y = c("no","yes")
    ),
    weights = rnorm(30,0,1),
    design = formula("~.")
  )

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25","26-34","35+"),30))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25","26-34","35+")
    ),
    weights = rnorm(90,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_mapping2 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping2$validate()
  ex_mapping2$mapping()
  expect_setequal(levels(ex_mapping2$.__enclos_env__$private$samp_obj_$mapped_data()$age),
                  c('18-25', "26+"))
  expect_setequal(ex_mapping2$.__enclos_env__$private$popn_obj_$mapped_data()$age,
                  c('18-25', "26+"))
})

test_that("validate creates revels the data correctly (example2)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26+"),15)),
                      y = rbinom(30,1,.5)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26+"),
      y = c("no","yes")
    ),
    weights = rnorm(30,0,1),
    design = formula("~.")
  )

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25","26-34","35+"),30))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25","26-34","35+")
    ),
    weights = rnorm(90,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-25", "26+" = "26-34","26+" = "35+")
  )
  expect_silent(ex_mapping2 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping2$validate()
  ex_mapping2$mapping()
  expect_setequal(ex_mapping2$.__enclos_env__$private$samp_obj_$mapped_data()$age[7:16],
                  c("18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+", "18-25", "26+"))
  expect_setequal(ex_mapping2$.__enclos_env__$private$popn_obj_$mapped_data()$age[7:16],
                  c("18-25", "26+", "26+", "18-25", "26+", "26+", "18-25", "26+", "26+", "18-25"))
})

test_that("validate creates correct levels (example3)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("18-25","26-30","31-40","41-55","56+"),20)),
                      y = rbinom(100,1,.5)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("18-25","26-30","31-40","41-55","56+"),
      y = c("no","yes")
    ),
    weights = rnorm(100,0,1),
    design = formula("~.")
  )

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep(c("18-25","26-35","36-45","46-55","56+"),40))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 = c("18-25","26-35","36-45","46-55","56+")
    ),
    weights = rnorm(200,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "18-25" = "18-25", "26-30" = "26-35","31-40" = "26-35","31-40" = "36-45",
      "41-55" = "36-45","41-55"="46-55","56+"="56+")
  )
  expect_silent(ex_mapping3 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping3$validate()
  ex_mapping3$mapping()
  expect_setequal(levels(ex_mapping3$.__enclos_env__$private$samp_obj_$mapped_data()$age),
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
  expect_setequal(ex_mapping3$.__enclos_env__$private$popn_obj_$mapped_data()$age,
                  c("18-25","26-30 + 31-40 + 41-55","56+"))
})

test_that("validate creates correct levels (example4)", {
  samp_obj <- SurveyData$new(
    data = data.frame(age1 = factor(rep(c("A","A","B","C","D","D","E"),20)),
                      y = rbinom(140,1,.5)),
    questions = list(
      age1 = "Please identify your age group",
      y = "Response"
    ),
    responses = list(
      age1 = c("A","A","B","C","D","D","E"),
      y = c("no","yes")
    ),
    weights = rnorm(140,0,1),
    design = formula("~.")
  )

  popn_obj <- SurveyData$new(
    data = data.frame(age2 = factor(rep( c("Z","Y","Y","C","X","Q","Q"),40))),
    questions = list(
      age2 = "Which age group are you?"
    ),
    responses = list(
      age2 =  c("Z","Y","Y","C","X","Q","Q")
    ),
    weights = rnorm(280,0,1),
    design = formula("~.")
  )

  q1 <- SurveyQuestion$new(
    name = "age",
    col_names = c("age1","age2"),
    values_map = list(
      "A" = "Z", "A" = "Y","B" = "Y","C" = "C",
      "D" = "X","D"="Q","E"="Q")
  )
  expect_silent(ex_mapping4 <- SurveyMap$new(samp_obj, popn_obj, q1))
  ex_mapping4$validate()
  ex_mapping4$mapping()
  expect_setequal(levels(ex_mapping4$.__enclos_env__$private$samp_obj_$mapped_data()$age),
                  c("A + B","C","D + E"))
  expect_setequal(ex_mapping4$.__enclos_env__$private$popn_obj_$mapped_data()$age,
                  c("A + B","C","D + E"))
})

# This test is currently failing. Should it be?
# Expected match: "Predictor variables not known in population."
# Actual message: "Not all variables available in the data. Missing vars:  gender"
test_that("Error if predictor vars not included in poststrat matrix",{
  ex_map1 <-  SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2)
  ex_map1$mapping()
  ex_map1$tabulate()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map1 <-  SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2)
  ex_map1$mapping()
  ex_map1$tabulate()
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
  ex_map1 <-  SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2)
  ex_map1$mapping()
  ex_map1$tabulate()
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Predictor variables not known in population.", fixed = TRUE
  )
})

test_that("Error if vars not included in data",{
  ex_map1 <-  SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1,q2,q3)
  ex_map1$mapping()
  ex_map1$tabulate()
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender) + pineapple,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender) +pineapple ,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data. ", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender) +pineapple,
      family = binomial(link="logit")
    ),
    "Not all variables available in the data.", fixed = TRUE
  )

  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = pineapple ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = pineapple ~ (1|age) + (1|gender) ,
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = pineapple ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Outcome variable not present in data. ", fixed = TRUE
  )
})

test_that("Error if not fitting a bernoulli/binomial model", {
  ex_map1 <- SurveyMap$new(samp_obj, popn_obj, q1, q2, q3)
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender)
    ),
    "Currently only binomial and bernoulli models are supported."
  )

})

# this test will fail until https://github.com/mitzimorris/mrp-kit/issues/59
# because the mapped data has already be
test_that("Error if data hasn't been mapped yet",{
  ex_map1 <-  SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2, q3)
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Mapped data not found.", fixed = TRUE
  )
})

test_that("Error if poststrat matrix hasn't been created yet",{
  ex_map1 <- SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2, q3)
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit")
    ),
    "Post-stratification data not found.", fixed = TRUE
  )
})

test_that("Error if data is given as input",{
  ex_map1 <- SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2, q3)
  skip_if_not_installed("lme4")
  expect_error(
    ex_map1$fit(
      fun = lme4::glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
  skip_if_not_installed("rstanarm")
  expect_error(
    ex_map1$fit(
      fun = rstanarm::stan_glmer,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    ex_map1$fit(
      fun = brms::brm,
      formula = y ~ (1|age) + (1|gender),
      family = binomial(link="logit"),
      data = feline_survey
    ),
    "The 'data' argument should not be specified.", fixed = TRUE
  )
})

test_that("Warning is given if fitting using packages that are not lme4, brms, rstanarm ", {
  ex_map1 <- SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2, q3)
  ex_map1$mapping()
  ex_map1$tabulate()
  expect_warning(
    ex_map1$fit(
      fun = stats::glm,
      formula = y ~ age + gender,
      family = "binomial"
    ),
    "Only rstanarm, brms and lme4 are supported natively.", fixed = TRUE
  )
})


test_that("Model fits do not cause errors if specified correctly",{
  ex_map1 <- SurveyMap$new(samp_obj = samp_obj, popn_obj = popn_obj, q1, q2, q3)
  ex_map1$mapping()
  ex_map1$tabulate()
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
    family = "bernoulli",
  )), regexp = NA)

  skip_if_not_installed("lme4")
  expect_error(ex_map$fit(
    fun = lme4::glmer,
    formula = y ~ (1|age) + (1|gender),
    family = "binomial"
  ), regexp = NA)
})

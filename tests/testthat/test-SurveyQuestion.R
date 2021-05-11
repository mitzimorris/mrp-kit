test_that("silent if specified correctly", {
  expect_silent(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    )
  )
})

test_that("error if 'name' specified incorrectly", {
  expect_error(
    SurveyQuestion$new(
      name = c("pet_own","pet_pref"),
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' must be a single string"
  )
  expect_error(
    SurveyQuestion$new(
      name = TRUE,
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' must be a single string"
  )
  expect_error(
    SurveyQuestion$new(
      name = NA_character_,
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' cannot be NA"
  )
})

test_that("error if 'col_names' specified incorrectly", {
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = "pet",
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'col_names' must be a character vector of length 2"
  )
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet", NA_character_),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "NAs not allowed in 'col_names'"
  )
})

test_that("error if 'values_map' specified incorrectly", {
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = NA, "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "NAs not allowed in 'values_map'"
  )
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "All elements of 'values_map' must have names"
  )
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = TRUE
    ),
    "'values_map' must be a list"
  )
})

test_that("warning if duplicated values in 'values_map'", {
  expect_warning(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "cat" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "Duplicated values in map, removing duplicates"
  )
})

test_that("error if not many to one mapping", {
  expect_error(
    SurveyQuestion$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "dog" = "cat", "dog" = "dog")
    ),
    "Package can only handle many to one mappings"
  )
})



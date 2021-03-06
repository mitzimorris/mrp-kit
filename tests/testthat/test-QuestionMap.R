test_that("object has correct R6 class", {
  expect_r6_class(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "QuestionMap"
  )
})

test_that("print output hasn't changed", {
  q <- QuestionMap$new(
    name = "pet",
    col_names = c("pet_own","pet_pref"),
    values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
  )
  expect_known_output(
    print(q),
    file = test_path("answers/QuestionMap-print")
  )
})

test_that("error if 'name' specified incorrectly", {
  expect_error(
    QuestionMap$new(
      name = c("pet_own","pet_pref"),
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' must be a single string"
  )
  expect_error(
    QuestionMap$new(
      name = TRUE,
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' must be a single string"
  )
  expect_error(
    QuestionMap$new(
      name = NA_character_,
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'name' cannot be NA"
  )
})

test_that("error if 'col_names' specified incorrectly", {
  expect_error(
    QuestionMap$new(
      name = "pet",
      col_names = "pet",
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "'col_names' must be a character vector of length 2"
  )
  expect_error(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet", NA_character_),
      values_map = list("cat" = "cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "NAs not allowed in 'col_names'"
  )
})

test_that("error if 'values_map' specified incorrectly", {
  expect_error(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = NA, "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "NAs not allowed in 'values_map'"
  )
  expect_error(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat", "kitten" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "All elements of 'values_map' must have names"
  )
  expect_error(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = TRUE
    ),
    "'values_map' must be a list"
  )
})

test_that("warning if duplicated values in 'values_map'", {
  expect_warning(
    QuestionMap$new(
      name = "pet",
      col_names = c("pet_own","pet_pref"),
      values_map = list("cat" = "cat", "cat" = "cat","dog" = "dog","puppy" = "dog")
    ),
    "Duplicated values in map, removing duplicates"
  )
})


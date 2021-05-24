expect_r6_class <- function(x, class) {
  expect_true(R6::is.R6(x))
  expect_true(inherits(x, class))
}

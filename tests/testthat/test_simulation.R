test_that("rescaled sample sums to 1,000,000", {
  sample = abs(rnorm(10000))
  assert("sums to milltion", round(sum(scale_to_million(sample))) == 1e6)
})
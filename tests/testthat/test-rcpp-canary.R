# Day 15 toolchain canary, lifted to a test.

test_that("RcppArmadillo canary reports correct trace", {
  msg <- cevcmm_rcpp_check()
  expect_type(msg, "character")
  expect_match(msg, "^OK")
  expect_match(msg, "trace.* 3")
})

# Day 21: validation-error tests for the explicit-spec API. Every
# stop() / warning() in vcmm(), vcmm_control(), build_vcmm_design(),
# build_penalty_matrix(), build_kronecker_precision(),
# estimate_kronecker_components(), compute_sufficient_stats(),
# init_accumulator(), accumulate_stats(), and fit_from_summaries()
# gets one test. Cheap to run, big lift for R/control.R, R/vcmm.R,
# R/design.R, R/sufficient_stats.R, and R/distributed.R coverage.

# ---- vcmm_control ----------------------------------------------------------

test_that("vcmm_control rejects bad input", {
  expect_error(vcmm_control(max_iter = -1),         "max_iter")
  expect_error(vcmm_control(max_iter = 1.5),        "max_iter")
  expect_error(vcmm_control(tol_beta = 0),          "tol_beta")
  expect_error(vcmm_control(tol_alpha = -1),        "tol_alpha")
  expect_error(vcmm_control(sigma_eps = -1),        "sigma_eps")
  expect_error(vcmm_control(sigma_alpha = 0),       "sigma_alpha")
  expect_error(vcmm_control(update_variance = NA),  "update_variance")
  expect_error(vcmm_control(verbose = "yes"),       "verbose")
})

# ---- vcmm ------------------------------------------------------------------

test_that("vcmm rejects malformed y/X/Z", {
  set.seed(401)
  N <- 50L
  Z <- matrix(rnorm(N * 2), N, 2)
  t <- runif(N); x <- runif(N)
  y_bad_na <- rnorm(N); y_bad_na[3] <- NA

  expect_error(vcmm(y_bad_na,     X = x, Z = Z, t = t),         "NA")
  expect_error(vcmm("not numeric", X = x, Z = Z, t = t),        "numeric")
  expect_error(vcmm(rnorm(N), X = x, Z = matrix(rnorm(40), 20, 2),
                    t = t),                                     "nrow")
  expect_error(vcmm(rnorm(N), X = x, Z = rnorm(N), t = t),      "matrix")
})

test_that("vcmm requires n_groups for kronecker and q_left for separable", {
  set.seed(402)
  N <- 60L; G <- 4L; q <- 2L * G
  Z <- matrix(rnorm(N * q), N, q)
  t <- runif(N); x <- runif(N); y <- rnorm(N)
  expect_error(vcmm(y, X = x, Z = Z, t = t, re_cov = "kronecker"),
               "n_groups")
  expect_error(vcmm(y, X = x, Z = Z, t = t, re_cov = "separable",
                    n_groups = G),
               "q_left")
})

test_that("vcmm rejects mismatched Z dimensions for kron", {
  set.seed(403)
  N <- 60L; G <- 4L; q <- 17L
  Z <- matrix(rnorm(N * q), N, q)
  t <- runif(N); x <- runif(N); y <- rnorm(N)
  expect_error(
    vcmm(y, X = x, Z = Z, t = t, re_cov = "kronecker", n_groups = G),
    "q_left \\* n_groups")
})

test_that("vcmm rejects multiple aliases for Sigma_left", {
  set.seed(404)
  N <- 60L; G <- 4L; q <- 2L * G
  Z <- matrix(rnorm(N * q), N, q)
  t <- runif(N); x <- runif(N); y <- rnorm(N)
  expect_error(
    vcmm(y, X = x, Z = Z, t = t, re_cov = "kronecker", n_groups = G,
         Sigma_2x2_init = diag(2), Sigma_left_init = diag(2)),
    "Pass only ONE")
})

test_that("vcmm rejects mis-shaped Sigma_left_init / Sigma_right_init", {
  set.seed(405)
  N <- 60L; G <- 4L; q <- 2L * G
  Z <- matrix(rnorm(N * q), N, q)
  t <- runif(N); x <- runif(N); y <- rnorm(N)
  expect_error(
    vcmm(y, X = x, Z = Z, t = t, re_cov = "kronecker", n_groups = G,
         Sigma_2x2_init = diag(3)),
    "2 by 2")
  expect_error(
    vcmm(y, X = x, Z = Z, t = t, re_cov = "kronecker", n_groups = G,
         Sigma_spatial_init = diag(3)),
    sprintf("%d by %d", G, G))
})

# ---- build_vcmm_design -----------------------------------------------------

test_that("build_vcmm_design validates inputs", {
  expect_error(build_vcmm_design(X = "bad", t = 1:5),         "matrix or vector")
  expect_error(build_vcmm_design(X = matrix(c(1, NA), 2, 1),
                                  t = c(0, 1)),               "NA")
  expect_error(build_vcmm_design(X = 1:10, t = 1:5),          "length")
  expect_error(build_vcmm_design(X = 1:10, t = 1:10, degree = 0),
               "degree")
  expect_error(build_vcmm_design(X = 1:10, t = 1:10, n_basis = 2),
               "degree \\+ 1")
})

test_that("build_vcmm_design errors on zero-range t", {
  expect_error(build_vcmm_design(X = 1:10, t = rep(0.5, 10)),
               "zero range")
})

test_that("build_vcmm_design warns when normalize_t = FALSE and t out of range", {
  set.seed(406)
  N <- 50L
  t <- runif(N, 0, 2)   # values > 1
  x <- runif(N)
  expect_warning(
    build_vcmm_design(X = x, t = t, normalize_t = FALSE),
    "extrapolate")
})

# ---- build_penalty_matrix --------------------------------------------------

test_that("build_penalty_matrix validates inputs", {
  expect_error(build_penalty_matrix(n_basis = 2, lambda = 1),  ">= 3")
  expect_error(build_penalty_matrix(n_basis = 5, lambda = -1), "non-negative")
  expect_error(build_penalty_matrix(n_basis = 5, lambda = 1, n_blocks = 0),
               ">= 1")
})

# ---- build_kronecker_precision + estimate_kronecker_components -------------

test_that("build_kronecker_precision rejects bad inputs", {
  expect_error(build_kronecker_precision(
    Sigma_left = matrix(1:6, 2, 3), Sigma_right = diag(3),
    sigma_eps = 0.5),
    "square")
  expect_error(build_kronecker_precision(
    Sigma_left = diag(2), Sigma_right = matrix(1:6, 2, 3),
    sigma_eps = 0.5),
    "square")
  expect_error(build_kronecker_precision(
    Sigma_left = diag(2), Sigma_right = diag(3), sigma_eps = -1),
    "positive")
})

test_that("estimate_kronecker_components validates input", {
  expect_error(estimate_kronecker_components(alpha = "bad", n_groups = 4),
               "numeric")
  expect_error(estimate_kronecker_components(alpha = rnorm(8), n_groups = -1),
               "positive integer")
  expect_error(estimate_kronecker_components(alpha = rnorm(7),
                                              n_groups = 4, q_left = 2),
               "does not match")
})

# ---- compute_sufficient_stats / init_accumulator / accumulate_stats --------

test_that("compute_sufficient_stats validates inputs", {
  expect_error(compute_sufficient_stats("y", matrix(1, 2, 1),
                                         matrix(1, 2, 1)),
               "numeric")
  expect_error(compute_sufficient_stats(1:5, matrix(1, 4, 1),
                                         matrix(1, 5, 1)),
               "nrow")
  expect_error(compute_sufficient_stats(c(1, NA, 3),
                                         matrix(1, 3, 1),
                                         matrix(1, 3, 1)),
               "NA")
})

test_that("init_accumulator validates p, q", {
  expect_error(init_accumulator(p =  0, q = 3), "positive integer")
  expect_error(init_accumulator(p =  3, q = -1), "positive integer")
})

test_that("accumulate_stats catches class and dim mismatches", {
  acc <- init_accumulator(p = 3, q = 2)
  ss_wrong <- compute_sufficient_stats(1:4, matrix(1, 4, 5),
                                        matrix(1, 4, 2))
  expect_error(accumulate_stats(acc, ss_wrong),  "Dimension mismatch")
  expect_error(accumulate_stats("bad", ss_wrong), "vcmm_accumulator")
  expect_error(accumulate_stats(acc, "bad"),      "vcmm_ss")
})

# ---- fit_from_summaries ----------------------------------------------------

test_that("fit_from_summaries rejects malformed inputs", {
  d <- gen_diag_data(N = 60L, q = 2L)
  design <- build_vcmm_design(X = d$x, t = d$t)
  ss <- compute_sufficient_stats(d$y, design$X_design, d$Z)

  expect_error(fit_from_summaries(list(), penalty = diag(3)),
               "empty list")
  expect_error(fit_from_summaries(list("not a summary"),
                                   penalty = diag(3)),
               "vcmm_ss")
  expect_error(fit_from_summaries(42, penalty = diag(3)),
               "vcmm_ss")
  expect_error(fit_from_summaries(ss, penalty = diag(3)),
               "by")
  expect_error(fit_from_summaries(ss, penalty = design$penalty,
                                   re_cov = "kronecker"),
               "n_groups")
  expect_error(fit_from_summaries(ss, penalty = design$penalty,
                                   re_cov = "separable", n_groups = 1),
               "q_left")
  expect_error(fit_from_summaries(ss, penalty = design$penalty,
                                   method = "ss",
                                   rowsum_constant = "not numeric"),
               "rowsum_constant")
})

# ---- +.vcmm_ss dimension mismatch ------------------------------------------

test_that("+.vcmm_ss rejects mismatched operands", {
  ss_a <- compute_sufficient_stats(rnorm(20), matrix(rnorm(20), 20, 1),
                                    matrix(rnorm(40), 20, 2))
  ss_b <- compute_sufficient_stats(rnorm(20), matrix(rnorm(40), 20, 2),
                                    matrix(rnorm(40), 20, 2))
  expect_error(ss_a + ss_b, "Dimension mismatch")
  expect_error(ss_a + 1,    "vcmm_ss")
})

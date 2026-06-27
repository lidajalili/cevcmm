#===============================================================================
# Day 13: plot.vcmm_fit validation
#
# Strategy
# --------
# Plotting is inherently visual; we cannot programmatically "assert" a plot
# looks right. So this script does two things:
#
#   (a) FUNCTIONAL check: call plot.vcmm_fit() with every (which, re_cov)
#       combination, plus the "missing data" path for panel 2, inside
#       tryCatch() so any error fails the run loudly. We test that
#       nothing throws.
#
#   (b) VISUAL artefact: render all panels to a single PDF at
#       inst/validation/day13_plots.pdf. The user can open it and
#       confirm the panels look sensible: CI bands hug a sin curve, residuals
#       are mean-zero around 0, ranef heatmap looks like the truth, etc.
#
# We also assert that the PDF file gets created and is non-trivially large.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

.expect_true <- function(cond, label) {
  if (!isTRUE(cond)) stop(sprintf("[FAIL] %s", label), call. = FALSE)
  cat(sprintf("  [PASS] %s\n", label))
}
.expect_no_error <- function(expr, label) {
  result <- tryCatch(expr, error = function(e) e)
  if (inherits(result, "error")) {
    stop(sprintf("[FAIL] %s -- threw: %s", label,
                 conditionMessage(result)),
         call. = FALSE)
  }
  cat(sprintf("  [PASS] %s\n", label))
}

# ----- Shared simulation harness (matches earlier validation scripts) ----
make_data <- function(N, G, q_left, re_cov_type, sigma_eps = 0.5, seed = 13L) {
  set.seed(seed)
  q <- if (re_cov_type == "diag") max(5L, G) else q_left * G
  t_vec <- runif(N); x_vec <- runif(N)
  beta0_true <- 2.0
  beta1_fn   <- function(tt) sin(2 * pi * tt)

  if (re_cov_type == "kronecker") {
    origin_id <- sample.int(G, N, replace = TRUE)
    dest_id   <- sample.int(G, N, replace = TRUE)
    Z <- matrix(0, N, q)
    Z[cbind(seq_len(N), origin_id)]   <- 1
    Z[cbind(seq_len(N), G + dest_id)] <- 1
    Sigma_left  <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
    Sigma_right <- outer(seq_len(G), seq_len(G),
                         function(i, j) exp(-abs(i - j) / 5))
    alpha_true  <- as.vector(crossprod(chol(kronecker(Sigma_left, Sigma_right)),
                                       rnorm(q)))
  } else if (re_cov_type == "separable") {
    g_id <- sample.int(G, N, replace = TRUE)
    W    <- matrix(rnorm(N * q_left), N, q_left)
    Z    <- matrix(0, N, q)
    for (kk in seq_len(q_left)) {
      Z[cbind(seq_len(N), (kk - 1L) * G + g_id)] <- W[, kk]
    }
    Sigma_left  <- outer(seq_len(q_left), seq_len(q_left),
                         function(i, j) 0.5 * 0.3^abs(i - j))
    Sigma_right <- outer(seq_len(G), seq_len(G),
                         function(i, j) exp(-abs(i - j) / 6))
    alpha_true  <- as.vector(crossprod(chol(kronecker(Sigma_left, Sigma_right)),
                                       rnorm(q)))
  } else { # diag
    Z <- matrix(rnorm(N * q), N, q)
    alpha_true <- rnorm(q, sd = 0.4)
    Sigma_left <- Sigma_right <- NULL
  }

  y <- beta0_true + beta1_fn(t_vec) * x_vec +
       as.vector(Z %*% alpha_true) + rnorm(N, sd = sigma_eps)

  list(y = y, x = x_vec, Z = Z, t = t_vec, q = q,
       Sigma_left = Sigma_left, Sigma_right = Sigma_right,
       sigma_eps_true = sigma_eps)
}

# ----- Fit one model per re_cov mode -------------------------------------
cat("[Day 13] Fitting three models (diag, kronecker, separable)...\n")

d <- make_data(N = 500L,  G = 5L,  q_left = 1L, re_cov_type = "diag")
fit_d <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t, method = "csl",
              re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                                     update_variance = TRUE))

k <- make_data(N = 2000L, G = 20L, q_left = 2L, re_cov_type = "kronecker")
fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t, method = "csl",
              re_cov = "kronecker", n_groups = 20L,
              Sigma_spatial_init = k$Sigma_right,
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

s <- make_data(N = 3000L, G = 15L, q_left = 4L, re_cov_type = "separable")
fit_s <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t, method = "csl",
              re_cov = "separable", n_groups = 15L, q_left = 4L,
              Omega_G_init = s$Sigma_right,
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

# ----- Functional checks: every (which, re_cov) combination ---------------
cat("\n========== Functional checks (no errors thrown) ==========\n")

# Suppress on-screen graphical output in non-interactive contexts.
pdf_path <- "inst/validation/day13_plots.pdf"
if (!dir.exists(dirname(pdf_path))) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
}

grDevices::pdf(pdf_path, width = 10, height = 5)

for (mode in c("diag", "kron", "sep")) {
  fit  <- switch(mode, diag = fit_d, kron = fit_k, sep = fit_s)
  dlst <- switch(mode,
                 diag = list(y = d$y, X = d$x, Z = d$Z, t = d$t),
                 kron = list(y = k$y, X = k$x, Z = k$Z, t = k$t),
                 sep  = list(y = s$y, X = s$x, Z = s$Z, t = s$t))

  graphics::par(oma = c(0, 0, 2, 0))   # outer margin for the mode title

  # Panel 1
  .expect_no_error(
    {plot(fit, which = 1L, ask = FALSE);
     graphics::mtext(sprintf("[%s] panel 1: varying coefficient(s)", mode),
                     side = 3, line = 0.5, outer = TRUE, cex = 1.1)},
    sprintf("plot(fit, which = 1) works for re_cov = '%s'", mode))

  # Panel 2 -- requires data
  .expect_no_error(
    {plot(fit, which = 2L, data = dlst, ask = FALSE);
     graphics::mtext(sprintf("[%s] panel 2: residual diagnostics", mode),
                     side = 3, line = 0.5, outer = TRUE, cex = 1.1)},
    sprintf("plot(fit, which = 2, data = ...) works for re_cov = '%s'", mode))

  # Panel 3
  .expect_no_error(
    {plot(fit, which = 3L, ask = FALSE);
     graphics::mtext(sprintf("[%s] panel 3: random-effects diagnostics", mode),
                     side = 3, line = 0.5, outer = TRUE, cex = 1.1)},
    sprintf("plot(fit, which = 3) works for re_cov = '%s'", mode))
}

grDevices::dev.off()

# ----- Edge cases --------------------------------------------------------
cat("\n========== Edge cases ==========\n")

# (a) Missing data for panel 2 -> warning, not error
.expect_no_error(
  withCallingHandlers(
    plot(fit_k, which = 2L, ask = FALSE),
    warning = function(w) invokeRestart("muffleWarning")
  ),
  "plot(fit, which = 2) without data warns and continues")

# (b) Invalid which -> hard error
result_bad <- tryCatch(plot(fit_k, which = 4L, ask = FALSE),
                       error = function(e) e)
.expect_true(inherits(result_bad, "error"),
             "plot(fit, which = 4) throws an error")

# (c) varying_coef with se.fit returns list of right shape
vc_se <- varying_coef(fit_k, t_new = seq(0, 1, length.out = 11L),
                      se.fit = TRUE)
.expect_true(is.list(vc_se) &&
             identical(dim(vc_se$fit),    c(11L, 1L)) &&
             identical(dim(vc_se$se.fit), c(11L, 1L)) &&
             all(vc_se$se.fit >= 0),
             "varying_coef(..., se.fit = TRUE) returns list(fit, se.fit) of right shape")

# (d) The PDF was written and is non-trivially large
.expect_true(file.exists(pdf_path),
             sprintf("PDF exists at %s", pdf_path))
.expect_true(file.size(pdf_path) > 5000,
             "PDF size > 5 KB (non-empty)")

cat(sprintf("\n[Day 13] Wrote diagnostic plots to %s (%.1f KB).\n",
            pdf_path, file.size(pdf_path) / 1024))
cat("Open it to visually verify the panels look correct.\n")
cat("\nAll plot functional checks passed.\n")
invisible(list(fit_d = fit_d, fit_k = fit_k, fit_s = fit_s))

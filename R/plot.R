#===============================================================================
# plot.vcmm_fit -- three diagnostic panels in the plot.lm style
#
# panel 1: each varying coefficient hat_beta_k(t) with pointwise CI band.
# panel 2: residual diagnostics (residuals vs fitted, residuals vs t).
#          Requires `data` because the fit object does NOT carry the raw
#          training (y, X, Z, t) -- the package is communication-efficient
#          and only stores sufficient statistics.
# panel 3: random-effects diagnostics:
#            - re_cov = "diag":      Normal Q-Q plot of hat_alpha.
#            - re_cov = "kronecker" or "separable":
#                heatmap of the (G x q_left) ranef matrix +
#                heatmap of the estimated Sigma_left.
#
# Usage mimics plot.lm:
#   plot(fit)                         # all three panels in order
#   plot(fit, which = 1)              # just varying coefficients
#   plot(fit, which = 2, data = ...)  # just residuals
#   plot(fit, which = 3)              # just random effects
#
# No ggplot2 dependency: base R graphics only.
#===============================================================================

#' Diagnostic plots for a vcmm fit
#'
#' Three panels, each requestable via \code{which}:
#' \describe{
#'   \item{1}{Varying-coefficient curves
#'     \eqn{\hat\beta_k(t)} with pointwise Wald confidence bands.}
#'   \item{2}{Residual diagnostics: residuals vs fitted, residuals vs
#'     \eqn{t}. \strong{Requires \code{data}} since a vcmm fit does not
#'     carry the raw training data.}
#'   \item{3}{Random-effects diagnostics: Normal Q-Q for
#'     \code{re_cov = "diag"}, or a heatmap of the
#'     \eqn{G \times q_{\mathrm{left}}} random-effects matrix together
#'     with a heatmap of the estimated \eqn{\Sigma_{\mathrm{left}}} for
#'     \code{"kronecker"} / \code{"separable"}.}
#' }
#'
#' Uses base R graphics, no \pkg{ggplot2} dependency. Each requested
#' panel is drawn on its own figure (or set of subfigures via
#' \code{par(mfrow)}); call \code{par(mfrow = c(2, 2))} or similar
#' before \code{plot()} to combine panels in one figure.
#'
#' @param x A \code{vcmm_fit} object.
#' @param which Integer vector subset of \code{1:3}. Default
#'   \code{1:3}.
#' @param data Optional list with components \code{y}, \code{X},
#'   \code{Z}, \code{t} (the training data, or any data on which to
#'   compute residuals). Required when panel 2 is requested.
#' @param t_grid Numeric. Grid of \eqn{t} values for panel 1. Default
#'   \code{NULL} means an evenly-spaced grid over the stored training
#'   range \code{[t_min, t_max]}.
#' @param n_grid Integer. Number of grid points if \code{t_grid} is
#'   \code{NULL}. Default 200.
#' @param conf_level Numeric in (0, 1). Confidence level for panel-1
#'   bands. Default 0.95.
#' @param ask Logical. Passed to \code{devAskNewPage()} when multiple
#'   panels are requested in an interactive session.
#' @param ... Further arguments passed to base graphics calls.
#'
#' @return Invisibly \code{NULL}. Called for side effects (plots).
#'
#' @seealso \code{\link{varying_coef}}, \code{\link{ranef.vcmm_fit}},
#'   \code{\link{predict.vcmm_fit}}.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
#' a <- rnorm(3, sd = 0.5)
#' y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% a) + rnorm(n, sd = 0.5)
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' \dontrun{
#' plot(fit)                                         # all three panels
#' plot(fit, which = 1)                              # only varying coefs
#' plot(fit, which = 2, data = list(y = y, X = x, Z = Z, t = t))
#' plot(fit, which = 3)                              # ranef diagnostics
#' }
plot.vcmm_fit <- function(x,
                          which      = 1:3,
                          data       = NULL,
                          t_grid     = NULL,
                          n_grid     = 200L,
                          conf_level = 0.95,
                          ask        = (length(which) > 1L) && interactive(),
                          ...) {
  if (!inherits(x, "vcmm_fit")) {
    stop("`x` must be a 'vcmm_fit'.", call. = FALSE)
  }
  show <- intersect(as.integer(which), 1:3)
  if (length(show) == 0L) {
    stop("`which` must contain at least one of 1, 2, 3.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric in (0, 1).", call. = FALSE)
  }

  if (isTRUE(ask)) {
    old_ask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)
  }

  for (k in show) {
    if (identical(k, 1L)) {
      .plot_varying_coef(x, t_grid = t_grid, n_grid = n_grid,
                         conf_level = conf_level, ...)
    } else if (identical(k, 2L)) {
      if (is.null(data)) {
        warning("plot.vcmm_fit: panel 2 requires `data` (a list with y, X, Z, t). Skipping.",
                call. = FALSE)
      } else {
        .plot_residuals(x, data = data, ...)
      }
    } else if (identical(k, 3L)) {
      .plot_ranef(x, ...)
    }
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------
# Panel 1: varying-coefficient curves with CI bands.
# -----------------------------------------------------------------------
.plot_varying_coef <- function(fit, t_grid, n_grid, conf_level, ...) {
  ds <- fit$design
  if (is.null(ds)) {
    stop("Fit has no design metadata; was it produced by vcmm()?",
         call. = FALSE)
  }
  K <- as.integer(ds$K)
  if (K == 0L) {
    message("No varying coefficients to plot.")
    return(invisible(NULL))
  }

  if (is.null(t_grid)) {
    t_grid <- seq(ds$t_min, ds$t_max, length.out = as.integer(n_grid))
  } else {
    if (!is.numeric(t_grid)) {
      stop("`t_grid` must be numeric.", call. = FALSE)
    }
    t_grid <- sort(unique(as.numeric(t_grid)))
  }

  vc <- varying_coef(fit, t_new = t_grid, se.fit = TRUE)
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  lower  <- vc$fit - z_crit * vc$se.fit
  upper  <- vc$fit + z_crit * vc$se.fit

  par_old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(par_old), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(K))

  for (k in seq_len(K)) {
    yl <- range(c(lower[, k], upper[, k], 0), finite = TRUE)
    plot(t_grid, vc$fit[, k], type = "n", ylim = yl,
         xlab = "t",
         ylab = bquote(hat(beta)[.(k)] * (t)),
         main = sprintf("Varying coefficient beta_%d(t) (%.0f%% band)",
                        k, 100 * conf_level), ...)
    graphics::polygon(c(t_grid, rev(t_grid)),
                      c(upper[, k], rev(lower[, k])),
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.25),
                      border = NA)
    graphics::lines(t_grid, vc$fit[, k], col = "steelblue", lwd = 2)
    graphics::abline(h = 0, lty = 3, col = "gray60")
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------
# Panel 2: residual diagnostics (residuals vs fitted, residuals vs t).
# Requires `data` because vcmm_fit does not carry the raw training data.
# -----------------------------------------------------------------------
.plot_residuals <- function(fit, data, ...) {
  if (!is.list(data)) {
    stop("`data` must be a named list with components y, X, Z, t.",
         call. = FALSE)
  }
  required <- c("y", "t")
  missing_fields <- setdiff(required, names(data))
  if (length(missing_fields) > 0L) {
    stop(sprintf("`data` is missing: %s.",
                 paste(missing_fields, collapse = ", ")),
         call. = FALSE)
  }

  y    <- as.numeric(data$y)
  yhat <- predict.vcmm_fit(fit, newdata = data, include_random = TRUE)
  resid <- y - yhat

  par_old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(par_old), add = TRUE)
  graphics::par(mfrow = c(1L, 2L))

  # Residuals vs fitted
  plot(yhat, resid,
       xlab = "Fitted",
       ylab = "Residual",
       main = "Residuals vs Fitted",
       pch = 16, cex = 0.6,
       col = grDevices::adjustcolor("black", alpha.f = 0.4),
       ...)
  graphics::abline(h = 0, lty = 2, col = "red")
  graphics::lines(stats::lowess(yhat, resid), col = "red", lwd = 1.5)

  # Residuals vs t
  plot(data$t, resid,
       xlab = "t",
       ylab = "Residual",
       main = "Residuals vs t",
       pch = 16, cex = 0.6,
       col = grDevices::adjustcolor("black", alpha.f = 0.4),
       ...)
  graphics::abline(h = 0, lty = 2, col = "red")
  graphics::lines(stats::lowess(data$t, resid), col = "red", lwd = 1.5)
  invisible(NULL)
}

# -----------------------------------------------------------------------
# Panel 3: random-effects diagnostics.
# - diag      : Normal Q-Q of hat_alpha.
# - kron / sep: heatmaps of the (G x q_left) ranef matrix and Sigma_left.
# -----------------------------------------------------------------------
.plot_ranef <- function(fit, ...) {
  state <- fit$re_cov_state
  re    <- ranef.vcmm_fit(fit)

  if (is.null(state) || identical(state$type, "diag")) {
    par_old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(par_old), add = TRUE)
    graphics::par(mfrow = c(1L, 1L))

    stats::qqnorm(as.numeric(re),
                  main = "Random effects: Normal Q-Q",
                  pch = 16, cex = 0.7,
                  col = grDevices::adjustcolor("black", alpha.f = 0.6), ...)
    stats::qqline(as.numeric(re), col = "red", lty = 2)
    return(invisible(NULL))
  }

  # kronecker or separable: re is G x q_left
  G  <- nrow(re)
  kL <- ncol(re)
  SL <- state$Sigma_left

  par_old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(par_old), add = TRUE)
  graphics::par(mfrow = c(1L, 2L))

  pal <- grDevices::colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(64)

  # Left panel: ranef heatmap
  rng_re <- max(abs(re), na.rm = TRUE)
  if (rng_re == 0) rng_re <- 1
  graphics::image(seq_len(G), seq_len(kL),
                  matrix(re, nrow = G, ncol = kL),
                  col = pal,
                  zlim = c(-rng_re, rng_re),
                  xlab = "Group g",
                  ylab = if (identical(state$type, "kronecker") &&
                             identical(as.integer(kL), 2L))
                           "O/D" else "Effect type k",
                  main = sprintf("Random effects (%d x %d)", G, kL),
                  axes = TRUE, ...)
  if (identical(state$type, "kronecker") &&
      identical(as.integer(kL), 2L)) {
    graphics::axis(2, at = 1:2, labels = c("origin", "dest"))
  }

  # Right panel: Sigma_left heatmap
  rng_S <- max(abs(SL), na.rm = TRUE)
  if (rng_S == 0) rng_S <- 1
  left_lbl <- if (identical(state$type, "separable")) "Sigma_q" else
              if (identical(as.integer(kL), 2L))      "Sigma_2x2" else
                                                       "Sigma_left"
  graphics::image(seq_len(kL), seq_len(kL),
                  matrix(SL, nrow = kL, ncol = kL),
                  col = pal,
                  zlim = c(-rng_S, rng_S),
                  xlab = "k", ylab = "l",
                  main = sprintf("Estimated %s (%d x %d)",
                                 left_lbl, kL, kL),
                  axes = TRUE, ...)
  # Overlay numeric values for small matrices
  if (kL <= 6L) {
    for (i in seq_len(kL)) for (j in seq_len(kL)) {
      graphics::text(i, j, sprintf("%.2f", SL[i, j]), cex = 0.8)
    }
  }
  invisible(NULL)
}

#===============================================================================
# vcmm: main user-facing fit function for the Communication-Efficient
# Varying Coefficient Mixed-Effects Model
#
# Wraps:
#   build_vcmm_design()        -> X_design, penalty
#   compute_sufficient_stats() -> aggregated SS
#   fit_ss() or fit_csl()      -> the chosen estimator
#
# Three re_cov modes supported:
#   "diag"      -- alpha ~ N(0, sigma_alpha^2 I_q)        (simplest)
#   "kronecker" -- alpha ~ N(0, Sigma_left ⊗ Sigma_right) (OD-style)
#   "separable" -- alpha ~ N(0, Sigma_q ⊗ Omega_G)        (group-shared dense)
#
# "kronecker" and "separable" share the same internal Kronecker machinery
# (in covariance.R) with q_left = 2 the default for OD-style and required
# (no default) for separable. See ?vcmm for the column-stacking convention
# expected of Z and the initial covariance matrices.
#===============================================================================

#' Fit a Varying Coefficient Mixed-Effects Model
#'
#' The main user-facing fit function. Builds the B-spline design and
#' penalty for one or more varying coefficients in \code{t}, computes
#' the aggregated sufficient statistics, and fits the model using either
#' the one-step CSL estimator (default) or the iterative SS estimator.
#'
#' \strong{Model.} For observations \eqn{i = 1, \ldots, n} the fitted
#' model is
#' \deqn{
#'   y_i
#'   = \beta_0(t_i) + \sum_{k=1}^{K} x_{ik}\, \beta_k(t_i)
#'   + z_i^\top \alpha + \varepsilon_i,
#' }
#' where each \eqn{\beta_k(t)} is a cubic B-spline with \code{n_basis}
#' basis functions and a second-order difference penalty, and
#' \eqn{\alpha \sim N(0, \Sigma_\alpha)} with structure chosen by
#' \code{re_cov}.
#'
#' \strong{Method selection.} The default is \code{method = "auto"},
#' which picks \code{"csl"} when \eqn{N \cdot q > 10^5} or \eqn{q > 50}
#' and \code{"ss"} otherwise.
#'
#' \strong{Random-effects covariance.} Three structures:
#' \itemize{
#'   \item \code{re_cov = "diag"}: \eqn{\alpha \sim N(0,
#'     \sigma_\alpha^2 I_q)}.
#'   \item \code{re_cov = "kronecker"}: \eqn{\alpha \sim N(0,
#'     \Sigma_{\mathrm{left}} \otimes \Sigma_{\mathrm{right}})} with
#'     \eqn{\Sigma_{\mathrm{left}}} of size \code{q_left x q_left}
#'     (default \code{q_left = 2}; OD-style with origin / destination
#'     blocks). User-facing names \code{Sigma_2x2_init},
#'     \code{Sigma_spatial_init} are accepted as aliases.
#'   \item \code{re_cov = "separable"}: \eqn{\alpha \sim N(0,
#'     \Sigma_q \otimes \Omega_G)} with \code{Sigma_q} of size
#'     \code{q_left x q_left} (required, no default) and \code{Omega_G}
#'     of size \eqn{G \times G}. User-facing names \code{Sigma_q_init},
#'     \code{Omega_G_init} are accepted as aliases for
#'     \code{Sigma_left_init}, \code{Sigma_right_init}.
#' }
#'
#' \strong{Column-stacking convention.} For \code{re_cov = "kronecker"}
#' or \code{"separable"}, the random-effects vector is
#' \eqn{\alpha = \mathrm{vec}_{\mathrm{col}}(M)} where
#' \eqn{M \in \mathbb R^{G \times q_{\mathrm{left}}}}, i.e.\
#' \code{alpha[(k - 1) * G + g] = M[g, k]}. The \code{Z} matrix must be
#' constructed so that \code{Z \%*\% alpha} gives the correct random-effect
#' contribution. \code{ncol(Z)} must equal \code{q_left * n_groups}.
#'
#' \strong{Identifiability of the right component.} The right-side
#' covariance (\code{Sigma_spatial} / \code{Omega_G}) is not separately
#' identifiable from a single \eqn{\hat\alpha}, so it is held fixed at
#' the user-supplied initial value (default \eqn{I_G}). Supply a
#' parametric kernel via \code{Sigma_right_init} (e.g., \code{exp(-D /
#' phi)} for OD; AR(1) for separable). The left-side covariance
#' (\code{Sigma_2x2} / \code{Sigma_q}) is updated every iteration via the
#' EM-style estimator (Theorem 1 \eqn{M_\eta} rule) when
#' \code{control$update_variance = TRUE}.
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param X Numeric \eqn{n \times K} matrix (or length-\eqn{n} vector if
#'   \eqn{K = 1}) of covariates that get varying coefficients in
#'   \code{t}.
#' @param Z Numeric \eqn{n \times q} random-effects design matrix.
#'   For \code{re_cov = "kronecker"} or \code{"separable"}, \eqn{q} must
#'   equal \code{q_left * n_groups}.
#' @param t Numeric vector of length \eqn{n} in which the coefficients
#'   vary smoothly.
#' @param method Character: \code{"auto"} (default), \code{"csl"}, or
#'   \code{"ss"}.
#' @param re_cov Character: \code{"diag"} (default), \code{"kronecker"},
#'   or \code{"separable"}. See Details.
#' @param n_groups Integer \eqn{G}. Required for \code{"kronecker"} and
#'   \code{"separable"}.
#' @param q_left Integer. The left (within) dimension of the Kronecker
#'   factor. For \code{re_cov = "kronecker"} defaults to 2 (OD setting);
#'   for \code{re_cov = "separable"} this is the per-group random-effect
#'   dimension and \strong{must be supplied}.
#' @param Sigma_left_init Optional \eqn{k \times k} initial left
#'   covariance (\eqn{k = q_{\mathrm{left}}}). Aliases accepted:
#'   \code{Sigma_2x2_init} (for \code{"kronecker"} with \code{q_left =
#'   2}) and \code{Sigma_q_init} (for \code{"separable"}). Default is
#'   \code{sigma_alpha^2 * I_k}.
#' @param Sigma_right_init Optional \eqn{G \times G} initial right
#'   covariance. Aliases accepted: \code{Sigma_spatial_init} (for
#'   \code{"kronecker"}) and \code{Omega_G_init} (for \code{"separable"}).
#'   Default is \eqn{I_G}.
#' @param Sigma_2x2_init Legacy alias for \code{Sigma_left_init} when
#'   \code{re_cov = "kronecker"} and \code{q_left = 2}.
#' @param Sigma_spatial_init Legacy alias for \code{Sigma_right_init}
#'   when \code{re_cov = "kronecker"}.
#' @param Sigma_q_init Alias for \code{Sigma_left_init} when
#'   \code{re_cov = "separable"}.
#' @param Omega_G_init Alias for \code{Sigma_right_init} when
#'   \code{re_cov = "separable"}.
#' @param n_basis Integer or \code{NULL}. Number of B-spline basis
#'   functions per varying coefficient. \code{NULL} (default) auto-picks
#'   \code{max(floor(n^(1/3)) + 4, 10)}.
#' @param degree Integer. B-spline degree (default 3 = cubic).
#' @param lambda Non-negative numeric. Smoothing parameter (default 1).
#' @param control A \code{vcmm_control} object.
#' @param normalize_t Logical. If \code{TRUE} (default), \code{t} is
#'   linearly mapped to \code{[0, 1]} before building the basis.
#' @param ... Further arguments passed to \code{fit_csl()} when CSL is
#'   used.
#'
#' @return A \code{vcmm_fit} object as returned by \code{fit_ss()} or
#'   \code{fit_csl()}, augmented with \code{design} (basis metadata) and
#'   \code{re_cov}. When \code{re_cov} is \code{"kronecker"} or
#'   \code{"separable"}, \code{re_cov_state} contains the canonical
#'   fields \code{Sigma_left}, \code{Sigma_right}, plus legacy aliases
#'   \code{Sigma_2x2}/\code{Sigma_spatial} (when \code{q_left = 2}) or
#'   \code{Sigma_q}/\code{Omega_G} (for separable).
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' t  <- runif(n)
#' x  <- runif(n)
#' Z  <- matrix(rnorm(n * 3), n, 3)
#' alpha_true <- rnorm(3, sd = 0.5)
#' y  <- 2 + sin(2 * pi * t) * x +
#'       as.vector(Z %*% alpha_true) + rnorm(n, sd = 0.5)
#'
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' fit
vcmm <- function(y,
                 X,
                 Z,
                 t,
                 method             = c("auto", "csl", "ss"),
                 re_cov             = c("diag", "kronecker", "separable"),
                 n_groups           = NULL,
                 q_left             = NULL,
                 Sigma_left_init    = NULL,
                 Sigma_right_init   = NULL,
                 Sigma_2x2_init     = NULL,
                 Sigma_spatial_init = NULL,
                 Sigma_q_init       = NULL,
                 Omega_G_init       = NULL,
                 n_basis            = NULL,
                 degree             = 3L,
                 lambda             = 1,
                 control            = vcmm_control(),
                 normalize_t        = TRUE,
                 ...) {

  # ----- Argument matching and basic validation ---------------------------
  method <- match.arg(method)
  re_cov <- match.arg(re_cov)

  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("NA values are not allowed in `y`.", call. = FALSE)
  }
  n <- length(y)

  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(Z) != n) {
    stop(sprintf("nrow(Z) = %d does not match length(y) = %d.",
                 nrow(Z), n),
         call. = FALSE)
  }
  if (anyNA(Z)) {
    stop("NA values are not allowed in `Z`.", call. = FALSE)
  }
  q_dim <- ncol(Z)

  # ----- re_cov-specific validation and initial state --------------------
  re_cov_state <- NULL

  if (re_cov == "kronecker" || re_cov == "separable") {
    # n_groups required for both
    if (is.null(n_groups) ||
        !is.numeric(n_groups) || length(n_groups) != 1L ||
        n_groups < 1L || n_groups != as.integer(n_groups)) {
      stop(sprintf(
        "re_cov = '%s' requires `n_groups` (a single positive integer).",
        re_cov),
        call. = FALSE)
    }
    n_groups <- as.integer(n_groups)

    # Resolve q_left
    if (re_cov == "kronecker") {
      if (is.null(q_left)) {
        q_left <- 2L                       # OD-style default
      }
    } else { # separable
      if (is.null(q_left)) {
        stop("re_cov = 'separable' requires `q_left` (per-group random-effect dimension).",
             call. = FALSE)
      }
    }
    if (!is.numeric(q_left) || length(q_left) != 1L ||
        q_left < 1L || q_left != as.integer(q_left)) {
      stop("`q_left` must be a single positive integer.", call. = FALSE)
    }
    q_left <- as.integer(q_left)

    # Shape check on Z
    if (q_dim != q_left * n_groups) {
      stop(sprintf(
        "For re_cov = '%s', ncol(Z) must equal q_left * n_groups = %d * %d = %d; got ncol(Z) = %d.",
        re_cov, q_left, n_groups, q_left * n_groups, q_dim),
        call. = FALSE)
    }

    # ----- Resolve Sigma_left_init from canonical + aliases --------------
    Sigma_left_aliases  <- list(
      Sigma_left_init = Sigma_left_init,
      Sigma_2x2_init  = Sigma_2x2_init,
      Sigma_q_init    = Sigma_q_init
    )
    supplied_left <- names(Sigma_left_aliases)[!vapply(
      Sigma_left_aliases, is.null, logical(1L))]
    if (length(supplied_left) > 1L) {
      stop(sprintf(
        "Pass only ONE of: %s.", paste(supplied_left, collapse = ", ")),
        call. = FALSE)
    }
    if (length(supplied_left) == 1L) {
      Sigma_left_init <- Sigma_left_aliases[[supplied_left]]
    }
    if (is.null(Sigma_left_init)) {
      Sigma_left_init <- (control$sigma_alpha^2) * diag(q_left)
    }
    if (!is.matrix(Sigma_left_init) ||
        !isTRUE(all.equal(dim(Sigma_left_init), c(q_left, q_left)))) {
      stop(sprintf("`Sigma_left_init` (or its alias) must be %d by %d.",
                   q_left, q_left), call. = FALSE)
    }

    # ----- Resolve Sigma_right_init from canonical + aliases -------------
    Sigma_right_aliases <- list(
      Sigma_right_init   = Sigma_right_init,
      Sigma_spatial_init = Sigma_spatial_init,
      Omega_G_init       = Omega_G_init
    )
    supplied_right <- names(Sigma_right_aliases)[!vapply(
      Sigma_right_aliases, is.null, logical(1L))]
    if (length(supplied_right) > 1L) {
      stop(sprintf(
        "Pass only ONE of: %s.", paste(supplied_right, collapse = ", ")),
        call. = FALSE)
    }
    if (length(supplied_right) == 1L) {
      Sigma_right_init <- Sigma_right_aliases[[supplied_right]]
    }
    if (is.null(Sigma_right_init)) {
      Sigma_right_init <- diag(n_groups)
    }
    if (!is.matrix(Sigma_right_init) ||
        !isTRUE(all.equal(dim(Sigma_right_init), c(n_groups, n_groups)))) {
      stop(sprintf("`Sigma_right_init` (or its alias) must be %d by %d.",
                   n_groups, n_groups), call. = FALSE)
    }

    re_cov_state <- .new_re_cov_state(
      type        = re_cov,
      Sigma_left  = Sigma_left_init,
      Sigma_right = Sigma_right_init,
      n_groups    = n_groups,
      q_left      = q_left
    )

    if (!isTRUE(control$update_variance)) {
      message(
        "Note: re_cov = '", re_cov, "' is typically run with ",
        "control$update_variance = TRUE so that Sigma_left is estimated ",
        "from the data. Pass vcmm_control(..., update_variance = TRUE) ",
        "to enable iterative estimation."
      )
    }
  }

  # ----- Resolve method = "auto" -----------------------------------------
  if (method == "auto") {
    if (n * q_dim > 1e5 || q_dim > 50L) {
      method_used <- "csl"
    } else {
      method_used <- "ss"
    }
    if (isTRUE(control$verbose)) {
      cat(sprintf("  [vcmm] method='auto' resolved to '%s' (N=%d, q=%d)\n",
                  method_used, n, q_dim))
    }
  } else {
    method_used <- method
  }

  # ----- Build design and penalty -----------------------------------------
  design <- build_vcmm_design(X           = X,
                              t           = t,
                              n_basis     = n_basis,
                              degree      = degree,
                              lambda      = lambda,
                              normalize_t = normalize_t)

  if (nrow(design$X_design) != n) {
    stop(sprintf("Internal: design rows (%d) != length(y) (%d).",
                 nrow(design$X_design), n),
         call. = FALSE)
  }

  # ----- Aggregated sufficient statistics ---------------------------------
  stats_obj <- compute_sufficient_stats(y, design$X_design, Z)

  # ----- Dispatch to the chosen estimator ---------------------------------
  fit <- switch(
    method_used,
    ss  = fit_ss (stats_obj, design$penalty, control,
                  re_cov_state = re_cov_state),
    csl = fit_csl(stats_obj, design$penalty, control,
                  re_cov_state = re_cov_state, ...)
  )

  # ----- Attach design metadata for downstream predict/plot ---------------
  fit$design <- design[c("internal_knots", "boundary_knots",
                         "degree", "n_basis", "K", "lambda",
                         "normalize_t", "t_min", "t_max")]
  fit$re_cov <- re_cov

  # ----- Identifiability post-processing ----------------------------------
  # When rowSums(Z) is constant (OD-style designs, ANOVA designs, etc.),
  # the model y = beta_0 + Z alpha + ... is unidentified along
  #   (beta_0, alpha) -> (beta_0 - s c, alpha + c 1).
  # Re-distribute the slack so sum(alpha_hat) = 0 and beta_0 absorbs the
  # shift. Zero-cost adjustment along the unidentified direction.
  row_sums <- rowSums(Z)
  if (max(row_sums) - min(row_sums) < 1e-10) {
    s_const <- row_sums[1L]
    if (abs(s_const) > 1e-10) {
      alpha_shift <- mean(fit$alpha)
      fit$alpha   <- fit$alpha - alpha_shift
      fit$beta[1L] <- fit$beta[1L] + s_const * alpha_shift
    }
  }

  fit$call <- match.call()
  fit
}

# Shared data generators used across the testthat suite.
#
# Each generator returns a list with consistent field names so the
# individual tests stay focused on what they're checking rather than
# data prep. All use deterministic seeds so failures are reproducible.

#' Small dataset with diagonal RE structure
gen_diag_data <- function(N = 300L, q = 3L, seed = 1L) {
  set.seed(seed)
  t <- runif(N)
  x <- runif(N)
  Z <- matrix(rnorm(N * q), N, q)
  alpha_true <- rnorm(q, sd = 0.4)
  y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% alpha_true) +
       rnorm(N, sd = 0.5)
  list(y = y, x = x, t = t, Z = Z, N = N, q = q, alpha_true = alpha_true)
}

#' Small OD-kron dataset
gen_kron_data <- function(N = 600L, G = 8L, seed = 2L) {
  set.seed(seed)
  q <- 2L * G
  origin_id <- sample.int(G, N, replace = TRUE)
  dest_id   <- sample.int(G, N, replace = TRUE)
  Z <- matrix(0, N, q)
  Z[cbind(seq_len(N), origin_id)]   <- 1
  Z[cbind(seq_len(N), G + dest_id)] <- 1
  Sigma_2x2     <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
  Sigma_spatial <- outer(seq_len(G), seq_len(G),
                         function(i, j) exp(-abs(i - j) / 4))
  alpha_true <- as.vector(crossprod(chol(kronecker(Sigma_2x2, Sigma_spatial)),
                                    rnorm(q)))
  t <- runif(N); x <- runif(N)
  y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% alpha_true) +
       rnorm(N, sd = 0.5)
  list(y = y, x = x, t = t, Z = Z,
       Sigma_2x2 = Sigma_2x2, Sigma_spatial = Sigma_spatial,
       G = G, q = q, N = N, alpha_true = alpha_true)
}

#' Small separable dataset
gen_separable_data <- function(N = 800L, G = 10L, q_left = 3L, seed = 3L) {
  set.seed(seed)
  q <- G * q_left
  group_id <- sample.int(G, N, replace = TRUE)
  Z <- matrix(0, N, q)
  for (k in seq_len(q_left)) {
    Z[cbind(seq_len(N), (k - 1L) * G + group_id)] <- rnorm(N)
  }
  Omega_G <- outer(seq_len(G), seq_len(G),
                   function(i, j) exp(-abs(i - j) / 5))
  t <- runif(N); x <- runif(N)
  y <- 2 + sin(2 * pi * t) * x + rnorm(N, sd = 0.5)
  list(y = y, x = x, t = t, Z = Z, Omega_G = Omega_G,
       G = G, q_left = q_left, q = q, N = N)
}

#' Symmetric positive-definite matrix at a given dimension
make_spd <- function(p, seed = 17L) {
  set.seed(seed + p)
  M <- matrix(rnorm(p * p), p, p)
  A <- crossprod(M) / p + diag(p)
  (A + t(A)) / 2
}

#' Symmetric matrix with negative eigenvalues (forces Cholesky failure)
make_indef <- function(p = 10L, seed = 4L) {
  set.seed(seed)
  M <- matrix(rnorm(p * p), p, p)
  (M + t(M)) / 2
}

#' Rank-deficient SPD matrix with controlled effective rank
make_rank_deficient_spd <- function(q, r, seed = 18L, cond_kept = 1e6) {
  set.seed(seed + q + r)
  U <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
  d <- numeric(q)
  d[seq_len(r)] <- exp(seq(0, -log(cond_kept), length.out = r))
  d[seq_len(q - r) + r] <- d[r] * 1e-14
  A <- U %*% (d * t(U))
  (A + t(A)) / 2
}

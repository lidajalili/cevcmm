#===============================================================================
# Day 10: Distributed sufficient-statistics aggregation
#
# Verifies that splitting the data into K mock nodes, computing summaries
# locally with node_summary(), aggregating server-side with `+`, and
# fitting with fit_from_summaries() reproduces the single-node vcmm()
# fit to floating-point tolerance.
#
# This is the strongest test we can write: bit-equivalence, not just
# statistical agreement. If the elementwise summation in `+.vcmm_ss` is
# correct, the aggregated stats are mathematically identical to the
# pooled stats (sums of crossprod blocks). The subsequent fit_ss /
# fit_csl call therefore sees the same input to ~1e-12 in every entry,
# and the resulting fit agrees on every estimated quantity.
#
# Uses the Day 8 OD setting so the comparison is anchored to a known
# realistic example: G = 20, q = 40 (2 per district), N = 3000.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

# ----- Truth (Day 8 OD setting, single replicate) ---------------------------
set.seed(10L)
G              <- 20L
q_total        <- 2L * G
N              <- 3000L
sigma_eps_true <- 0.5
beta0_true     <- 2.0
beta1_true_fn  <- function(tt) sin(2 * pi * tt)

Sigma_2x2_true     <- matrix(c(0.50, 0.20,
                               0.20, 0.50), 2L, 2L)
Sigma_spatial_true <- outer(seq_len(G), seq_len(G),
                            function(i, j) exp(-abs(i - j) / 5))
Sigma_alpha_true   <- kronecker(Sigma_2x2_true, Sigma_spatial_true)
L_alpha            <- chol(Sigma_alpha_true)

alpha_true <- as.vector(crossprod(L_alpha, rnorm(q_total)))
origin_id  <- sample.int(G, N, replace = TRUE)
dest_id    <- sample.int(G, N, replace = TRUE)

Z_mat <- matrix(0, N, q_total)
Z_mat[cbind(seq_len(N), origin_id)]      <- 1
Z_mat[cbind(seq_len(N), G + dest_id)]    <- 1

t_vec <- runif(N)
x_vec <- runif(N)
y_vec <- beta0_true +
         beta1_true_fn(t_vec) * x_vec +
         as.vector(Z_mat %*% alpha_true) +
         rnorm(N, sd = sigma_eps_true)

ctrl <- vcmm_control(sigma_eps       = sigma_eps_true,
                     sigma_alpha     = sqrt(0.5),
                     update_variance = TRUE,
                     max_iter        = 50L)

# ============================================================================
# (1) Reference fit: single-node vcmm()
# ============================================================================
cat("[Day 10] Reference fit: single-node vcmm()...\n")
t_ref <- system.time({
  fit_ref <- vcmm(y_vec,
                  X                  = x_vec,
                  Z                  = Z_mat,
                  t                  = t_vec,
                  method             = "csl",
                  re_cov             = "kronecker",
                  n_groups           = G,
                  Sigma_spatial_init = Sigma_spatial_true,
                  control            = ctrl)
})[["elapsed"]]
cat(sprintf("  done in %.2fs\n", t_ref))

# ============================================================================
# (2) Build the shared basis exactly as vcmm() did internally,
#     so each mock node uses an IDENTICAL X_design.
# ============================================================================
design_full <- build_vcmm_design(X = x_vec, t = t_vec)
X_design    <- design_full$X_design
penalty_mat <- design_full$penalty

# ============================================================================
# (3) Mock distributed fit: split into 3 nodes, summarize locally
# ============================================================================
n_nodes <- 3L
node_id <- sample.int(n_nodes, N, replace = TRUE)

cat(sprintf("\n[Day 10] Splitting N=%d obs across %d mock nodes:\n",
            N, n_nodes))
for (s in seq_len(n_nodes)) {
  cat(sprintf("  Node %d: %d obs\n", s, sum(node_id == s)))
}

summaries <- lapply(seq_len(n_nodes), function(s) {
  idx <- which(node_id == s)
  node_summary(y = y_vec[idx],
               X = X_design[idx, , drop = FALSE],
               Z = Z_mat[idx,    , drop = FALSE])
})

cat("\n[Day 10] One node summary (this is the only thing transmitted):\n")
print(summaries[[1L]])

# Show the size, to demonstrate communication efficiency
node_bytes <- object.size(summaries[[1L]])
raw_bytes  <- object.size(list(y_vec[node_id == 1L],
                               X_design[node_id == 1L, ],
                               Z_mat[node_id == 1L, ]))
cat(sprintf(
  "\n[Day 10] Node 1 transmission cost: %.1f KB summary vs %.1f KB raw data (%.0fx reduction)\n",
  as.numeric(node_bytes) / 1024,
  as.numeric(raw_bytes)  / 1024,
  as.numeric(raw_bytes)  / as.numeric(node_bytes)))

# Server-side aggregation: `Reduce("+", ...)` uses the +.vcmm_ss method
agg <- Reduce(`+`, summaries)
cat("\n[Day 10] Aggregated summary (server side):\n")
print(agg)

# ============================================================================
# (4) Server fits from aggregated summary
# ============================================================================
cat("\n[Day 10] fit_from_summaries (list path)...\n")
t_dist <- system.time({
  fit_dist <- fit_from_summaries(
    summaries,
    penalty            = penalty_mat,
    control            = ctrl,
    method             = "csl",
    re_cov             = "kronecker",
    n_groups           = G,
    Sigma_spatial_init = Sigma_spatial_true,
    rowsum_constant    = 2     # vcmm() applies this re-centering for
                               # constant-rowSum Z; tell the server too.
  )
})[["elapsed"]]
cat(sprintf("  done in %.2fs\n", t_dist))

# Also verify the "already-aggregated single summary" path
cat("\n[Day 10] fit_from_summaries (single-aggregate path)...\n")
fit_dist_agg <- fit_from_summaries(
  agg,
  penalty            = penalty_mat,
  control            = ctrl,
  method             = "csl",
  re_cov             = "kronecker",
  n_groups           = G,
  Sigma_spatial_init = Sigma_spatial_true,
  rowsum_constant    = 2
)

# ============================================================================
# (5) Bit-equivalence checks
# ============================================================================
beta_diff       <- max(abs(fit_ref$beta      - fit_dist$beta))
alpha_diff      <- max(abs(fit_ref$alpha     - fit_dist$alpha))
sigma_eps_diff  <- abs   (fit_ref$sigma_eps  - fit_dist$sigma_eps)
Sigma_left_diff <- max(abs(fit_ref$re_cov_state$Sigma_left -
                           fit_dist$re_cov_state$Sigma_left))

list_vs_agg <- max(c(
  max(abs(fit_dist$beta      - fit_dist_agg$beta)),
  max(abs(fit_dist$alpha     - fit_dist_agg$alpha)),
  abs   (fit_dist$sigma_eps  - fit_dist_agg$sigma_eps)
))

cat("\n========== Reference vs distributed (max abs diff) ==========\n")
cat(sprintf("  beta              : %.3e\n", beta_diff))
cat(sprintf("  alpha             : %.3e\n", alpha_diff))
cat(sprintf("  sigma_eps         : %.3e\n", sigma_eps_diff))
cat(sprintf("  Sigma_left        : %.3e\n", Sigma_left_diff))
cat(sprintf("  list path vs agg path : %.3e\n", list_vs_agg))

# ----- Pass / fail ----------------------------------------------------------
tol  <- 1e-7
pass <- c(
  beta_match        = beta_diff       < tol,
  alpha_match       = alpha_diff      < tol,
  sigma_eps_match   = sigma_eps_diff  < tol,
  Sigma_left_match  = Sigma_left_diff < tol,
  list_eq_agg_path  = list_vs_agg     < tol
)

cat("\n========== Day 10 pass / fail ==========\n")
for (nm in names(pass)) {
  cat(sprintf("  %-22s %s\n", nm, if (pass[[nm]]) "PASS" else "FAIL"))
}
cat(sprintf("\nOverall: %s\n",
            if (all(pass)) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"))

invisible(list(fit_ref      = fit_ref,
               fit_dist     = fit_dist,
               fit_dist_agg = fit_dist_agg,
               summaries    = summaries,
               agg          = agg))

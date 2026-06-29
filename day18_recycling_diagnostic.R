#===============================================================================
# Diagnose: find which file emits the "d * t(U)" recycling warning.
#
# Wraps source(day18_rspectra_validation.R) in a calling handler that
# intercepts ANY warning containing "longer object length", prints the
# full R call stack at the moment it fires, and continues. The first
# such warning will reveal the exact function and line.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

caught_any <- FALSE

withCallingHandlers(
  source("inst/validation/day18_rspectra_validation.R"),
  warning = function(w) {
    if (grepl("longer object length", conditionMessage(w))) {
      caught_any <<- TRUE
      cat("\n\n========================================\n")
      cat("RECYCLING WARNING CAUGHT:\n")
      cat("  '", conditionMessage(w), "'\n", sep = "")
      cat("\nCALL STACK (most recent last):\n")
      calls <- sys.calls()
      for (i in seq_along(calls)) {
        # Skip the source() and handler frames at top of stack
        txt <- deparse(calls[[i]], width.cutoff = 200, nlines = 2)
        cat(sprintf("  [%2d] %s\n", i, paste(txt, collapse = " ")))
      }
      cat("========================================\n\n")
      invokeRestart("muffleWarning")
    }
  }
)

if (!caught_any) {
  cat("\n\nNo recycling warning fired. The patch is working.\n")
} else {
  cat("\n\nWarning fired -- see traceback above for source.\n")
}

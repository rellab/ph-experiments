# script for experiment1

source("script.R")

# n set from args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Missing argument: n. Usage: Rscript experiment1_100.R <n>")
}
n <- as.integer(args[1])
if (is.na(n)) {
    stop("Invalid integer for n: ", args[1])
}

mlist <- c(5L, 10L, 20L, 30L, 50L, 70L, 100L)

method <- "UNI"
steps <- 10L
abstol <- 1.0e-3
reltol <- 1.0e-5
maxiter <- 5000L

# for (uni_epsilon in c(1.0e-8, 1.0e-6)) {
for (uni_epsilon in c(1.0e-8)) {
  cat(sprintf("Running matrixdist for general PH and right_censored_sample_%d; tolerance of UNI = %e\n", n, uni_epsilon))
  data_file <- paste("data/right_censored_sample_", n, ".txt", sep = "")
  for (m in mlist) {
    run_fit(data_file, paste("params/general_phase", m, "_", n, ".txt", sep = ""), method, uni_epsilon, steps, abstol, reltol, maxiter)
  }
  cat("Done.\n")

  cat(sprintf("Running matrixdist for bidiagonal PH and right_censored_sample_%d; tolerance of UNI = %e\n", n, uni_epsilon))
  data_file <- paste("data/right_censored_sample_", n, ".txt", sep = "")
  for (m in mlist) {
    run_fit(data_file, paste("params/bidiagonal_phase", m, "_", n, ".txt", sep = ""), method, uni_epsilon, steps, abstol, reltol, maxiter)
  }
  cat("Done.\n")
}

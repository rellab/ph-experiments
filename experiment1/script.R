# install.packages("matrixdist")
library(matrixdist)

# Read right-censored data
read_rc_data <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (any(grepl("^\\s*-?1\\s*$", lines))) {
    lines <- lines[seq_len(which(grepl("^\\s*-?1\\s*$", lines))[1] - 1)]
  }
  lines <- lines[!grepl("^\\s*$", lines)]
  df <- read.table(text = paste0(lines, collapse = "\n"),
                   col.names = c("delta", "time"))
  list(
    y    = df$time[df$delta == 1],  # observed times
    rcen = df$time[df$delta == 0]   # right-censored times
  )
}

# Read initial PH (Coxian assumption) parameters (m rows × (m+1) columns: alpha|S)
# Example:
#  0.9074368532077339  -2.0088593068970515   2.0088593068970515
#  0.09256314679226413  0.0                  -2.0224913991942643
# → m=2, alpha=c(0.9074, 0.09256),
#   S = matrix(c(-2.0088593,  2.0088593,
#                 0.0,       -2.0224914), nrow=2, byrow=TRUE)

read_alpha_S <- function(path) {
  M <- as.matrix(read.table(path, check.names = FALSE))
  m <- nrow(M)
  stopifnot(ncol(M) == m + 1)
  alpha <- as.numeric(M[, 1])
  S     <- t(matrix(as.numeric(M[, -1, drop = FALSE]), nrow = m, byrow = TRUE))
  list(alpha = alpha, S = S)
}

# Fit PH model and measure time
# Return elapsed time and log-likelihood
# method: "RK", "UNI", "PADE"
# uni_eps: UNI method epsilon (ignored for non-UNI)
# em_steps: number of EM steps
time_fit <- function(init_ph, rc, method, uni_eps = NA_real_, em_steps = 100L) {
  fit_result <- NULL
  t <- system.time({
    fit_result <- fit(
      x        = init_ph,
      y        = rc$y,
      rcen     = rc$rcen,
      stepsEM  = em_steps,
      methods  = c(method, method), # "RK" / "UNI" / "PADE"
      every    = 100L,               # Display progress every step
      uni_epsilon = uni_eps         # Tolerance for UNI (ignored for non-UNI)
    )
  })

  return(list(
    elapsed = unname(t["elapsed"]),
    loglik = logLik(fit_result),
    fit_result = fit_result
  ))
}

# Run fitting until convergence
# method: "RK", "UNI", "PADE"
# uni_eps: UNI method epsilon (ignored for non-UNI)
# em_steps: number of EM steps between checkpoints of convergence
# (e.g., 10L means check convergence every 10 EM steps)
# abstol: absolute tolerance for convergence
# reltol: relative tolerance for convergence
# maxiter: maximum number of iterations
fit_until_convergence <- function(init_ph, rc, method, uni_eps = NA_real_,
                                  em_steps = 10L, abstol = 1e-2, reltol = 1e-6, maxiter = 5000L) {
  prev_llf <- -Inf
  current_llf <- -Inf
  total_steps <- 0L
  rerror <- NA_real_
  capture.output({
    result <- fit(
        x        = init_ph,
        y        = rc$y,
        rcen     = rc$rcen,
        stepsEM  = em_steps,
        methods  = c(method, method), # "RK" / "UNI" / "PADE"
        every    = 100L,               # Display progress every step
        uni_epsilon = uni_eps        # Tolerance for UNI (ignored for non-UNI)
    )
  }, file = NULL)
  while (total_steps < maxiter) {
    total_steps <- total_steps + em_steps
    current_llf <- logLik(result)
    aerror <- abs(current_llf - prev_llf)
    rerror <- aerror / abs(prev_llf)
    # cat(sprintf("After %d steps, log-likelihood: %.10f\n", total_steps, current_llf))
    if (abs(current_llf - prev_llf) < reltol * abs(prev_llf) || aerror < abstol) {
      cat("\nConverged.\n")
      return(list(status = "Converged", ph = result, aerror = aerror, rerror = rerror, total_steps = total_steps))
    }
    prev_llf <- current_llf
    capture.output({
      result <- fit(
        x        = result,
        y        = rc$y,
        rcen     = rc$rcen,
        stepsEM  = em_steps,
        methods  = c(method, method), # "RK" / "UNI" / "PADE"
        every    = 100L,               # Display progress every step
        uni_epsilon = uni_eps        # Tolerance for UNI (ignored for non-UNI)
      )
    }, file = NULL)
    # cat("#")
  }
  cat("\nReached maximum iterations without convergence.\n")
  return(list(status = "MaxIter", ph = result, aerror = aerror, rerror = rerror, total_steps = total_steps))
}

time_fit_until_convergence <- function(init_ph, rc, method, uni_eps = NA_real_,
                                  em_steps = 10L, abstol = 1e-2, reltol = 1e-6, maxiter = 5000L) {
  fit_result <- NULL
  t <- system.time({
    fit_result <- fit_until_convergence(
      init_ph, rc, method, uni_eps,
      em_steps, abstol, reltol, maxiter
    )
  })

  return(list(
    elapsed = unname(t["elapsed"]),
    loglik = logLik(fit_result$ph),
    result = fit_result
  ))
}

# Write final alpha and S to file
write_alpha_S <- function(final_model, output_file) {
  # Get parameters
  final_alpha <- coef(final_model)$alpha
  final_S <- coef(final_model)$S

  cat("Saving final parameters to:", output_file, "\n")

  # Write in alpha | S format (same as initial file)
  m <- length(final_alpha)
  output_matrix <- cbind(final_alpha, final_S)

  write.table(output_matrix, file = output_file, 
              row.names = FALSE, col.names = FALSE, 
              sep = "  ", quote = FALSE)
}

# ==== コマンドライン引数の取得 ====
run_fit <- function(data_file, phase_file, method="UNI", uni_epsilon=1e-8, steps=10L, abstol=1e-2, reltol=1e-6, maxiter=5000L) {
  # args <- commandArgs(trailingOnly = TRUE)
  # data_file <- args[1]   # 例: "unweighted10000"
  # phase_file <- args[2]  # 例: "phase30"
  # method <- if (length(args) >= 3) args[3] else "UNI"  # デフォルト: UNI
  # uni_epsilon <- if (length(args) >= 4) as.numeric(args[4]) else 1e-8  # デフォルト: 1e-8
  # steps <- if (length(args) >= 5) as.integer(args[5]) else 10L  # デフォルト: 10
  # abstol <- if (length(args) >= 6) as.numeric(args[6]) else 1e-2  # デフォルト: 1e-2
  # reltol <- if (length(args) >= 6) as.numeric(args[6]) else 1e-6  # デフォルト: 1e-6
  # maxiter <- if (length(args) >= 7) as.integer(args[7]) else 5000L  # デフォルト: 5000

  # GC
  gc()

  cat("Starting fit with parameters:\n")
  cat(sprintf("  data_file: %s\n", data_file))
  cat(sprintf("  phase_file: %s\n", phase_file))

  rc <- read_rc_data(data_file)
  y_obs    <- rc$y
  rcen_obs <- rc$rcen
  cat(sprintf("n=%d (events=%d, r-censored=%d)\n",
              length(y_obs)+length(rcen_obs), length(y_obs), length(rcen_obs)))

  init <- read_alpha_S(phase_file)
  m_dim <- length(init$alpha)
  cat(sprintf("m=%d\n", m_dim))

  init_ph <- ph(alpha = init$alpha, S = init$S)

  if (!(method %in% c("RK", "UNI", "PADE"))) {
    cat("Error: method must be one of: RK, UNI, PADE\n")
    quit(status = 1)
  }

  cat(sprintf("Using method: %s\n", method))
  if (method == "UNI") {
    cat(sprintf("UNI epsilon: %e\n", uni_epsilon))
  }

  fit_results <- time_fit_until_convergence(init_ph, rc, method, uni_epsilon,
    em_steps=steps, abstol=abstol, reltol=reltol, maxiter=maxiter)
  # print(fit_results)

  # final_model <- fit_results$fit_result
  # output_file <- paste0("results/result_R_", gsub(".*/(.*)", "\\1", data_file), "_", gsub(".*/(.*)", "\\1", phase_file), "_", method, ".txt")
  # write_alpha_S(final_model, output_file)

  # 結果のllfをllf_paramX.txtに保存
  # llf_file <- paste0("llf_param_", gsub(".*/(.*)", "\\1", data_file), "_", gsub(".*/(.*)", "\\1", phase_file), "_", method, ".txt")
  # cat("Saving final log-likelihood to:", llf_file, "\n")
  # write.table(data.frame(loglik = fit_results$loglik), file = llf_file,
  #             row.names = FALSE, col.names = TRUE,
  #             sep = "  ", quote = FALSE)

  results <- data.frame(
    method = method,
    steps  = steps,
    reltol = reltol,
    maxiter = maxiter,
    uni_epsilon = if (method == "UNI") uni_epsilon else NA_real_,
    elapsed_sec = fit_results$elapsed,
    final_loglik = fit_results$loglik,
    aerror = fit_results$result$aerror,
    rerror = fit_results$result$rerror,
    status = fit_results$result$status,
    total_steps = fit_results$result$total_steps,
    row.names = NULL
  )

  print(results)
}


source("script.R")

makeinit_general_ph <- function(m, rc, output_file, seed = 1234) {
  set.seed(seed)
  init_ph <- ph(structure="general", dimension=m)
  result <- time_fit(init_ph, rc, method="UNI", uni_eps=1e-8, em_steps=5L)
  write_alpha_S(result$fit_result, output_file)
}

## gen bidiagonal PH: alpha is a random vector and S is a bidiagonal matrix
generate_bidiagonal_ph <- function(m, seed = 1234) {
    set.seed(seed)
    # Generate random initial probability vector alpha
    alpha <- runif(m)
    alpha <- alpha / sum(alpha)
    # Generate bidiagonal matrix S
    S <- matrix(0, nrow = m, ncol = m)
    for (i in 1:m) {
        rate <- runif(1, min = 0.1, max = 1)
        S[i, i] <- -rate
        if (i < m) {
            S[i, i + 1] <- rate
        }
    }
    list(alpha = alpha, S = S)
}

makeinit_bidiagonal_ph <- function(m, rc, output_file, seed = 1234) {
  set.seed(seed)
  bidiagonal_params <- generate_bidiagonal_ph(m, seed)
  init_ph <- ph(alpha = bidiagonal_params$alpha, S = bidiagonal_params$S)
  result <- time_fit(init_ph, rc, method="UNI", uni_eps=1e-8, em_steps=5L)
  write_alpha_S(result$fit_result, output_file)
}

# initial general PH parameters for unweighted100
rc <- read_rc_data("data/unweighted100")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/general_phase", m, "_100.txt")
    result <- makeinit_general_ph(m, rc, output_file)
}

# initial bidiagonal parameters for unweighted100
rc <- read_rc_data("data/unweighted100")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/bidiagonal_phase", m, "_100.txt")
    result <- makeinit_bidiagonal_ph(m, rc, output_file)
}

# initial general PH parameters for unweighted1000
rc <- read_rc_data("data/unweighted1000")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/general_phase", m, "_1000.txt")
    result <- makeinit_general_ph(m, rc, output_file)
}

# initial bidiagonal parameters for unweighted1000
rc <- read_rc_data("data/unweighted1000")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/bidiagonal_phase", m, "_1000.txt")
    result <- makeinit_bidiagonal_ph(m, rc, output_file)
}

# initial general PH parameters for unweighted10000
rc <- read_rc_data("data/unweighted10000")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/general_phase", m, "_10000.txt")
    result <- makeinit_general_ph(m, rc, output_file)
}

# initial bidiagonal parameters for unweighted10000
rc <- read_rc_data("data/unweighted10000")
for (m in c(2, 5, 10, 20, 30, 50, 70, 100)) {
    output_file <- paste0("params/bidiagonal_phase", m, "_10000.txt")
    result <- makeinit_bidiagonal_ph(m, rc, output_file)
}



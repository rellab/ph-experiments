# script for experiment1

source("script.R")

method <- "UNI"
uni_epsilon <- 1.0e-8
steps <- 10L
abstol <- 1.0e-3
reltol <- 1.0e-5
maxiter <- 5000L

cat("Running matrixdist for general PH and unweighted1000; tolerance of UNI = 1e-6\n")
data_file <- "data/unweighted1000"
uni_epsilon <- 1.0e-6
run_fit(data_file, "params/general_phase5_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase10_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase20_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase30_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase50_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase70_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase100_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
cat("Done.\n")

cat("Running matrixdist for bidiagonal PH and unweighted1000; tolerance of UNI = 1e-6\n")
data_file <- "data/unweighted1000"
uni_epsilon <- 1.0e-6
run_fit(data_file, "params/bidiagonal_phase5_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase10_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase20_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase30_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase50_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase70_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase100_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
cat("Done.\n")

### tolerance of UNI = 1e-8

cat("Running matrixdist for general PH and unweighted1000; tolerance of UNI = 1e-8\n")
data_file <- "data/unweighted1000"
uni_epsilon <- 1.0e-8
run_fit(data_file, "params/general_phase5_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase10_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase20_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase30_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase50_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase70_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase100_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
cat("Done.\n") 

cat("Running matrixdist for bidiagonal PH and unweighted1000; tolerance of UNI = 1e-8\n")
data_file <- "data/unweighted1000"
uni_epsilon <- 1.0e-8
run_fit(data_file, "params/bidiagonal_phase5_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase10_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase20_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase30_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase50_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase70_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase100_1000.txt", method, uni_epsilon, steps, abstol, reltol, maxiter)
cat("Done.\n")

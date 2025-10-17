# script for experiment1

include("phlib.jl")

steps = 10
abstol = 1.0e-3
reltol = 1.0e-5
maxiter = 5000

println("Running proposed method for general PH and unweighted1000")
data_file = "data/unweighted1000"
run_fit(data_file, "params/general_phase2_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase5_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase10_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase20_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase30_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase50_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase70_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase100_1000.txt", steps, abstol, reltol, maxiter)
println("Done.")

println("Running proposed method for bidiagonal PH and unweighted1000")
data_file = "data/unweighted1000"
run_fit(data_file, "params/bidiagonal_phase2_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase5_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase10_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase20_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase30_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase50_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase70_1000.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase100_1000.txt", steps, abstol, reltol, maxiter)
println("Done.")

println("Running proposed method for CF1 and unweighted1000")
data_file = "data/unweighted1000"
run_fit_cf1(data_file, "params/bidiagonal_phase2_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase5_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase10_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase20_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase30_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase50_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase70_1000.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase100_1000.txt", steps, abstol, reltol, maxiter)
println("Done.")

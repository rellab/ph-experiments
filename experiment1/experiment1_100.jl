# script for experiment1

include("phlib.jl")

steps = 10
abstol = 1.0e-3
reltol = 1.0e-5
maxiter = 5000

println("Running proposed method for general PH and unweighted100")
data_file = "data/unweighted100"
run_fit(data_file, "params/general_phase2_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase5_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase10_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase20_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase30_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase50_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase70_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/general_phase100_100.txt", steps, abstol, reltol, maxiter)
println("Done.")

println("Running proposed method for bidiagonal PH and unweighted100")
data_file = "data/unweighted100"
run_fit(data_file, "params/bidiagonal_phase2_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase5_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase10_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase20_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase30_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase50_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase70_100.txt", steps, abstol, reltol, maxiter)
run_fit(data_file, "params/bidiagonal_phase100_100.txt", steps, abstol, reltol, maxiter)
println("Done.")

println("Running proposed method for CF1 and unweighted100")
data_file = "data/unweighted100"
run_fit_cf1(data_file, "params/bidiagonal_phase2_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase5_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase10_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase20_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase30_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase50_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase70_100.txt", steps, abstol, reltol, maxiter)
run_fit_cf1(data_file, "params/bidiagonal_phase100_100.txt", steps, abstol, reltol, maxiter)
println("Done.")

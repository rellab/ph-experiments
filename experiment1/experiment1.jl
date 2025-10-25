# script for experiment1

include("phlib.jl")

# get n from args
n = parse(Int, ARGS[1])

steps = 10
abstol = 1.0e-3
reltol = 1.0e-5
maxiter = 5000

mlist = [2, 5, 10, 20, 30, 50, 70, 100]

println("Running proposed method for general PH; n=$n")
data_file = "data/right_censored_sample_$n.txt"
for m in mlist
    run_fit(data_file, "params/general_phase$(m)_$n.txt", steps, abstol, reltol, maxiter)
end
println("Done.")

println("Running proposed method for bidiagonal PH n=$n")
data_file = "data/right_censored_sample_$n.txt"
for m in mlist
    run_fit(data_file, "params/bidiagonal_phase$(m)_$n.txt", steps, abstol, reltol, maxiter)
end
println("Done.")

println("Running proposed method for CF1 n=$n")
data_file = "data/right_censored_sample_$n.txt"
for m in mlist
    run_fit_cf1(data_file, "params/bidiagonal_phase$(m)_$n.txt", steps, abstol, reltol, maxiter)
end
println("Done.")

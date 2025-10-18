using JLD2
using Random
using LinearAlgebra
BLAS.set_num_threads(1)

include("phlib.jl")

maxiter = 5000
steps = 10
reltol = 1e-5
bsample = 200

mlist = [2, 5, 10, 20, 30, 50, 70, 100]
#mlist = [2, 5]
results = Dict{Tuple{String, Int}, Any}()

println("=============================================")
println("Compute EIC for Lognormal sample (n=1000)")
rng = MersenneTwister(1234)
dat = read_rc_data("data/lognormal_sample.txt")
for m in mlist
    GC.gc() # garbage collection
    println("Starting EIC computation for m=$m....")
    t = @elapsed begin
        ph0, llf0, _, _, _, d0 = phfit(CF1(m), to_data(dat); initialize=true, maxiter=maxiter, steps=steps, reltol=reltol)
        value, upper, lower = eic(rng, ph0, llf0, d0, dat, bsample=bsample, maxiter=maxiter, steps=steps, reltol=reltol)
    end
    results[("Lognormal", m)] = ph0
    @printf("\nm=%d, LLF=%.12f, AIC=%.12f, EIC=%.12f, 95%% CI=(%.12f, %.12f), time=%.2f\n", m, llf0, -2*(llf0 - (2*m-1)), value, lower, upper, t)
    println("--------------------------------------------------")
end

println("=============================================")
println("Compute EIC for Weibull sample (n=1000)")
rng = MersenneTwister(1234)
dat = read_rc_data("data/weibull_sample.txt")
for m in mlist
    GC.gc() # garbage collection
    println("Starting EIC computation for m=$m....")
    t = @elapsed begin
        ph0, llf0, _, _, _, d0 = phfit(CF1(m), to_data(dat); initialize=true, maxiter=maxiter, steps=steps, reltol=reltol)
        value, upper, lower = eic(rng, ph0, llf0, d0, dat, bsample=bsample, maxiter=maxiter, steps=steps, reltol=reltol)
    end
    results[("Weibull", m)] = ph0
    @printf("\nm=%d, LLF=%.12f, AIC=%.12f, EIC=%.12f, 95%% CI=(%.12f, %.12f), time=%.2f\n", m, llf0, -2*(llf0 - (2*m-1)), value, lower, upper, t)
    println("--------------------------------------------------")
end

println("=============================================")
println("Compute EIC for Mixture model sample (n=1000)")
rng = MersenneTwister(1234)
dat = read_rc_data("data/mixturemodel_sample.txt")
for m in mlist
    GC.gc() # garbage collection
    println("Starting EIC computation for m=$m....")
    t = @elapsed begin
        ph0, llf0, _, _, _, d0 = phfit(CF1(m), to_data(dat); initialize=true, maxiter=maxiter, steps=steps, reltol=reltol)
        value, upper, lower = eic(rng, ph0, llf0, d0, dat, bsample=bsample, maxiter=maxiter, steps=steps, reltol=reltol)
    end
    results[("Mixture", m)] = ph0
    @printf("\nm=%d, LLF=%.12f, AIC=%.12f, EIC=%.12f, 95%% CI=(%.12f, %.12f), time=%.2f\n", m, llf0, -2*(llf0 - (2*m-1)), value, lower, upper, t)
    println("--------------------------------------------------")
end

# save results to a file: we require the file to draw graphs later
@save "results/phparams200.jld2" results

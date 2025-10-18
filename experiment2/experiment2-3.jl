using JLD2
using Plots
using Distributions
using SpecialFunctions

include("phlib.jl")

@load "results/phparams200.jld2" results

println("=================================================")
println("Lognormal distribution")

mu, sigma = 1.0, 0.5
dist = LogNormal(mu, sigma)

phaic = results[("Lognormal", 5)]
pheic = results[("Lognormal", 10)]
phcross = results[("Lognormal", 20)]

xx = LinRange(0, 10.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, xlabel="t", ylabel="Probability density", lw=2, label="Original", ls=:solid)
plot!(xx, yyaic, lw=2, label="AIC (m=5)", ls=:dash)
plot!(xx, yyeic, lw=2, label="EIC (m=10)", ls=:dot)
p = plot!(xx, yycross, lw=2, label="Cross (m=20)", ls=:dash)
savefig(p, "resultslognormal-phfit.pdf")

println("=================================================")
println("Weibull distribution")

alpha = 5.0
beta = 8.0
dist = Weibull(alpha, beta)

phaic = results[("Weibull", 10)]
pheic = results[("Weibull", 30)]
phcross = results[("Weibull", 50)]

xx = LinRange(0, 15.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, xlabel="t", ylabel="Probability density", lw=2, label="Original", ls=:solid)
plot!(xx, yyaic, lw=2, label="AIC (m=10)", ls=:dash)
plot!(xx, yyeic, lw=2, label="EIC (m=30)", ls=:dot)
p = plot!(xx, yycross, lw=2, label="Cross (m=50)", ls=:dash)
savefig(p, "resultsweibull-phfit.pdf")

println("=================================================")
println("Mixture model")

mu = 1.0
sigma = 0.5
alpha = 5.0
beta = 8.0
p = 0.3
dist = MixtureModel([LogNormal(mu, sigma), Weibull(alpha, beta)], [p, 1-p])

phaic = results[("Mixture", 10)]
pheic = results[("Mixture", 20)]
phcross = results[("Mixture", 70)]

xx = LinRange(0, 15.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, xlabel="t", ylabel="Probability density", lw=2, label="Original", ls=:solid)
plot!(xx, yyaic, lw=2, label="AIC (m=10)", ls=:dash)
plot!(xx, yyeic, lw=2, label="EIC (m=20)", ls=:dot)
p = plot!(xx, yycross, lw=2, label="Cross (m=70)", ls=:dash)
savefig(p, "results/mixture-phfit.pdf")


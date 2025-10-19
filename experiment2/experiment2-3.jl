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
# plot(xx, yy, xlabel="t", ylabel="Probability density", lw=2, label="Original", ls=:solid)
# plot!(xx, yyaic, lw=2, label="AIC (m=5)", ls=:dash)
# plot!(xx, yyeic, lw=2, label="EIC (m=10)", ls=:dot)
# p = plot!(xx, yycross, lw=2, label="Cross (m=20)", ls=:dash)
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=5)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=10)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=20)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)
savefig(p, "results/lognormal-phfit.pdf")

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
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=10)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=30)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=50)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)
savefig(p, "results/weibull-phfit.pdf")

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
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=10)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=20)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=70)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)

savefig(p, "results/mixture-phfit.pdf")


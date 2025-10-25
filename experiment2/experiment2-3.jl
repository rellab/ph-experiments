using JLD2
using Plots
using Distributions
using SpecialFunctions

include("phlib.jl")

# get n, (maic, meic, mcross), (maic, meic, mcross), (maic, meic, mcross) from args
n = parse(Int, ARGS[1])
maic1, meic1, mcross1 = parse(Int, ARGS[2]), parse(Int, ARGS[3]), parse(Int, ARGS[4])
maic2, meic2, mcross2 = parse(Int, ARGS[5]), parse(Int, ARGS[6]), parse(Int, ARGS[7])
maic3, meic3, mcross3 = parse(Int, ARGS[8]), parse(Int, ARGS[9]), parse(Int, ARGS[10])

@load "results/phparams_$n.jld2" results

println("=================================================")
println("Lognormal distribution")

mu, sigma = 1.0, 0.5
dist = LogNormal(mu, sigma)

# please rewrite m to select the model having the smallest criterion
maic, meic, mcross = maic1, meic1, mcross1
phaic = results[("Lognormal", maic)]
pheic = results[("Lognormal", meic)]
phcross = results[("Lognormal", mcross)]

xx = LinRange(0, 10.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=$maic)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=$meic)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=$mcross)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)
savefig(p, "results/lognormal-phfit_$n.pdf")

println("=================================================")
println("Weibull distribution")

alpha = 5.0
beta = 8.0
dist = Weibull(alpha, beta)

# please rewrite m to select the model having the smallest criterion
maic, meic, mcross = maic2, meic2, mcross2
phaic = results[("Weibull", maic)]
pheic = results[("Weibull", meic)]
phcross = results[("Weibull", mcross)]

xx = LinRange(0, 15.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=$maic)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=$meic)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=$mcross)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)
savefig(p, "results/weibull-phfit_$n.pdf")

println("=================================================")
println("Mixture model")

mu = 1.0
sigma = 0.5
alpha = 5.0
beta = 8.0
p = 0.3
dist = MixtureModel([LogNormal(mu, sigma), Weibull(alpha, beta)], [p, 1-p])

# please rewrite m to select the model having the smallest criterion
maic, meic, mcross = maic3, meic3, mcross3
phaic = results[("Mixture", maic)]
pheic = results[("Mixture", meic)]
phcross = results[("Mixture", mcross)]

xx = LinRange(0, 15.0, 10000)
yy = pdf.(dist, xx)
yyaic = phpdf(phaic, xx)
yyeic = phpdf(pheic, xx)
yycross = phpdf(phcross, xx)
plot(xx, yy, lw=2, color=:black, ls=:solid,
     label="Original")

plot!(xx, yyaic, lw=2, color=:red, ls=:dash,
      label="AIC (m=$maic)")

plot!(xx, yyeic, lw=2, color=:blue, ls=:dot,
      label="EIC (m=$meic)")

plot!(xx, yycross, lw=2, color=:green, ls=:dashdot,
      label="Cross (m=$mcross)")

p = plot!(legend=:topright, grid=false,
      xlabel="t", ylabel="Probability density",
      background_color=:white, framestyle=:box)

savefig(p, "results/mixture-phfit_$n.pdf")


using JLD2
using Plots
using Distributions
using SpecialFunctions

include("phlib.jl")

@load "results/phparams200.jld2" results

mlist = [2, 5, 10, 20, 30, 50, 70, 100]

println("=================================================")
println("Lognormal distribution")

mu, sigma = 1.0, 0.5
dist = LogNormal(mu, sigma)
m1 = exp(1*mu + 0.5*1^2*sigma^2)
m2 = exp(2*mu + 0.5*2^2*sigma^2)
m3 = exp(3*mu + 0.5*3^2*sigma^2)
@printf("First three moments (original): %.8f, %.8f, %.8f\n", m1, m2, m3)

for m in mlist
    phmodel = results[("Lognormal", m)]
    m1 = phmean(phmodel, 1)
    m2 = phmean(phmodel, 2)
    m3 = phmean(phmodel, 3)
    cr = cross(dist, phmodel)
    @printf("First three moments (m=%d): %.8f, %.8f, %.8f\n", m, m1, m2, m3)
    @printf("Cross-moment (m=%d): %.8f\n", m, cr)
end

println("=================================================")
println("Weibull distribution")

alpha = 5.0
beta = 8.0
dist = Weibull(alpha, beta)
m1 = beta * gamma(1 + 1/alpha)
m2 = beta^2 * gamma(1 + 2/alpha)
m3 = beta^3 * gamma(1 + 3/alpha)
@printf("First three moments (original): %.8f, %.8f, %.8f\n", m1, m2, m3)

for m in mlist
    phmodel = results[("Weibull", m)]
    m1 = phmean(phmodel, 1)
    m2 = phmean(phmodel, 2)
    m3 = phmean(phmodel, 3)
    cr = cross(dist, phmodel)
    @printf("First three moments (m=%d): %.8f, %.8f, %.8f\n", m, m1, m2, m3)
    @printf("Cross-moment (m=%d): %.8f\n", m, cr)
end

println("=================================================")
println("Mixture model")

mu = 1.0
sigma = 0.5
alpha = 5.0
beta = 8.0
p = 0.3
dist = MixtureModel([LogNormal(mu, sigma), Weibull(alpha, beta)], [p, 1-p])

m1 = p * exp(1*mu + 0.5*1^2*sigma^2) + (1-p) * beta * gamma(1 + 1/alpha)
m2 = p * exp(2*mu + 0.5*2^2*sigma^2) + (1-p) * beta^2 * gamma(1 + 2/alpha)
m3 = p * exp(3*mu + 0.5*3^2*sigma^2) + (1-p) * beta^3 * gamma(1 + 3/alpha)
@printf("First three moments (original): %.8f, %.8f, %.8f\n", m1, m2, m3)

for m in mlist
    phmodel = results[("Mixture", m)]
    m1 = phmean(phmodel, 1)
    m2 = phmean(phmodel, 2)
    m3 = phmean(phmodel, 3)
    cr = cross(dist, phmodel)
    @printf("First three moments (m=%d): %.8f, %.8f, %.8f\n", m, m1, m2, m3)
    @printf("Cross-moment (m=%d): %.8f\n", m, cr)
end

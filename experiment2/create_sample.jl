using PhaseTypeDistributions
using PhaseTypeDistributions.Phfit
# using CSV
# using DataFrames
# using Plots
using Distributions
using Random
using Deformula
using Test
using Printf
# using HDF5
# using Logging

function makeSample(rng, dist, trunc_prob, trunc_dist, censoring_dist)
    if trunc_prob > rand(rng, Float64)
        tau = rand(rng, trunc_dist)
    else
        tau = 0.0
    end

    nn = 0
    x = rand(rng, dist)
    while x < tau
        x = rand(rng, dist)
        nn += 1
        if nn > 10000
            @warn("Fail to make sample1")
            return nothing
        end
    end

    nn = 0
    y = rand(rng, censoring_dist)
    while y < tau
        y = rand(rng, censoring_dist)
        nn += 1
        if nn > 10000
            @warn("Fail to make sample2")
            return nothing
        end
    end

    if x < y
        (x, tau, true)
    else
        (y, tau, false)
    end
end

function makeSample(rng, n, dist, trunc_prob, trunc_dist, censoring_dist)
    [makeSample(rng, dist, trunc_prob, trunc_dist, censoring_dist) for i = 1:n]
end

function lognormalSample(rng; mu=1.0, sig=0.5)
    ## Log normal
    dist = LogNormal(mu, sig)
    m1 = Distributions.mean(dist)
    ires = deint(0.0, Inf64) do x
        pdf(dist, x)
    end
    @test isapprox(m1, sum(@. ires.w * ires.x) * ires.h)
    m2 = sum(@. ires.w * ires.x^2) * ires.h
    m3 = sum(@. ires.w * ires.x^3) * ires.h
    @printf("# (LogNormal(%.6f, %.6f)) moments: (m1,m2,m3) = (%.6f, %.6f, %.6f)\n", mu, sig, m1, m2, m3)

    sample = makeSample(rng, 1000, dist, 0.3, Exponential(2.0), Exponential(5.0))
    (ntruncated, ncensored, nboth) = (sum([x[2] > 0.0 for x = sample]), sum(x[3] == 0 for x = sample), sum(x[2] > 0.0 && x[3] == 0 for x = sample))
    @printf("# Sample: N=%d (truncated=%d, censored=%d, both=%d)\n", length(sample), ntruncated, ncensored, nboth)
    sample
end

function weibullSample(rng; shape=5.0, scale=8.0)
    dist = Weibull(shape, scale)
    m1 = Distributions.mean(dist)
    ires = deint(0.0, Inf64) do x
        pdf(dist, x)
    end
    @test isapprox(m1, sum(@. ires.w * ires.x) * ires.h)
    m2 = sum(@. ires.w * ires.x^2) * ires.h
    m3 = sum(@. ires.w * ires.x^3) * ires.h
    println((m1, m2, m3))

    sample = makeSample(rng, 1000, dist, 0.3, Exponential(2.0), Exponential(5.0))
    (ntruncated, ncensored, nboth) = (sum([x[2] > 0.0 for x = sample]), sum(x[3] == 0 for x = sample), sum(x[2] > 0.0 && x[3] == 0 for x = sample))
    @printf("Sample: N=%d (truncated=%d, censored=%d, both=%d)\n", length(sample), ntruncated, ncensored, nboth)
    sample
end

function mixturemodelSample(rng; mu=1.0, sig=0.5, shape=5.0, scale=8.0, p=0.3)
    dist = MixtureModel([LogNormal(mu, sig), Weibull(shape, scale)], [p, 1-p])
    m1 = Distributions.mean(dist)
    ires = deint(0.0, Inf64) do x
        pdf(dist, x)
    end
    @test isapprox(m1, sum(@. ires.w * ires.x) * ires.h)
    m2 = sum(@. ires.w * ires.x^2) * ires.h
    m3 = sum(@. ires.w * ires.x^3) * ires.h
    println((m1, m2, m3))

    sample = makeSample(rng, 1000, dist, 0.3, Exponential(2.0), Exponential(5.0))
    (ntruncated, ncensored, nboth) = (sum([x[2] > 0.0 for x = sample]), sum(x[3] == 0 for x = sample), sum(x[2] > 0.0 && x[3] == 0 for x = sample))
    @printf("Sample: N=%d (truncated=%d, censored=%d, both=%d)\n", length(sample), ntruncated, ncensored, nboth)
    sample
end

function to_data(data)
    t = [x[1] for x = data]
    tau = [x[2] for x = data]
    delta = [x[3] for x = data]
    LeftTruncRightCensoredSample(t, tau, delta)
end

### generate lognormal sample

seed = 20220810
rng = MersenneTwister(seed)

sample = lognormalSample(rng, mu=1.0, sig=0.5)
open("data/lognormal_sample.txt", "w") do f
    for (t, tau, delta) in sample
        @printf(f, "%.12f %.12f %d\n", t, tau, delta)
    end
end
dat = to_data(sample)

##### test phfitting
phfit(CF1(3), dat, verbose = false)

### generate weibull sample

seed = 20220810
rng = MersenneTwister(seed)

sample = weibullSample(rng, shape=5.0, scale=8.0)
open("data/weibull_sample.txt", "w") do f
    for (t, tau, delta) in sample
        @printf(f, "%.12f %.12f %d\n", t, tau, delta)
    end
end
dat = to_data(sample)

####### test phfitting
phfit(CF1(3), dat, verbose = false)

### generate mixture model sample
seed = 20220810
rng = MersenneTwister(seed)

sample = mixturemodelSample(rng, mu=1.0, sig=0.5, shape=5.0, scale=8.0, p=0.3)
open("data/mixturemodel_sample.txt", "w") do f
    for (t, tau, delta) in sample
        @printf(f, "%.12f %.12f %d\n", t, tau, delta)
    end
end
dat = to_data(sample)
####### test phfitting
phfit(CF1(3), dat, verbose = false)


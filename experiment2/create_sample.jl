using PhaseTypeDistributions
using PhaseTypeDistributions.Phfit
using Distributions
using Random
using Deformula
using Test
using Printf

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
            throw("Fail to make sample1 x: $(x), tau: $(tau)")
            return nothing
        end
    end

    y = rand(rng, censoring_dist)
    if x < y
        (x, tau, true)
    else
        (y, tau, false)
    end
end

function makeSample(rng, n, dist, trunc_prob, trunc_dist, censoring_dist)
    [makeSample(rng, dist, trunc_prob, trunc_dist, censoring_dist) for i = 1:n]
end

function lognormalSample(rng; n, mu, sig, p, lambda1, lambda2)
    ## Log normal
    dist = LogNormal(mu, sig)
    sample = makeSample(rng, n, dist, p, Exponential(lambda1), Exponential(lambda2))
    # println("Generated samples: ", sample)
    (ntruncated, ncensored, nboth) = (sum([x[2] > 0.0 for x = sample]), sum(x[3] == 0 for x = sample), sum(x[2] > 0.0 && x[3] == 0 for x = sample))
    @printf("# Sample: N=%d (truncated=%d, censored=%d, both=%d)\n", length(sample), ntruncated, ncensored, nboth)
    sample
end

function weibullSample(rng; n, shape, scale, p, lambda1, lambda2)
    dist = Weibull(shape, scale)
    sample = makeSample(rng, n, dist, p, Exponential(lambda1), Exponential(lambda2))
    # println("Generated samples: ", sample)
    (ntruncated, ncensored, nboth) = (sum([x[2] > 0.0 for x = sample]), sum(x[3] == 0 for x = sample), sum(x[2] > 0.0 && x[3] == 0 for x = sample))
    @printf("Sample: N=%d (truncated=%d, censored=%d, both=%d)\n", length(sample), ntruncated, ncensored, nboth)
    sample
end

function mixturemodelSample(rng; n, mu, sig, shape, scale, mixturep, p, lambda1, lambda2)
    dist = MixtureModel([LogNormal(mu, sig), Weibull(shape, scale)], [mixturep, 1-mixturep])
    sample = makeSample(rng, n, dist, p, Exponential(lambda1), Exponential(lambda2))
    # println("Generated samples: ", sample)
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

seed = 20220810
nlist = [200, 1000]

### generate lognormal sample

for n in nlist
    rng = MersenneTwister(seed)
    println("Generating lognormal sample n=$(n)")
    sample = lognormalSample(rng, n=n, mu=1.0, sig=0.5, p=0.3, lambda1=1.0, lambda2=8.0)
    open("data/lognormal_sample_$(n).txt", "w") do f
        for (t, tau, delta) in sample
            @printf(f, "%.12f %.12f %d\n", t, tau, delta)
        end
    end
    ##### test phfitting
    dat = to_data(sample)
    phfit(CF1(3), dat, verbose = false)
end


### generate weibull sample

for n in nlist
    rng = MersenneTwister(seed)
    println("Generating weibull sample n=$(n)")
    sample = weibullSample(rng, n=n, shape=5.0, scale=8.0, p=0.3, lambda1=1.0, lambda2=8.0)
    open("data/weibull_sample_$(n).txt", "w") do f
        for (t, tau, delta) in sample
            @printf(f, "%.12f %.12f %d\n", t, tau, delta)
        end
    end
    ##### test phfitting
    dat = to_data(sample)
    phfit(CF1(3), dat, verbose = false)
end

### generate mixture model sample

for n in nlist
    rng = MersenneTwister(seed)
    println("Generating mixture model sample n=$(n)")
    sample = mixturemodelSample(rng, n=n, mu=1.0, sig=0.5, shape=5.0, scale=8.0, mixturep=0.3, p=0.3, lambda1=1.0, lambda2=8.0)
    open("data/mixturemodel_sample_$(n).txt", "w") do f
        for (t, tau, delta) in sample
            @printf(f, "%.12f %.12f %d\n", t, tau, delta)
        end
    end
    ##### test phfitting
    dat = to_data(sample)
    phfit(CF1(3), dat, verbose = false)
end


using Distributions
using Random
using Printf

function makeSample(rng, dist, censoring_dist)
    x = rand(rng, dist)
    y = rand(rng, censoring_dist)
    if x < y
        (x, true)
    else
        (y, false)
    end
end

function makeSample(rng, n, dist, censoring_dist)
    [makeSample(rng, dist, censoring_dist) for i = 1:n]
end

seed = 1234
nlist = [200, 1000]

for n in nlist
    rng = MersenneTwister(seed)
    println("Generating right-censored sample n=$(n)")
    sample = makeSample(rng, n, Weibull(2.0, 1.0), Exponential(2.0))
    @printf("obs: %d, censored: %d\n", count(x -> x[2], sample), count(x -> !x[2], sample))
    open("data/right_censored_sample_$(n).txt", "w") do f
        for (t, delta) in sample
            @printf(f, "%.12f %d\n", t, delta)
        end
    end
end

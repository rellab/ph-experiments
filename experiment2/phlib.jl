
using PhaseTypeDistributions
using PhaseTypeDistributions.Phfit
using SparseArrays
using Printf
using Random
using Deformula

# ---------- I/O ----------
function read_rc_data(path::AbstractString)
    # read data in format:
    # 2.814430454718 0.000000000000 1
    # 2.473600738771 0.329502387196 1
    # 1.716287505774 0.000000000000 1
    # 3.600494144950 2.818997659998 1
    # 2.773293824760 0.000000000000 1
    data = []
    event_count = 0
    cens_count = 0
    trunc_count = 0
    both_count = 0
    open(path, "r") do f
        for x in eachline(f)
            x = strip(x)
            if !isempty(x)
                a = split(x, " ")
                if length(a) == 3
                    t = parse(Float64, a[1])
                    tau = parse(Float64, a[2])
                    delta = parse(Int, a[3])
                    push!(data, (t, tau, delta))
                    if delta == 0 && tau > 0.0
                        both_count += 1
                    elseif delta == 1 && tau > 0.0
                        trunc_count += 1
                    elseif delta == 0
                        cens_count += 1
                    elseif delta == 1
                        event_count += 1
                    else
                        error("Invalid data format")
                    end
                end
            end
        end
    end

    @printf("Loaded data: N=%d (events=%d, right-censored=%d, left-truncated=%d, both=%d)\n", length(data), event_count, cens_count, trunc_count, both_count)

    return data
    # # generate LeftTruncRightCensoredSample
    # times = [x[1] for x in data]
    # left_trunc = [x[2] for x in data]  # left truncation time (0.0 if none)
    # is_event = [x[3] == 1 for x in data]  # true = event, false = right-censored

    # dat = LeftTruncRightCensoredSample(times, left_trunc, is_event)
    # return dat
end

function calcllf(model, dat::LeftTruncRightCensoredSample)
    result = 0.0
    cumtime = 0.0
    for i in 1:dat.length
        cumtime += dat.tdat[i]
        if dat.nu[i] == 0
            tmp = log(phpdf(model, cumtime))
            # @printf("k=%d observed value = %.12f\n", i, tmp)
            result += tmp
        elseif dat.nu[i] == 1
            tmp = log(phccdf(model, cumtime))
            # @printf("k=%d censored value = %.12f\n", i, tmp)
            result += tmp
        end
    end
    return result
end

function bootstrap(rng, samples)
    n = length(samples)
    b = []
    for i = 1:n
        j = rand(rng, 1:n)
        push!(b, samples[j])
    end
    b
end

function to_data(data)
    t = [x[1] for x = data]
    tau = [x[2] for x = data]
    delta = [x[3]==1 for x = data]
    LeftTruncRightCensoredSample(t, tau, delta)
end

function cross(dist, ph)
    ires = deint(0.0, Inf64) do x
        pdf(dist, x)
    end
#     sum(ires.w .* log.(pdf.(dist, ires.x) ./ phpdf(ph, ires.x))) * ires.h
    s = 0.0
    xx = phpdf(ph, ires.x)
    for (i,w) = enumerate(ires.w)
        if xx[i] != 0.0
            s -= w * log(xx[i])
        end
    end
    s * ires.h
end

function eic(rng, ph0, llf0, d0, data; bsample, maxiter, steps, reltol)
#     b1 = Vector{Float64}(undef, bsample)
#     b2 = Vector{Float64}(undef, bsample)
#     b3 = Vector{Float64}(undef, bsample)
#     b4 = Vector{Float64}(undef, bsample)
    bias = Vector{Union{Float64,Nothing}}(undef, bsample)
    d1 = [to_data(bootstrap(rng, data)) for k = 1:bsample]
    Threads.@threads for i = 1:bsample
        try
            ph1, llf1, conv1, _ = phfit(ph0, d1[i]; initialize=false, maxiter=maxiter, steps=steps, reltol=reltol)
            if conv1 == false
                @warn("Did not convergence")
            end
            b3 = llf0
            b2 = calcllf(ph0, d1[i])
            b1 = llf1
            b4 = calcllf(ph1, d0)
            bias[i] = b1 - b2 + b3 - b4
        catch e
            bias[i] = nothing
        end
    end
    v = [x for x = bias if x !== nothing]
    mu = sum(v) / length(v)
    sigma = sqrt(sum((x - mu) .^ 2 for x in v) / (length(v) - 1))
    # return EIC and 95% CI of -2*(llf0 - V)
    return -2*(llf0 - mu), -2*(llf0 - (mu + 1.96*sigma/sqrt(length(v)))), -2*(llf0 - (mu - 1.96*sigma/sqrt(length(v))))
end

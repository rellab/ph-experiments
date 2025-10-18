
using PhaseTypeDistributions
using PhaseTypeDistributions.Phfit
using SparseArrays
using Printf

using Logging

struct ErrorLogger <: AbstractLogger end
Logging.min_enabled_level(::ErrorLogger) = Logging.Debug
Logging.shouldlog(::ErrorLogger, level, _module, group, id) = true

function Logging.handle_message(::ErrorLogger, level, message, _module, group, id, file, line; kwargs...)
    if level == Logging.Warn
        error("Warning treated as error: $message")
    else
        println("[$(level)] $message")
    end
end

global_logger(ErrorLogger())

# ---------- I/O ----------
function read_rc_data(path::AbstractString)
    # read delta time format. delta=1: event, delta=0: right-censored
    # until -1 or 1 line is found.
    data = []
    event_count = 0
    cens_count = 0
    open(path, "r") do f
        for x in eachline(f)
            x = strip(x)
            if x == "-1" || x == "1"
                break
            end
            if !isempty(x)
                a = split(x, " ")
                if length(a) == 2
                    delta = parse(Int, a[1])
                    t = parse(Float64, a[2])
                    push!(data, (t, 0.0, delta))
                    if delta == 1
                        event_count += 1
                    else
                        cens_count += 1
                    end
                end
            end
        end
    end

    @printf("Loaded data: N=%d (events=%d, right-censored=%d)\n", event_count + cens_count, event_count, cens_count)

    # generate LeftTruncRightCensoredSample
    times = [x[1] for x in data]
    left_trunc = [x[2] for x in data]  # left truncation time (0.0 if none)
    is_event = [x[3] == 1 for x in data]  # true = event, false = right-censored

    dat = LeftTruncRightCensoredSample(times, left_trunc, is_event)
    return dat
end

"Generate GPH model from init file (m x (m+1)) where 1st column is alpha and rest are Q."
function read_alpha_Q(path::AbstractString)
    al = Float64[]
    Q_rows = Vector{Float64}[]
    
    open(path, "r") do f
        for x in eachline(f)
            x = replace(x, r"^\s+" => "")  # remove leading whitespace
            if !isempty(strip(x))
                a = split(x, r"\s+")  # split by whitespace
                if !isempty(a)
                    alpha_val = parse(Float64, a[1])
                    push!(al, alpha_val)
                    
                    q_row = [parse(Float64, v) for v in a[2:end]]
                    push!(Q_rows, q_row)
                end
            end
        end
    end
    
    n = length(al)
    Q = zeros(n, n)
    for (i, row) in enumerate(Q_rows)
        for (j, val) in enumerate(row)
            Q[i, j] = val
        end
    end
    
    xi = -Q * ones(n)
    
    model = GPH(al, SparseMatrixCSC(Q), xi)
    return model
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

# ---------- one run ----------
"Run EM for fixed steps"
function run_em_once(model::GPH, dat::LeftTruncRightCensoredSample; steps::Int, eps::Float64=1e-8)
    # E-step workspace
    eres = PhaseTypeDistributions.Phfit.Estep(model)

    llf = 0.0
    t = @elapsed begin
        for k in 1:steps
            llf = PhaseTypeDistributions.Phfit.estep!(model, dat, eres, eps=eps)
            PhaseTypeDistributions.Phfit.mstep!(model, eres)
            # @printf(" EM step %3d: loglik = %.6f\n", k, llf)
        end
    end
    return t, llf, model
end

"EM until convergence"
function run_em_until(model, eres, dat::LeftTruncRightCensoredSample; steps::Int,
    abstol::Float64, reltol::Float64, maxsteps::Int, eps::Float64=1e-8)

    prev_llf = -Inf
    current_llf = 0.0
    total_steps = 0

    ## run first EMsteps
    for _ in 1:steps
        current_llf = PhaseTypeDistributions.Phfit.estep!(model, dat, eres, eps=eps)
        PhaseTypeDistributions.Phfit.mstep!(model, eres)
    end

    while total_steps < maxsteps
        total_steps += steps
        prev_llf = current_llf
        for _ in 1:steps
            current_llf = PhaseTypeDistributions.Phfit.estep!(model, dat, eres, eps=eps)
            PhaseTypeDistributions.Phfit.mstep!(model, eres)
        end
        # @printf("After %d steps, log-likelihood: %.10f\n", total_steps, current_llf)
        if abs(current_llf - prev_llf) < reltol * abs(prev_llf) || abs(current_llf - prev_llf) < abstol
            @printf("\nConverged.\n")
            return (status = "Converged", ph = model, aerror = abs(current_llf - prev_llf),
                rerror = abs(current_llf - prev_llf) / abs(prev_llf), total_steps = total_steps)
        end
    end
    @printf("\nReached max steps (%d) without convergence.\n", maxsteps)
    return (status = "MaxSteps", ph = model, aerror = abs(current_llf - prev_llf),
        rerror = abs(current_llf - prev_llf) / abs(prev_llf), total_steps = total_steps)
end

function run_fit(data_file::AbstractString, init_file::AbstractString,
    steps::Int, abstol::Float64=1e-3, reltol::Float64=1e-6, maxsteps::Int=5000)
    GC.gc()  # run garbage collection to free memory for benchmarking
    println("Starting run_fit with:", " data_file=", data_file, " init_file=", init_file)
    dat = read_rc_data(data_file)
    model = read_alpha_Q(init_file)
    eres = PhaseTypeDistributions.Phfit.Estep(model)
    m = length(model.alpha)
    t = @elapsed begin
        emresult = run_em_until(model, eres, dat; steps=steps, maxsteps=maxsteps, abstol=abstol, reltol=reltol)
    end
    llf = calcllf(emresult.ph, dat)
    result = (status=emresult.status, m=m, t=t, llf=llf, aerror = emresult.aerror, rerror = emresult.rerror, total_steps=emresult.total_steps)
    @printf("m=%d Status=%s\n", m, emresult.status)
    @printf(" elapsed=%.6f[s] final loglik=%.6f total_steps=%d\n", t, llf, emresult.total_steps)
    @printf("  aerror=%.6e rerror=%.6e\n\n\n", emresult.aerror, emresult.rerror)
    return result
end

function run_fit_cf1(data_file::AbstractString, init_file::AbstractString,
    steps::Int, abstol::Float64=1e-3, reltol::Float64=1e-6, maxsteps::Int=5000)
    GC.gc()  # run garbage collection to free memory for benchmarking
    println("Starting run_fit with:", " data_file=", data_file, " init_file=", init_file)
    dat = read_rc_data(data_file)
    model = read_alpha_Q(init_file)
    eres = PhaseTypeDistributions.Phfit.Estep(model)
    m = length(model.alpha)

    rate = Vector{Float64}(undef, length(model.alpha))
    for i in 1:length(model.alpha)
        rate[i] = -model.T[i,i]
    end
    cf1model = CF1(model.alpha, rate)

    t = @elapsed begin
        emresult = run_em_until(cf1model, eres, dat; steps=steps, maxsteps=maxsteps, abstol=abstol, reltol=reltol)
    end
    llf = calcllf(emresult.ph, dat)
    result = (status=emresult.status, m=m, t=t, llf=llf, aerror = emresult.aerror, rerror = emresult.rerror, total_steps=emresult.total_steps)
    @printf("m=%d Status=%s\n", m, emresult.status)
    @printf(" elapsed=%.6f[s] final loglik=%.6f total_steps=%d\n", t, llf, emresult.total_steps)
    @printf("  aerror=%.6e rerror=%.6e\n\n\n", emresult.aerror, emresult.rerror)
    return result
end


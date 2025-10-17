#!/usr/bin/env julia
# ph.jl — LTRC データ + 固定初期値から EM を固定回数だけ実行し、時間を計測
#
# 使用例:
#   julia ph.jl --data data/unweighted10000 --init params/phase30 --steps 100
#
# データ形式 (right-censored; left-trunc は 0 とする):
#   <delta> <time>
#   ...
#   -1            # 以降は無視（あれば）
#
# 初期値ファイル (m 行, m+1 列):
#   a1  s11 s12 ... s1m
#   a2  s21 s22 ... s2m
#   ...
#
# 依存:
#   using PhaseTypeDistributions, PhaseTypeDistributions.Phfit, SparseArrays

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

# ---------- tiny CLI ----------
function getarg(flag::String; default::Union{Nothing,String}=nothing)
    i = findfirst(==(flag), ARGS)
    if i === nothing
        return default
    end
    if i == length(ARGS)
        error("Missing value for $flag")
    end
    return ARGS[i+1]
end

function hasflag(flag::String)
    any(==(flag), ARGS)
end

# ---------- I/O ----------
function read_rc_data(path::AbstractString)
    # delta time 形式を読む。-1 行があればそこまで。
    data = []
    event_count = 0
    cens_count = 0
    open(path, "r") do f
        for x in eachline(f)
            x = strip(x)
            if x == "-1" || x == "1"  # 終了マーカー
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

    # LeftTruncRightCensoredSample を正しい方法で作成
    times = [x[1] for x in data]
    left_trunc = [x[2] for x in data]  # 左打ち切り時間（すべて0.0）
    is_event = [x[3] == 1 for x in data]  # true = イベント, false = 右打ち切り
    
    dat = LeftTruncRightCensoredSample(times, left_trunc, is_event)
    return dat
end

"init ファイル (m x (m+1)) から (alpha, Q) を読み、GPHモデルを作成する。1列目=alpha、残り=Q。"
function read_alpha_Q(path::AbstractString)
    al = Float64[]
    Q_rows = Vector{Float64}[]
    
    open(path, "r") do f
        for x in eachline(f)
            x = replace(x, r"^\s+" => "")  # 先頭の空白を除去
            if !isempty(strip(x))
                a = split(x, r"\s+")  # 空白で分割
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
    
    # xi を計算（Qの各行の合計の負値）
    xi = -Q * ones(n)
    
    # GPH モデルを作成
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
"GPHモデルとデータ dat から、EM を steps 回だけ実行。所要時間と最終 loglik を返す。"
function run_em_once(model::GPH, dat::LeftTruncRightCensoredSample; steps::Int, eps::Float64=1e-8)
    # E-step ワーク領域を確保
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

# ---------- main ----------
function main()
    data_path = getarg("--data")
    init_path = getarg("--init")
    
    if data_path === nothing || init_path === nothing
        error("Usage: julia ph.jl --data <rc_data.txt> --init <init_ph.txt> [--steps 100]")
    end
    
    steps = parse(Int, getarg("--steps", default="100"))
    llf = parse(Float64, getarg("--llf", default="1e10"))

    dat = read_rc_data(data_path)
    # # データ構造を確認
    # totalN = dat.length
    
    # # イベントと打ち切りの数を元のデータから計算
    # # 元のファイルを再読み込みして正確な数を取得
    # event_count = 0
    # cens_count = 0
    # open(data_path, "r") do f
    #     for x in eachline(f)
    #         x = strip(x)
    #         if x == "-1" || x == "1"  # 終了マーカー
    #             break
    #         end
    #         if !isempty(x)
    #             a = split(x, " ")
    #             if length(a) == 2
    #                 delta = parse(Int, a[1])
    #                 if delta == 1
    #                     event_count += 1
    #                 else
    #                     cens_count += 1
    #                 end
    #             end
    #         end
    #     end
    # end
    
    # @printf("Loaded data: N=%d (events=%d, right-censored=%d)\n", totalN, event_count, cens_count)
    @printf("EM fixed steps: %d\n", steps)

    # 単一の初期値ファイルから実行
    model = read_alpha_Q(init_path)
    m = length(model.alpha)
    @printf("Initial llf: %.6f\n", calcllf(model, dat))
    #t, llf, final_model = run_em_once(model, dat; steps=steps)
    t, llf, final_model, actual_steps = run_em_until(model, dat; targetllf=llf, maxsteps=steps)
    llf = calcllf(final_model, dat)

    # パラメータをファイルに出力
    data_name = replace(basename(data_path), r".*/" => "")
    init_name = replace(basename(init_path), r".*/" => "")
    output_file = "result_Julia_$(data_name)_$(init_name).txt"
    
    println("Saving final parameters to: $output_file")
    
    # alpha | Q の形式で出力（初期ファイルと同じ形式）
    # GPHモデルから正しくパラメータを取得
    final_alpha = final_model.alpha
    final_Q = Array(final_model.T)  # .TフィールドがQに対応、SparseMatrixCSCから密行列へ変換
    
    open(output_file, "w") do f
        for i in 1:m
            # alpha値を出力
            @printf(f, "%.16f", final_alpha[i])
            # Q行列の各行を出力
            for j in 1:m
                if final_Q[i, j] == 0.0
                    @printf(f, "  0")
                else
                    @printf(f, "  %.16f", final_Q[i, j])
                end
            end
            println(f)
        end
    end

    println("\n=== Timing (same data; fixed EM steps) ===")
    @printf("%-12s  %6s  %12s  %14s  %6s\n", "source", "m", "elapsed[s]", "final loglik", "steps")
    @printf("%-12s  %6d  %12.6f  %14.6f  %6d\n", "init-file", m, t, llf, actual_steps)
end

# main()

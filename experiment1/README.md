# Experiment 1: PH fitting with censored data (Julia vs R)

This folder reproduces Experiment 1: a comparative analysis of EM algorithms for Phase-Type (PH) distributions on right‑censored data. We run:

- Julia implementation using PhaseTypeDistributions.Phfit (proposed method)
- R implementation using the CRAN package matrixdist (UNI method; RK/PADE optional)

Both variants iterate EM until convergence based on absolute/relative log‑likelihood improvements.

## Directory layout

- `experiment1_100.jl`, `experiment1_1000.jl` – Julia batch scripts for N=100 and N=1000 datasets
- `phlib.jl` – Julia helpers: data/param readers, EM loop with convergence checks
- `experiment1_100.R`, `experiment1_1000.R` – R batch scripts for N=100 and N=1000
- `script.R` – R helpers: data/param readers, wrappers around matrixdist::fit with convergence loop
- `create_initph.R` – R script to generate initial PH parameter files under `params/`
- `Dockerfile-julia` – Julia image with PhaseTypeDistributions
- `Dockerfile-matrixdist` – R image with `matrixdist`
- `data/` – input datasets (right‑censored): `unweighted100`, `unweighted1000`, …
- `params/` – initial PH parameters for general/bidiagonal structures, various dimensions (m)
- `results/` – reference spreadsheets from prior runs (not overwritten here)
- `log_julia_*.txt`, `log_r_*.txt` – example logs from previous runs

## Quick start (Docker)

Requirements:

- Docker Desktop (macOS)

From this `experiment1/` directory:

1) Build images

```bash
docker build -f Dockerfile-julia -t ph-julia .
docker build -f Dockerfile-matrixdist -t ph-matrixdist .
```

2) Generate initial parameter files (if not already present)

```bash
mkdir -p params
docker run --rm -v "$PWD":/work -w /work ph-matrixdist Rscript generate_initial_params.R create_initph.R
```

3) Run the Julia experiments

- N = 100

```bash
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment1_100.jl | tee log_julia_100.txt
```

- N = 1000

```bash
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment1_1000.jl | tee log_julia_1000.txt
```

4) Run the R (matrixdist) experiments

- N = 100

```bash
docker run --rm -v "$PWD":/work -w /work ph-matrixdist Rscript experiment1_100.R | tee log_r_100.txt
```

- N = 1000

```bash
docker run --rm -v "$PWD":/work -w /work ph-matrixdist Rscript experiment1_1000.R | tee log_r_1000.txt
```

Notes

- All scripts print timing, final log‑likelihood, convergence status, errors, and total EM steps.
- The example uses `tee` to save a copy of the console output as a log file.

## What the scripts do

Shared settings (default; adjustable in scripts):

- EM batch steps per check: `steps = 10`
- Convergence: `abstol = 1e-3`, `reltol = 1e-5`
- Maximum steps: `maxiter = 5000`

Julia (`experiment1_*.jl` via `phlib.jl`)

- Loads right‑censored data (`data/unweighted100` or `data/unweighted1000`).
- For multiple initializations in `params/` (general and bidiagonal, various m), runs EM until convergence using PhaseTypeDistributions.Phfit.
- Also evaluates CF1 (Coxian‑form obtained from diagonal rates of the fitted GPH) for bidiagonal initializations via `run_fit_cf1`.

R (`experiment1_*.R` via `script.R`)

- Uses `matrixdist::fit` with method = `"UNI"` (changeable to `"RK"` or `"PADE"`).
- Runs EM in blocks of `steps`, checks convergence, and reports timing/log‑likelihood.
- Two tolerance settings are evaluated for UNI (1e‑6 and 1e‑8) in the provided scripts.

## Input formats

Right‑censored data (`data/unweighted100`, `data/unweighted1000`)

```
<delta> <time>
...
-1           # optional terminator (lines after are ignored)
```

- `delta = 1` observed event at `time`
- `delta = 0` right‑censored at `time`
- Julia treats left‑truncation as 0 for all rows (LeftTruncRightCensoredSample).

Initial PH parameters (`params/*.txt`)

- m rows × (m+1) columns in the form: `alpha | S` (or `Q`) per row
	- Column 1: initial probability `alpha_i`
	- Columns 2..(m+1): row i of the transient generator matrix (S/Q)
- Example (m = 2):

```
0.90   -2.00   2.00
0.10    0.00  -2.02
```

`phlib.jl` reads these files and builds a GPH model. The R path (`script.R`) builds a PH with the same convention.

## Customization

- Change EM settings at the top of `experiment1_*.jl` and `experiment1_*.R`:
	- `steps`, `abstol`, `reltol`, `maxiter`
- To try different dimensions/structures, point to other files under `params/`.
- For matrixdist, switch method by editing `method <- "UNI"` to `"RK"` or `"PADE"`.

## Troubleshooting

- If a run stops with warnings promoted to errors (Julia), `phlib.jl` sets a custom logger that treats warnings as errors to surface numerical issues. Inspect the log for the offending step/initialization.

## Reproducibility notes

- No RNG is used in these batch scripts; results depend only on the data and 
	initial parameter files.
- Console logs in this folder (`log_*.txt`) show example outputs from prior runs.


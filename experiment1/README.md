# Experiment 1: PH fitting with censored data (Julia vs R)

This folder reproduces Experiment 1: a comparative analysis of EM algorithms for Phase-Type (PH) distributions on right‑censored data. We run:

- Julia implementation using PhaseTypeDistributions.Phfit (proposed method)
- R implementation using the CRAN package matrixdist (UNI method; RK/PADE optional)

Both variants iterate EM until convergence based on absolute/relative log‑likelihood improvements.

## Directory layout

- `experiment1.jl` – Julia scripts
- `phlib.jl` – Julia helpers: data/param readers, EM loop with convergence checks
- `experiment1.R` – R scripts
- `script.R` – R helpers: data/param readers, wrappers around matrixdist::fit with convergence loop
- `create_sample.jl` – Julia script to generate right-censored sample data under `data/`
- `create_initph.R` – R script to generate initial PH parameter files under `params/`
- `Dockerfile-julia` – Julia image with PhaseTypeDistributions
- `Dockerfile-matrixdist` – R image with `matrixdist`
- `data/` – input datasets (right‑censored): `unweighted100`, `unweighted1000`, …
- `params/` – initial PH parameters for general/bidiagonal structures, various dimensions (m)
- `julia_200.txt`, `julia_1000.txt` – Julia experiment logs
- `r_200.txt`, `r_1000.txt` – R experiment logs

## Quick start (Docker)

Requirements:

- Docker Desktop (macOS, probably Windows and Linux run similarly)

From this `experiment1/` directory:

1) Build images

```bash
docker build -f Dockerfile-julia -t ph-julia .
docker build -f Dockerfile-matrixdist -t ph-matrixdist .
```

2) Generate data files (if not already present)

```bash
mkdir -p data
docker run --rm -v "$PWD":/work -w /work ph-julia julia create_sample.jl
```

3) Generate initial parameter files (if not already present)

```bash
mkdir -p params
docker run --rm -v "$PWD":/work -w /work matrixdist Rscript create_initph.R
```

4) Run the Julia experiments

```bash
caffeinate -i docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment1.jl 200 | tee julia_200.txt
caffeinate -i docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment1.jl 1000 | tee julia_1000.txt
```

5) Run the R (matrixdist) experiments

```bash
caffeinate -i docker run --rm -v "$PWD":/work -w /work matrixdist Rscript experiment1.R 200 | tee r_200.txt
caffeinate -i docker run --rm -v "$PWD":/work -w /work matrixdist Rscript experiment1.R 1000 | tee r_1000.txt
```

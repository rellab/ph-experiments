# Experiment 2: Demonstrate EIC model selection

This experiment showcases the use of the EIC (Extended Information Criterion) for model selection in the context of phase-type distributions. The goal is to compare the performance of different phase-type models in capturing the characteristics of the underlying data.

## Files

- `Dockerfile-julia`: Dockerfile for building the Julia environment.
- `phlib.jl`: Library of phase-type functions.
- `create_sample.jl`: Script for generating synthetic data.
- `experiment2-1.jl`: Script for computing AIC and EIC.
- `experiment2-2.jl`: Script for computing moments and cross entropy.
- `experiment2-3.jl`: Script for drawing plots.

## Steps

1) Build the Julia Docker image

```bash
docker build -f Dockerfile-julia -t ph-julia .
```
2) Run the script to generate synthetic data

```bash
mkdir -p data
docker run --rm -v "$PWD":/work -w /work ph-julia julia create_sample.jl
```

3) Run the experiment to compute AIC and EIC

```bash
mkdir -p results
caffeinate -i docker run --rm -v "$PWD":/work -w /work ph-julia julia -t 4 experiment2-1.jl 200 | tee eic200.txt
caffeinate -i docker run --rm -v "$PWD":/work -w /work ph-julia julia -t 4 experiment2-1.jl 1000 | tee eic1000.txt
```

4) Run the experiment to compute moments and cross entropy

```bash
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment2-2.jl 200 | tee mom200.txt
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment2-2.jl 1000 | tee mom1000.txt
```

5) Run the script to draw plots

- Modify the experiment2-3_XXX.jl script to indicate the selected models by AIC, EIC, and cross entropy.

```bash
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment2-3.jl 200 20 30 50 ...
docker run --rm -v "$PWD":/work -w /work ph-julia julia experiment2-3.jl 1000 10 20 30 ...
```

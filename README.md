# Experiments for Phase-Type Distributions

This workspace contains two experiments for fitting phase-type (PH) distributions to right‑censored data.

Summary

- Experiment 0 — PH distribution estimation using EM algorithm (EMpht: C implementation)
	- Location: `experiment0/` (see `experiment0/README.md` for full details).

- Experiment 1 — Comparative study (Julia vs R)
	- Purpose: compare EM implementations (Julia `PhaseTypeDistributions.Phfit` vs R `matrixdist`).
	- Location: `experiment1/` (see `experiment1/README.md` for full details).

- Experiment 2 — EIC model selection (Julia only)
	- Purpose: compute EIC (Extended Information Criterion), compare moments, and produce plots.
	- Location: `experiment2/` (see `experiment2/README.md` for full details).

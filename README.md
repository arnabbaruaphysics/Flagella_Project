# Flagella_Project
 Regulatory memory and growth-coupled inheritance shape nutrient-dependent flagella number variation in  Salmonella
# Regulatory memory and growth-coupled inheritance shape nutrient-dependent flagella number variation in *Salmonella*

Code and data accompanying the manuscript by **Arnab Barua<sup>†</sup>, Maria Giralt-Zuñiga<sup>†</sup>, Marc Erhardt, and Haralampos Hatzikirou** (<sup>†</sup> equal contribution).

**Affiliations:** Khalifa University (Mathematics Department, Abu Dhabi); TU Dresden (Center for Information Services and High Performance Computing); Humboldt-Universität zu Berlin (Institute of Biology / Molecular Microbiology); Max Planck Unit for the Science of Pathogens, Berlin.

---

## Overview

The number of flagella per cell in *Salmonella* varies from cell to cell and depends on nutrient availability. This repository contains a Jupyter notebook that:

1. Holds the raw experimental data: per-cell flagella-count histograms and OD600 growth curves for wild type (WT) and the Δ*rflP* mutant.
2. Fits a **mean-field model of flagella density** coupled to a latent regulatory-memory variable ρ(t) and to growth-driven dilution.
3. Quantifies fit quality, parameter identifiability, and sloppiness (Fisher Information Matrix).
4. Analyses the **signal-to-noise ratio (SNR)** of flagella-number variability and the **information-geometric speed limit** on flagellar reprogramming.

The notebook reproduces the data-analysis figures of the paper (Figs. 3–9 plus SI identifiability analysis).

## Model summary

The mean flagella count per cell, `X(t)`, relaxes toward a target `k` at a rate set by growth:

```
dX/dt = (μ0 − a·X) · (k − X)
k     = K0 + c·Γ(ρ),    Γ(ρ) = 1/√(1−ρ²) − 1  ≈ ρ²/2   (small ρ)
```

where ρ(t) is a latent correlation (regulatory memory) that decays on a timescale `τ_ρ`:

| Strain | ρ(t) |
|---|---|
| Δ*rflP* mutant | `ρ0_mut · exp(−t/τ_ρ)` |
| WT | `ρ* + (ρ0_wt − ρ*) · exp(−t/τ_ρ)` |

The notebook uses a **frozen-rate** analytical approximation, `λ = μ0 − a·k` (WT: `k = K0 + Γρ*²`; mutant: `k = K0`), which gives closed-form trajectories. These are cross-checked against the full nonlinear ODE (RK45) throughout ("hybrid" analytical + numerical mode).

**Fixed constants:** `a = 0.000265` (per-flagellum growth-cost coefficient), `c = 100` (sensing coupling; only `c·Γ(ρ)` is identifiable).

**Fitted parameters**

- *Global (shared across conditions):* `τ_ρ`, `ρ0_wt`, `ρ0_mut`, `X0_wt`, `X0_mut`
- *Per WT concentration:* `ρ*`, `K0` (with `μ0` taken from OD600 growth data)
- *Per mutant concentration:* `K0`

## Experimental conditions

| | WT | Δ*rflP* |
|---|---|---|
| Nutrient concentrations | 0.2%, 0.4%, 0.5%, 1% | 0.2%, 1% |
| Flagella-count time points (min) | 60, 90, 120, 150, 180, 240 (0.4%: up to 180) | 60, 120, 180, 240 |
| Biological replicates (flagella counts) | 3 (R1–R3) | 2 (R1–R2) |
| OD600 growth curves | 3 replicates, 0–240 min (0.2, 0.5, 1%) | 2 replicates, 0–240 min (0.2, 1%) |

Notes:
- The 0.4% WT condition has no 240-min OD data, so it is not in the growth plot. Its μ0 and fit are handled separately, and its parameters are predicted/interpolated from the other WT concentrations where indicated in the figures.
- Mutant standard errors are computed across replicates as `σ_repl/√2` (WT: `σ_repl/√3`).

## Notebook contents

| Section / cell | Description | Figure |
|---|---|---|
| Shared model | Single source of truth for all model functions (analytical + ODE). **Run first.** | – |
| Raw flagella distributions | Histograms per concentration/time point (replicate-averaged), WT and Δ*rflP* | Fig. 3 |
| Growth curves | OD600 vs time, WT vs Δ*rflP* | Fig. 4 |
| Model fitting | Hybrid analytical + numerical fit (5-D outer optimisation with warm starts + basin-hopping; inner Nelder–Mead per condition) | – |
| Mean-field dynamics | Fitted trajectories with goodness-of-fit | Fig. 5 |
| Goodness-of-fit table | N, dof, MSE, RMSE, χ², reduced χ² (exports CSV/PNG/PDF) | Table |
| ρ(t) and K(ρ) | Latent correlation and target flagella number over time, with bootstrap 95% bands | Fig. 6 |
| Parameters vs nutrient | Fitted parameters vs concentration; SNR vs τ_relax | Figs. 7–8a |
| SNR vs time | SNR = ⟨x̃⟩²/Var(x̃), with x̃ = ln(1+x) | Fig. 8b |
| Identifiability | Pairwise log₁₀(MSE) cost contours over (ρ*/ρ0, K0, μ0, X0) | SI |
| Sloppiness | Global Fisher Information Matrix for the 16-parameter model | Fig. 9 |
| Speed limit | Cramér–Rao bound on \|d⟨x⟩/dt\|, information length per condition, per-interval saturation ratios | Fig. 7 / 8 (speed limit) |

## Getting started

### Requirements

- Python ≥ 3.9
- `numpy`, `scipy`, `matplotlib`, `jupyter`

```bash
pip install numpy scipy matplotlib jupyter
```

### Run

```bash
git clone <this-repo-url>
cd <repo-name>
jupyter notebook Regulatory_memory_flagella_SHARED_lambda_mu0_ak.ipynb
```

Run cells **in order from the top**:

1. The shared-model cell defines all model functions, constants, and the output location.
2. The fitting cell (`run_simulations`) is the slow step. It writes `sim_cache_hybrid.pkl`, which all downstream analysis cells load.
3. Remaining cells can then be run individually.

Iteration counts in the fit and bootstrap are deliberately reduced for speed (`N_STARTS_INNER=10`, `N_BH_STEPS=2`, `N_WARM_STARTS=5`, `N_BOOT=300`). Increase them for final results.

The notebook can also be run on Google Colab.

### Output location

All figures and tables are written to the directory given by the `FLAG_OUT_DIR` environment variable (default: `./outputs`):

```bash
export FLAG_OUT_DIR=/path/to/results
```

> **Note:** the last cell (SI identifiability, standalone version) has hardcoded Colab paths (`/content/outputs/...`). Change them to `OUT_DIR` / `CACHE` if you run locally.

### Outputs

`growth_WT_vs_DrflP.png`, `flagella_distribution.png`, `fit_quality_table.{csv,png,pdf}`, `rho_and_K_over_time.{png,pdf}`, `params_vs_conc.{png,pdf}`, `snr_vs_tau_t240.{png,pdf}`, `snr_vs_time.{png,pdf}`, `cost_contours_all_conc.{png,pdf}`, `global_fim_16.{png,pdf}`, `figure7_speed_limit.{png,pdf}`, and the fit cache `sim_cache_hybrid.pkl`.

## Repository structure

```
.
├── Regulatory_memory_flagella_SHARED_lambda_mu0_ak.ipynb   # all data, model, analysis
├── README.md
└── outputs/                                                # created on first run
```

## Data

All raw data (flagella-count histograms per replicate, concentration, and time point; OD600 growth curves) is embedded directly in the notebook as Python dictionaries, so no external data files are required.

## Citation

If you use this code or data, please cite:

> Barua A., Giralt-Zuñiga M., Erhardt M., Hatzikirou H. *Regulatory memory and growth-coupled inheritance shape nutrient-dependent flagella number variation in Salmonella.* (manuscript; add journal/DOI/preprint link here)

## License

*Add a license (e.g. MIT, CC-BY-4.0) before publishing.*

## Contact

- Arnab Barua: Khalifa University
- Haralampos Hatzikirou: Khalifa University / TU Dresden
- Marc Erhardt: Humboldt-Universität zu Berlin

*(Add email addresses or GitHub handles as appropriate.)*

# Simulation Study: Power and FDR of Bayes Factor Tests for Spatial Nonstationarity

## Overview

This simulation study estimates the **false discovery rate (FDR)** and **average power** associated with the decision rule that rejects stationarity in favour of nonstationarity when BF₁₀ > c. Two simulation settings mirror actual data analyses:

| Simulation | Dataset | Genes | Non-stationary (%) | Key characteristic |
|:----------:|---------|------:|-------------------:|-------------------|
| 1 | Human Breast Cancer (ILC) | 16,564 | 49.99% | High nonstationarity prevalence |
| 2 | Human Ovarian Cancer | 16,684 | 5.32% | Low nonstationarity prevalence |

## Simulation Design

For each dataset, genes are classified based on the actual analysis results:

- **H0 genes** (BF₁₀ ≤ 1): Simulated from a stationary Matérn process using the gene's estimated parameters (τ, κ, σ²_ε).
- **H1 genes** (BF₁₀ > 1): Simulated from a nonstationary Matérn process (Ingebrigtsen et al. 2014) where log(κ(s)) = θ₂ + θ₃·x̃(s) + θ₄·ỹ(s).

### H1 Simulation: Same Model, Same Parameters When Possible

**Ideal case (fitted thetas available):** When the main analysis has been re-run with theta extraction, we use the fitted coefficients θ₁, θ₂, θ₃, θ₄ from the nonstationary INLA fit. Simulation is done via `inla.spde2.precision()` and `inla.qsample()`—the **same** SPDE model and **same** basis-function coefficients as the fitted model. So we simulate from the exact nonstationary process that was estimated for that gene.

**Fallback (fitted thetas not available):** When results lack theta columns (e.g. older runs), we assign synthetic θ₃, θ₄ with magnitude varying by BF strength: Anecdotal (0.15–0.40), Moderate (0.30–0.70), Strong (0.55–1.00), Very strong (0.80–1.50), Extreme (1.20–3.00). This gives a range of nonstationarity from low to high across H1 genes.

**Why INLA for simulation (not rSPDE/BayesNSGP):** We use the INLA SPDE framework for simulation because it is the same model used for inference. rSPDE and BayesNSGP implement different nonstationary formulations; using INLA ensures model consistency and direct use of fitted parameters.

Each replicate simulates spatial expression for all genes, fits both stationary and nonstationary SPDE models via INLA, and computes BF₁₀. The full study runs 100 replicates per dataset.

### Outputs

- **FDR(c)** = E[FP / (FP + TP)] — averaged over 100 replicates
- **Average power(c)** = mean over H1 genes of Pr(BF₁₀ > c | H1) — averaged over 100 replicates

## File Structure

```
simulation/
├── simulation_functions.R      # Core functions: simulation, fitting, FDR/power
├── prepare_simulation.R        # One-time setup: extract params, create config
├── run_replicate.R             # Run one simulation replicate (SLURM array)
├── aggregate_results.R         # Compute FDR/power curves from replicates
├── plot_fdr_power.R            # Create publication figure
├── prepare_simulation.sh        # Login node: prepare configs
├── aggregate_simulation.sh     # Login node: aggregate and plot
├── submit_simulation.sbatch    # SLURM: run replicates (array=1-100)
├── configs/                    # [generated] simulation config RDS files
└── output/                     # [generated] replicate results and curves
    ├── HumanBreastCancerILC/   # per-replicate RDS files
    ├── HumanOvarianCancer/     # per-replicate RDS files
    ├── *_fdr_power_curves.rds  # aggregated FDR/power curves
    └── fdr_power_figure.pdf    # publication figure
```

## Running on HPC (Alliance Canada)

### Prerequisites

Ensure the main analysis has been completed (with theta extraction for H1 genes; re-run if using older results):
- `results/HumanBreastCancerILC_results.rds` exists
- `results/HumanOvarianCancer_results.rds` exists
- `data/spatial_datasets/HumanBreastCancerILC.rds` exists
- `data/spatial_datasets/HumanOvarianCancer.rds` exists
- INLA is installed (`bash hpc/install_inla.sh`)

### Step 1: Prepare Simulation Configurations

Run on the **login node** (lightweight):

```bash
bash simulation/prepare_simulation.sh
```

This creates `simulation/configs/HumanBreastCancerILC_sim_config.rds` and `simulation/configs/HumanOvarianCancer_sim_config.rds`.

### Step 2: Run Simulation Replicates

After Step 1 completes:

```bash
# Launch 100 replicates for each dataset
sbatch --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch
sbatch --export=DATASET=HumanOvarianCancer   simulation/submit_simulation.sbatch
```

Each replicate runs as an independent SLURM array task (~2–3 hours on 8 cores). All 200 tasks can run in parallel if resources are available.

To test with a small subset first:

```bash
sbatch --array=1-3 --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch
```

### Step 3: Aggregate and Plot

After all 100 replicates complete for both datasets, run on the **login node**:

```bash
bash simulation/aggregate_simulation.sh
```

This produces `simulation/output/fdr_power_figure.pdf` and `fdr_power_figure.png`.

## Running Locally (Testing)

A callable function `run_local_simulation_test()` runs the full pipeline end-to-end. Use it from R or via the CLI script.

**From R (interactive):**
```r
source("simulation/run_local_test.R")
run_local_simulation_test(n_genes = 10, n_reps = 2)
run_local_simulation_test(10, 2, "HumanOvarianCancer")
# Returns list with curves, replicate_results, fdr_at_c1, power_at_c1, output_files
```

**From command line:**
```bash
# Default: 20 genes, 3 replicates, HumanBreastCancerILC (~5-15 min)
Rscript simulation/test_simulation_local.R

# Custom: 10 genes, 2 replicates, HumanOvarianCancer
Rscript simulation/test_simulation_local.R 10 2 HumanOvarianCancer
```

Output: `simulation/output/test_spatial_expression.png`, `test_fdr_power.png`, `test_gene_bf_distributions.png`

For a larger local test using the full pipeline scripts:

```bash
Rscript simulation/prepare_simulation.R HumanBreastCancerILC
for r in 1 2 3; do
  Rscript simulation/run_replicate.R HumanBreastCancerILC $r 4 \
      simulation/configs simulation/output 100
done
Rscript simulation/aggregate_results.R HumanBreastCancerILC 3
```

## Computational Requirements

| Resource | Per replicate | Total (100 reps × 2 datasets) |
|----------|--------------|-------------------------------|
| CPU cores | 12 | 1200 task-slots |
| Memory | 96 GB | — |
| Wall time | ~1.5–2 hours | ~1.5–2 hours (parallel) |
| Storage | ~2 MB/replicate | ~400 MB total |
| INLA fits | ~33,000 | ~3.3 million |

**Time estimate** (from local test with reduced INLA output): 30 genes in ~1.6 min (1 core) → full replicate (16,564 genes) ≈ 90 min on 12 cores.

## Key Quantities of Interest

- **FDR at c = 1**: The false discovery rate when using BF₁₀ > 1 as the decision rule
- **Average power at c = 1**: The average probability of detecting nonstationarity when it exists
- **FDR and power curves**: Full characterisation of the FDR-power tradeoff across thresholds

# Separate nonstationarity and anisotropy comparisons

This **additional** workflow refits three spatial models for each eligible gene
and tissue section. It does not replace the original Phase I/II pipeline or
write to its results, figures or manuscript files.

| Model | Spatial covariance specification | Comparison |
|---|---|---|
| `S` | Stationary isotropic Matérn | Shared reference |
| `NS` | Nonstationary, locally isotropic Matérn | `log_bf_NS_vs_S = log_ml_NS - log_ml_S` |
| `A` | Stationary anisotropic Matérn | `log_bf_A_vs_S = log_ml_A - log_ml_S` |

Both alternatives are fitted for **every** eligible gene, without a Phase I
gate. Positive natural-log Bayes factors favor the named alternative; negative
values favor stationary isotropy. These are two separate evidence comparisons,
not mutually exclusive classes or a decomposition into independent biological
effects. Both alternatives can be favored. In particular, this does not fit a
joint nonstationary-anisotropic model or test one feature conditional on the
other. Bayes factors quantify model evidence, not the magnitude of an effect.

## Cluster commands

Run from the repository root, after pulling the code. The input objects are the
12 `.rds` files listed in `hpc/analysis_datasets.txt`. Existing download scripts
can supply them if necessary.

```bash
mkdir -p logs
export RUN_DIR="$PWD/hpc/results/covariance_comparison_v1"
sbatch --export=ALL,RUN_DIR="$RUN_DIR",GENES_PER_TASK=100 \
  hpc/covariance_comparison/prepare.sbatch
```

Wait for preparation to finish successfully, inspect `units.csv` and
`task_count.txt` in that directory, then submit the gene arrays:

```bash
bash hpc/covariance_comparison/submit.sh "$RUN_DIR"
```

The wrapper submits arrays of at most 1,000 tasks (one task per gene chunk),
with a dependent combine job that runs only after all arrays succeed. Slurm
array indices are converted to global task IDs using `CC_TASK_OFFSET`.
Default resources are 1 CPU / 8 GB / 12 hours per task, with **32 concurrent
tasks per array**. For example, three arrays can run up to 96 tasks in total.
Override resources with `sbatch` for manual submissions or edit the `.sbatch`
resource directives before launching. The preparation job uses 32 GB; the
combine job uses 16 GB. No workflow command submits a job until invoked by you.

```bash
# For a lower per-array concurrency or a cluster with a smaller array limit:
CC_MAX_CONCURRENT=16 CC_ARRAY_BATCH_SIZE=500 \
  bash hpc/covariance_comparison/submit.sh "$RUN_DIR"

squeue -u "$USER"
```

Optional preparation environment variables are `DATA_DIR`, `DATASET_LIST`
(one dataset basename per line), and `GENES_PER_TASK`. Set `GENES_PER_TASK=1`
for one gene per task, at the cost of many more scheduler tasks. A chunk of 100
genes amortizes R/INLA startup and data loading. All three fits for a gene stay
in the same task and reuse the stationary fit for the two Bayes factors.

The module setup in `environment.sh` follows the existing Alliance HPC scripts
(`StdEnv/2023`, R 4.4 and Bioconductor 3.20). It also needs working `INLA` and
`fmesher` packages. `CC_R_LIBS` can override the user package library. INLA,
OpenMP and BLAS are restricted to one thread per task; parallelism comes from
the Slurm array, not nested workers. Prepare and fit on the same cluster
software stack: exact R/package versions are recorded and checked.

## Data and preprocessing

- Keep only spots marked `in_tissue` when that annotation exists. Objects
  without it retain their supplied spots; `units.csv` records this distinction.
  Unexpected annotation values cause an error rather than silently keeping
  off-tissue spots.
- Analyze different `sample_id` values separately. The anterior and posterior
  sagittal mouse objects each contain two serial sections; they must not be
  fitted as one overlapping coordinate system. The current 12 input objects
  therefore produce 14 section-level units. Results retain both the source
  dataset name and section ID; do not pool these into a 12-row summary without
  an explicit aggregation decision.
- Within each section, retain genes expressed in at least `ceiling(0.005*n)`
  spots and with at least 3 total counts. Preserve unique feature IDs even if
  display symbols are duplicated.
- Apply `log1p(counts)`, average-rank normal scores
  `qnorm((rank - 0.5)/n)`, and quadratic spatial detrending with intercept,
  centered `x`, `y`, `x^2`, `y^2`, `x*y`. There is no library-size normalization
  and no percentile/display threshold in fitting. Degenerate residuals are
  marked failed for all models, never counted as evidence for stationarity.
- The three models use the same residual vector, observation projection and
  section-specific mesh. Mesh edges and cutoff follow the existing analysis:
  `max.edge = c(max_extent/5, max_extent/3)`, `cutoff = min_extent/20`.

Do not combine these newly fitted comparisons with old Phase I results: the
run uses explicit shared priors and section-specific preparation throughout.

## Models and priors

All models have a Gaussian observation likelihood with an explicit intercept
and a Matérn field with SPDE order `alpha=2` (smoothness `nu=1` in 2D).

The stationary isotropic SPDE is

```text
(kappa^2 - Laplacian) (tau w) = spatial white noise.
```

`NS` retains `H=I` and the existing nonstationary specification:

```text
log(tau)      = theta1
log(kappa(s)) = theta2 + theta3 * scaled_x(s) + theta4 * scaled_y(s).
```

Coordinates in this basis are centered and standardized across mesh vertices.
This allows local range and marginal variance to vary; it does not introduce
an explicit directional parameter.

`A` keeps `tau`, `kappa` and the diffusion matrix `H` constant across space:

```text
(kappa^2 - div(H grad)) (tau w) = spatial white noise
H = expm([[v1, v2], [v2, -v1]]),   det(H) = 1.
```

With `r = sqrt(v1^2 + v2^2)`, the principal correlation-axis ratio is
`exp(r)` and the major-axis angle is `atan2(v2,v1)/2` modulo 180 degrees.
The angle is measured in the supplied x/y coordinates (not necessarily an
image's upward-positive coordinates). It is undefined at isotropy. Fixing
`det(H)=1` separates anisotropy from the overall range/variance scale.

The anisotropic finite-element precision is

```text
Q = tau^2 [kappa^4 C + 2 kappa^2 G(H) + G(H) C^-1 G(H)],
```

where `C` is the lumped mass matrix and `G(H)` is assembled from triangle basis
gradients on the same fixed mesh. INLA's `rgeneric` interface integrates over
the anisotropy parameters; this is not a plug-in maximization over directions.
The full Gaussian normalizing constant is retained. At `v1=v2=0`, this
precision agrees numerically with INLA's native stationary isotropic SPDE.

Priors are proper and explicit:

| Parameter | Prior |
|---|---|
| `log(tau)`, baseline `log(kappa)` | Independent normal, precision 0.1; the same mesh-dependent nominal INLA prior means in all three models |
| Nonstationarity slopes `theta3`, `theta4` | Independent normal, mean 0, SD `sqrt(10)` by default |
| Anisotropy `v1`, `v2` | Independent normal, mean 0, SD `log(4)/sqrt(-2*log(0.05))` by default |
| Observation precision | Gamma(shape 1, rate 0.00005) |
| Explicit intercept | Normal(mean 0, precision 0.001) |

The default anisotropy prior implies `Pr(axis ratio > 4) = 0.05` and a uniform
orientation. Set `CC_AXIS_RATIO_UPPER` and `CC_AXIS_RATIO_TAIL_PROB` before
preparation to change that calibration; `CC_SLOPE_SD` controls the two NS
slopes. Choose a **new run directory** for each prior-sensitivity run. Bayes
factors depend on these priors; the two alternatives are not automatically
equally flexible just because both add two parameters. Sensitivity checks
(for example, axis-ratio upper bounds 2 and 8 and alternative slope SDs) and
simulation calibration should precede paper-level claims or FDR claims.

All three fits use INLA's Gaussian latent-field strategy and CCD
hyperparameter integration. `log_ml_*` is `fit$mlik[1,1]`, the integrated
approximate log marginal likelihood. The code does not use fitted log
likelihoods or an information criterion as a substitute for Bayes factors.

Implementation references:
[INLA rgeneric tutorial](https://www.inla.r-inla-download.org/r-inla.org/doc/vignettes/rgeneric.html),
[fmesher finite-element documentation](https://inlabru-org.github.io/fmesher/reference/fm_fem.html),
and [Lindgren et al. (2011), the SPDE approach](https://doi.org/10.1111/j.1467-9868.2011.00777.x).

## Outputs and restart behavior

Preparation writes `manifest.rds`, `units.csv`, `tasks.csv`, `task_count.txt`
and one prepared `.rds` per section. The manifest records data/code hashes,
software versions, priors, mesh-derived prior means and ordered gene IDs.
Do not change the implementation, prepared files or software mid-run.

Each task writes only `chunks/task_XXXXXX.rds`, through an atomic rename. A
per-task directory lock prevents concurrent writers. A successful chunk has
exactly one row per scheduled gene, including any failed models. Its rows are
validated before writing and again before reuse/combination. Missing or
corrupt chunks, changed implementation, mismatched software and wrong gene
identities stop the workflow; incomplete data are not silently summarized.

The combine job writes:

- `comparisons.csv`: per-gene/section identifiers, both log Bayes factors,
  all three log marginal likelihoods, model-specific status/error/warnings,
  fit times, spatial hyperparameter means and anisotropy summaries.
- `comparison_summary.csv`: separate NS-vs-S and A-vs-S counts per section,
  valid and unavailable denominators, failed-model counts, and numbers above
  BF 1, 3 and 10. BF 1 is only a direction-of-preference threshold, not a
  significance or FDR threshold.
- `comparisons.rds`: the same rows and summaries plus the full run manifest.

`axis_ratio_at_posterior_mean` and `axis_angle_degrees_at_posterior_mean`
are transformations of the posterior means of `v1,v2`, **not** posterior means
or credible intervals of the ratio/angle themselves. The model-specific
`theta1_*` and `theta2_*` columns are posterior means of `log(tau)` and the
baseline `log(kappa)`. `theta3_NS`, `theta4_NS` are the two range-gradient
coefficients; `theta3_A`, `theta4_A` are `v1`, `v2`.

An S failure makes both comparisons unavailable. An NS failure leaves A-vs-S
usable if A and S succeeded, and conversely. Failed models have `NA` marginal
likelihoods/Bayes factors plus diagnostics; they are never negative evidence.
Nonzero optimizer status, non-finite output and reported convergence failure
are treated as failures, while later models/genes continue.

Rerunning `submit.sh` validates and skips completed chunks. To retry chunks
containing failed models, after previous jobs have stopped:

```bash
CC_RETRY_FAILED=1 bash hpc/covariance_comparison/submit.sh "$RUN_DIR"
```

A failed chunk is refitted in full. The old chunk is replaced only after the
new chunk passes validation. Successful chunks are not refitted. If a job was
killed, inspect its scheduler status before manually removing its specific
stale `.rds.lock` directory. Never remove a lock for a running task. Recombine
after any retries; old combined CSVs do not update until combination succeeds.

For a targeted global task retry and manual recombination, from the repo root:

```bash
sbatch --array=17 --export=ALL,RUN_DIR="$RUN_DIR",CC_TASK_OFFSET=0,CC_RETRY_FAILED=1 \
  hpc/covariance_comparison/array.sbatch
# After all required tasks finish:
sbatch --export=ALL,RUN_DIR="$RUN_DIR" hpc/covariance_comparison/combine.sbatch
```

## Tests and direct execution

```bash
Rscript hpc/covariance_comparison/test_workflow.R
Rscript hpc/covariance_comparison/test_workflow.R --fit
bash hpc/covariance_comparison/test_submit.sh

Rscript hpc/covariance_comparison/prepare.R data/spatial_datasets NEW_RUN_DIR 100
Rscript hpc/covariance_comparison/run_task.R NEW_RUN_DIR 1
Rscript hpc/covariance_comparison/combine.R NEW_RUN_DIR
```

The tests check determinant-one anisotropy, axes, an independent fmesher FEM
calculation, the isotropic precision limit, proper prior densities, independent
model failures, spot filtering, serial-section splitting, chunk validation,
restarts and separate valid-comparison denominators. `--fit` additionally runs
S, NS, A and an isotropic rgeneric model, checking that the latter reproduces
the native S marginal log likelihood within numerical tolerance. The shell
test uses a mock `sbatch` and submits no real jobs.

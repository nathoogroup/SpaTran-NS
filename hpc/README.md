# HPC workflow for the spatial nonstationarity analysis

These scripts run one SLURM array task per study dataset, retain an explicit
row for every attempted gene, and combine results only after every dataset job
has completed successfully.

## Files

| File | Purpose |
|---|---|
| `analysis_datasets.txt` | Stable array order for the 12 study datasets |
| `download_datasets.R` | Download SpatialExperiment resources from ExperimentHub |
| `analysis_functions.R` | INLA model fitting, safe per-gene rows, and integrity helpers |
| `run_analysis.R` | Analyze one array-selected dataset or all datasets sequentially |
| `combine_results.R` | Validate and combine completed per-dataset outputs |
| `submit_download.sbatch` | Download job |
| `submit_analysis.sbatch` | 12-task analysis array |
| `submit_combine.sbatch` | Post-array validation/combine job |
| `preflight_check.sh` | Modules, packages, files, syntax, and array-map checks |
| `test_result_handling.R` | Lightweight regression test for failed-gene handling |

## Failure and completeness semantics

An expected gene-level INLA failure, including Newton-Raphson
non-convergence, produces one row with model statistics set to `NA` and the
diagnostic in `error_message`. Later genes continue to run.

Before a per-dataset file is saved, the script requires:

- exactly one row per scheduled filtered gene;
- ordered, unique `gene_index` values;
- matching Ensembl `gene_id` values; and
- the required result/error columns.

An unexpected worker error or a completeness mismatch stops the R process and
returns a nonzero SLURM status. Results are written through a temporary file,
so an old result is replaced only after the complete new object has serialized.
Legacy or partial result files are not silently skipped; reusable outputs must
carry the current schema-v4 completion metadata plus fingerprints of the source
dataset and the exact analysis implementation.

Array tasks write only `Dataset_results.rds`. They do not write the shared
`all_results_combined.*` files, which avoids cross-task overwrite races.

## Cluster usage

Create the log directory before calling `sbatch`; Slurm opens its output files
before the submitted script begins.

```bash
mkdir -p logs
bash hpc/preflight_check.sh
```

Download datasets if needed:

```bash
sbatch hpc/submit_download.sbatch
```

Submit the full analysis and a dependent combine job:

```bash
analysis_job=$(sbatch --parsable hpc/submit_analysis.sbatch)
sbatch --dependency="afterok:${analysis_job}" hpc/submit_combine.sbatch
```

To refit only the two datasets affected by explicit off-tissue spots, submit
the dedicated dependency chain:

```bash
bash hpc/launch_in_tissue_refit.sh
```

This launches array tasks 11 and 12 with 32 gene-level workers each. A dependent
job then validates and migrates only the metadata for the other ten datasets
(whose source objects have no `in_tissue` field and therefore retain identical
spots), validates all 12 schema-v4 files, and rebuilds the combined result.

The full array order is version-controlled in `analysis_datasets.txt`; adding
or removing files from the data directory does not shift task indices.
Each selected SpatialExperiment is restricted to `in_tissue == 1` when that
annotation is present. Objects without an `in_tissue` field are treated as
already tissue-restricted and retain all supplied spots. This gives 3,639 spots
for `Visium_humanDLPFC` and 2,702 spots for `Visium_mouseCoronal`.
`MouseBrainCoronal`, a duplicate representation of the same coronal section,
is excluded so that this biological sample is counted once.

### Schema-v4 refresh and historical truncation indices

This change records and validates supplied, retained, and removed spot counts.
Use `launch_in_tissue_refit.sh` for the minimal refresh: only DLPFC and mouse
coronal are numerically refitted, while the dependent finalization job validates
and migrates the unchanged results for the other ten datasets. Alternatively,
the full array below recomputes all 12 datasets from scratch:

```bash
mkdir -p logs
sbatch --export=ALL,OVERWRITE=1 hpc/submit_analysis.sbatch
```

The three datasets affected by the historical truncation bug have stable
indices:

| Index | Dataset | Historical rows saved | Expected filtered genes |
|---:|---|---:|---:|
| 3 | HumanGlioblastoma | 13,368 | 17,963 |
| 5 | HumanLymphNode | 4,406 | 18,295 |
| 8 | MouseBrainSagittalAnterior | 1,880 | 16,431 |

After a complete schema-v4 result set exists, those indices can be used for a
targeted retry if any of the three tasks needs to be replaced:

```bash
mkdir -p logs
sbatch --array=3,5,8 --export=ALL,OVERWRITE=1 hpc/submit_analysis.sbatch
```

All schema-v3 files predate the in-tissue selection metadata and cannot be mixed
with schema-v4 outputs.

## Direct usage

```text
Rscript hpc/run_analysis.R [data_dir] [output_dir] [n_cores]

  data_dir   Directory containing the listed .rds datasets
  output_dir Directory for result files
  n_cores    Parallel workers (default: detected cores)
```

With no `SLURM_ARRAY_TASK_ID`, all datasets are processed sequentially and the
validated per-dataset results are also combined. In array mode the task ID
selects one entry from `analysis_datasets.txt`.

Set `OVERWRITE=1` to force a rerun. Without it, an existing output is reused
only if it passes the strict completeness check against the filtered source
dataset.

## Output columns

In addition to model estimates and Bayes factors, each per-gene row contains:

| Column | Meaning |
|---|---|
| `gene_index` | Unique position in the filtered dataset |
| `gene_id` | Stable feature/Ensembl identifier |
| `gene_name` | Display symbol; not assumed unique |
| `theta1_ns`–`theta4_ns` | Nonstationary SPDE hyperparameters |
| `error_message` | `NA` on success; diagnostic text on a failed fit |
| `dataset_name` | Source dataset basename |

Positive `log_bayes_factor` favors the nonstationary model; negative values
favor the stationary model.

## Monitoring

```bash
squeue -u "$USER"
tail -f logs/spatran_ns_JOBID_TASKID.out
```

If a task fails, its existing result file is left untouched. Fix the reported
cause and resubmit that task index; do not treat a nonempty legacy `.rds` file
as proof of completion.

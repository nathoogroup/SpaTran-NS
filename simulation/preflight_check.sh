#!/bin/bash
# =============================================================================
# Preflight Check Script for SpaTran-NS Simulation Study
# =============================================================================
# Run on a login node BEFORE submitting simulation jobs.
# Reports PASS / FAIL for each check.
#
# Prerequisites (must be completed before running this check):
#   1. bash hpc/install_inla.sh           (once, on login node)
#   2. sbatch hpc/submit_download.sbatch  (download datasets)
#   3. sbatch hpc/submit_analysis.sbatch  (run main analysis)
#
# Usage:  bash simulation/preflight_check.sh
# =============================================================================

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PASS=0; FAIL=0

ok()   { echo "[PASS] $*"; PASS=$((PASS+1)); }
fail() { echo "[FAIL] $*"; FAIL=$((FAIL+1)); }
info() { echo "[INFO] $*"; }
warn() { echo "[WARN] $*"; }

echo "========================================================"
echo " SpaTran-NS Simulation Preflight Check"
echo " Project: $PROJECT_DIR"
echo " Date:    $(date)"
echo "========================================================"
echo ""

# ---- 1. Core module environment ---------------------------------------------
echo "--- 1. Module environment ---"

if module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 \
              geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
              gsl/2.7 jags/4.3.2 \
              r-bundle-bioconductor/3.20 2>/dev/null; then
  ok "All modules loaded (r/4.4.0 + r-bundle-bioconductor/3.20)"
else
  fail "Failed to load one or more modules"
fi

if [[ -n "${EBVERSIONR:-}" ]]; then
  ok "R version: $EBVERSIONR"
else
  fail "EBVERSIONR not set — module load may have failed"
fi

# ---- 2. System R packages ---------------------------------------------------
echo ""
echo "--- 2. System R packages (r-bundle-bioconductor/3.20) ---"

check_sys_pkg() {
  local pkg=$1
  R --quiet --no-save -e "stopifnot(requireNamespace('$pkg', quietly=TRUE))" \
    >/dev/null 2>&1 \
    && ok "System package: $pkg" \
    || fail "System package MISSING: $pkg"
}

for pkg in SpatialExperiment SummarizedExperiment Matrix dplyr doParallel foreach; do
  check_sys_pkg "$pkg"
done

# ---- 3. User R library (INLA) -----------------------------------------------
echo ""
echo "--- 3. User R library (INLA) ---"

R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

if [[ -d "$R_LIBS_USER" ]]; then
  ok "User R_LIBS exists: $R_LIBS_USER"
else
  fail "User R_LIBS missing: $R_LIBS_USER  →  run: bash hpc/install_inla.sh"
fi

export R_LIBS="$R_LIBS_USER"

R --quiet --no-save -e "stopifnot(requireNamespace('INLA', quietly=TRUE))" \
  >/dev/null 2>&1 \
  && ok "INLA package found in user library" \
  || fail "INLA not installed  →  run: bash hpc/install_inla.sh"

# ---- 4. INLA binary ---------------------------------------------------------
echo ""
echo "--- 4. INLA binary ---"

R --quiet --no-save -e "
  library(INLA)
  INLA::inla.mesh.1d(seq(0,1,length.out=5))
" >/dev/null 2>&1 \
  && ok "INLA binary executes correctly" \
  || fail "INLA binary failed  →  run: bash hpc/install_inla.sh"

# ---- 5. Simulation script files ---------------------------------------------
echo ""
echo "--- 5. Simulation script files ---"

check_file() { [[ -f "$1" ]] && ok "File: $(basename $1)" || fail "MISSING: $1"; }

check_file "$PROJECT_DIR/simulation/simulation_functions.R"
check_file "$PROJECT_DIR/simulation/prepare_simulation.R"
check_file "$PROJECT_DIR/simulation/run_replicate.R"
check_file "$PROJECT_DIR/simulation/aggregate_results.R"
check_file "$PROJECT_DIR/simulation/submit_prepare.sbatch"
check_file "$PROJECT_DIR/simulation/submit_simulation.sbatch"
check_file "$PROJECT_DIR/simulation/submit_aggregate.sbatch"

check_file "$PROJECT_DIR/hpc/analysis_functions.R"

# ---- 6. SBATCH script content -----------------------------------------------
echo ""
echo "--- 6. SBATCH script content ---"

for sbatch in "$PROJECT_DIR/simulation/submit_prepare.sbatch" \
              "$PROJECT_DIR/simulation/submit_simulation.sbatch" \
              "$PROJECT_DIR/simulation/submit_aggregate.sbatch"; do
  bname=$(basename "$sbatch")
  errs=0

  grep -q "account=def-nathoo"    "$sbatch" || { fail "$bname: missing --account=def-nathoo"; errs=$((errs+1)); }
  grep -q "partition="            "$sbatch" && { fail "$bname: has --partition (remove it)";   errs=$((errs+1)); }
  grep -q "r-bundle-bioconductor" "$sbatch" || { fail "$bname: missing r-bundle-bioconductor module"; errs=$((errs+1)); }

  [[ $errs -eq 0 ]] && ok "SBATCH ok: $bname"
done

# ---- 7. Main analysis results (required by prepare step) --------------------
echo ""
echo "--- 7. Main analysis results ---"

RESULTS_DIR="$PROJECT_DIR/results"
DATASETS=("HumanBreastCancerILC" "HumanOvarianCancer")

if [[ -d "$RESULTS_DIR" ]]; then
  ok "Results directory exists: $RESULTS_DIR"
else
  fail "Results directory missing  →  run main analysis first (sbatch hpc/submit_analysis.sbatch)"
fi

for ds in "${DATASETS[@]}"; do
  rds="$RESULTS_DIR/${ds}_results.rds"
  if [[ -f "$rds" ]]; then
    sz=$(du -h "$rds" | cut -f1)
    ok "Results: ${ds}_results.rds ($sz)"
  else
    fail "Results MISSING: ${ds}_results.rds  →  main analysis not complete for $ds"
  fi
done

# ---- 8. Spatial datasets (required by prepare step) -------------------------
echo ""
echo "--- 8. Spatial datasets ---"

DATA_DIR="$PROJECT_DIR/data/spatial_datasets"

for ds in "${DATASETS[@]}"; do
  rds="$DATA_DIR/${ds}.rds"
  if [[ -f "$rds" ]]; then
    sz=$(du -h "$rds" | cut -f1)
    ok "Dataset: ${ds}.rds ($sz)"
  else
    fail "Dataset MISSING: ${ds}.rds  →  sbatch hpc/submit_download.sbatch"
  fi
done

# ---- 9. Output directories --------------------------------------------------
echo ""
echo "--- 9. Output directories ---"

check_or_create_dir() {
  if [[ -d "$1" ]]; then
    ok "Dir: $1"
  else
    mkdir -p "$1" 2>/dev/null \
      && ok "Dir: $1 (created)" \
      || fail "Cannot create: $1"
  fi
}

check_or_create_dir "$PROJECT_DIR/logs"
check_or_create_dir "$PROJECT_DIR/simulation/configs"
check_or_create_dir "$PROJECT_DIR/simulation/output"

# ---- 10. Existing simulation state ------------------------------------------
echo ""
echo "--- 10. Existing simulation state ---"

CONFIG_DIR="$PROJECT_DIR/simulation/configs"
OUTPUT_DIR="$PROJECT_DIR/simulation/output"

for ds in "${DATASETS[@]}"; do
  cfg="$CONFIG_DIR/${ds}_sim_config.rds"
  if [[ -f "$cfg" ]]; then
    sz=$(du -h "$cfg" | cut -f1)
    info "Config exists: ${ds}_sim_config.rds ($sz) — prepare step already done"
  else
    info "No config for $ds — run: sbatch simulation/submit_prepare.sbatch"
  fi

  rep_dir="$OUTPUT_DIR/$ds"
  if [[ -d "$rep_dir" ]]; then
    n_reps=$(find "$rep_dir" -maxdepth 1 -name "replicate_*.rds" 2>/dev/null | wc -l)
    if [[ $n_reps -gt 0 ]]; then
      info "$ds: $n_reps replicate(s) completed"
    else
      info "$ds: replicate directory exists but no results yet"
    fi
  fi
done

# ---- Summary ----------------------------------------------------------------
echo ""
echo "========================================================"
echo " Results: $PASS PASSED  |  $FAIL FAILED"
echo "========================================================"

if [[ $FAIL -eq 0 ]]; then
  echo ""
  echo " All checks passed. Submit simulation jobs in order:"
  echo ""
  echo "   1. sbatch simulation/submit_prepare.sbatch"
  echo "   2. sbatch --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch"
  echo "      sbatch --export=DATASET=HumanOvarianCancer   simulation/submit_simulation.sbatch"
  echo "   3. sbatch simulation/submit_aggregate.sbatch"
  echo "   4. Download results and plot locally"
else
  echo ""
  echo " Fix the FAILED items above, then re-run this check."
fi
echo ""
exit $FAIL

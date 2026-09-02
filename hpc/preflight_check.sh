#!/bin/bash
set -uo pipefail

# =============================================================================
# Preflight Check Script for SpaTran-NS Cluster Analysis
# =============================================================================
# Run on a login node BEFORE submitting jobs.
# Reports PASS / FAIL for each check.
#
# Package strategy:
#   Download job : r/4.4.0 + r-bundle-bioconductor/3.20 (system packages only)
#   Analysis job : r/4.4.0 + r-bundle-bioconductor/3.20 + user R_LIBS (INLA)
#
# Usage:  bash hpc/preflight_check.sh
# =============================================================================

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PASS=0; FAIL=0

ok()   { echo "[PASS] $*"; PASS=$((PASS+1)); }
fail() { echo "[FAIL] $*"; FAIL=$((FAIL+1)); }
info() { echo "[INFO] $*"; }

echo "========================================================"
echo " SpaTran-NS Preflight Check"
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

# ---- 2. System R packages (from r-bundle-bioconductor) ----------------------
echo ""
echo "--- 2. System R packages (r-bundle-bioconductor/3.20) ---"

check_sys_pkg() {
  local pkg=$1
  R --quiet --no-save -e "stopifnot(requireNamespace('$pkg', quietly=TRUE))" \
    >/dev/null 2>&1 \
    && ok "System package: $pkg" \
    || fail "System package MISSING: $pkg"
}

for pkg in SpatialExperiment SummarizedExperiment ExperimentHub \
           Matrix dplyr doParallel foreach; do
  check_sys_pkg "$pkg"
done

# ---- 3. User R library (INLA) -----------------------------------------------
echo ""
echo "--- 3. User R library (INLA from install_inla.sh) ---"

R_VERSION_SHORT="${EBVERSIONR:-}"
R_LIBS_USER=""
if [[ -n "$R_VERSION_SHORT" ]]; then
  R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/${R_VERSION_SHORT:0:3}"
fi

if [[ -n "$R_LIBS_USER" && -d "$R_LIBS_USER" ]]; then
  ok "User R_LIBS exists: $R_LIBS_USER"
else
  fail "User R_LIBS missing or R version unavailable  →  run: bash hpc/install_inla.sh"
fi

if [[ -n "$R_LIBS_USER" ]]; then
  export R_LIBS="$R_LIBS_USER"
fi

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

# ---- 5. Project file structure ----------------------------------------------
echo ""
echo "--- 5. Project file structure ---"

check_file() { [[ -f "$1" ]] && ok "File: $1" || fail "MISSING file: $1"; }
check_dir()  { [[ -d "$1" ]] && ok "Dir:  $1" || fail "MISSING dir:  $1 (created at runtime)"; }

check_file "$PROJECT_DIR/hpc/analysis_functions.R"
check_file "$PROJECT_DIR/hpc/analysis_datasets.txt"
check_file "$PROJECT_DIR/hpc/combine_results.R"
check_file "$PROJECT_DIR/hpc/download_datasets.R"
check_file "$PROJECT_DIR/hpc/run_analysis.R"
check_file "$PROJECT_DIR/hpc/install_inla.sh"
check_file "$PROJECT_DIR/hpc/submit_download.sbatch"
check_file "$PROJECT_DIR/hpc/submit_analysis.sbatch"
check_file "$PROJECT_DIR/hpc/submit_combine.sbatch"
check_file "$PROJECT_DIR/hpc/test_result_handling.R"
check_file "$PROJECT_DIR/SpaTran-NS.Rmd"

check_dir  "$PROJECT_DIR/data/spatial_datasets"
check_dir  "$PROJECT_DIR/results"
check_dir  "$PROJECT_DIR/logs"

# ---- 6. SBATCH script content -----------------------------------------------
echo ""
echo "--- 6. SBATCH script content ---"

for sbatch in "$PROJECT_DIR/hpc/submit_download.sbatch" \
              "$PROJECT_DIR/hpc/submit_analysis.sbatch" \
              "$PROJECT_DIR/hpc/submit_combine.sbatch"; do
  bname=$(basename "$sbatch")
  errs=0

  grep -q "account=def-nathoo"          "$sbatch" || { fail "$bname: missing --account=def-nathoo"; errs=$((errs+1)); }
  grep -q "partition="                  "$sbatch" && { fail "$bname: has --partition (remove it)";   errs=$((errs+1)); }
  grep -q "r-bundle-bioconductor"       "$sbatch" || { fail "$bname: missing r-bundle-bioconductor module"; errs=$((errs+1)); }

  [[ $errs -eq 0 ]] && ok "SBATCH ok: $bname"
done

# ---- 7. Script syntax --------------------------------------------------------
echo ""
echo "--- 7. Script syntax ---"

for shell_file in "$PROJECT_DIR/hpc/preflight_check.sh" \
                  "$PROJECT_DIR/hpc/install_inla.sh" \
                  "$PROJECT_DIR/hpc/submit_download.sbatch" \
                  "$PROJECT_DIR/hpc/submit_analysis.sbatch" \
                  "$PROJECT_DIR/hpc/submit_combine.sbatch"; do
  if bash -n "$shell_file"; then
    ok "Shell syntax: $(basename "$shell_file")"
  else
    fail "Shell syntax error: $(basename "$shell_file")"
  fi
done

for r_file in "$PROJECT_DIR/hpc/analysis_functions.R" \
              "$PROJECT_DIR/hpc/combine_results.R" \
              "$PROJECT_DIR/hpc/download_datasets.R" \
              "$PROJECT_DIR/hpc/run_analysis.R" \
              "$PROJECT_DIR/hpc/test_result_handling.R"; do
  if Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$r_file" \
      >/dev/null 2>&1; then
    ok "R syntax: $(basename "$r_file")"
  else
    fail "R syntax error: $(basename "$r_file")"
  fi
done

if Rscript "$PROJECT_DIR/hpc/test_result_handling.R" >/dev/null 2>&1; then
  ok "Per-gene failure/completeness regression test"
else
  fail "Per-gene failure/completeness regression test"
fi

# ---- 8. Downloaded datasets and array mapping -------------------------------
echo ""
echo "--- 8. Downloaded Visium datasets and array mapping ---"

DATA_DIR="$PROJECT_DIR/data/spatial_datasets"
DATASET_LIST="$PROJECT_DIR/hpc/analysis_datasets.txt"
if [[ -d "$DATA_DIR" && -f "$DATASET_LIST" ]]; then
  N=$(sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "$DATASET_LIST" | wc -l)
  if [[ $N -gt 0 ]]; then
    idx=0
    missing=0
    while IFS= read -r dataset_name; do
      idx=$((idx+1))
      dataset_file="$DATA_DIR/${dataset_name}.rds"
      if [[ -f "$dataset_file" ]]; then
        sz=$(du -h "$dataset_file" | cut -f1)
        info "  [$idx] ${dataset_name}.rds  ($sz)"
      else
        fail "Missing listed dataset: $dataset_file"
        missing=$((missing+1))
      fi
    done < <(sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "$DATASET_LIST")

    [[ $missing -eq 0 ]] && ok "All $N study datasets are present"

    ARRAY_MAX=$(sed -n \
      's/^#SBATCH --array=1-\([0-9][0-9]*\).*$/\1/p' \
      "$PROJECT_DIR/hpc/submit_analysis.sbatch" | head -n 1)
    if [[ -n "$ARRAY_MAX" && "$ARRAY_MAX" -eq "$N" ]]; then
      ok "SBATCH array 1-$ARRAY_MAX matches analysis_datasets.txt"
    else
      fail "SBATCH array bound (${ARRAY_MAX:-not found}) does not match $N listed datasets"
    fi
  else
    fail "No datasets are listed in $DATASET_LIST"
  fi
else
  [[ -d "$DATA_DIR" ]] || fail "Data directory missing  →  mkdir -p $DATA_DIR"
  [[ -f "$DATASET_LIST" ]] || fail "Dataset list missing: $DATASET_LIST"
fi

# ---- Summary ----------------------------------------------------------------
echo ""
echo "========================================================"
echo " Results: $PASS PASSED  |  $FAIL FAILED"
echo "========================================================"

if [[ $FAIL -eq 0 ]]; then
  echo " All checks passed. Ready to submit jobs."
else
  echo " Fix the FAILED items above, then re-run this check."
fi
echo ""
exit $FAIL

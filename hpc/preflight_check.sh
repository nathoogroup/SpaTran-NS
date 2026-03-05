#!/bin/bash
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

if module load StdEnv/2023 gcc/12.3 r/4.4.0 \
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

# ---- 5. Project file structure ----------------------------------------------
echo ""
echo "--- 5. Project file structure ---"

check_file() { [[ -f "$1" ]] && ok "File: $1" || fail "MISSING file: $1"; }
check_dir()  { [[ -d "$1" ]] && ok "Dir:  $1" || fail "MISSING dir:  $1 (created at runtime)"; }

check_file "$PROJECT_DIR/hpc/analysis_functions.R"
check_file "$PROJECT_DIR/hpc/download_datasets.R"
check_file "$PROJECT_DIR/hpc/run_analysis.R"
check_file "$PROJECT_DIR/hpc/install_inla.sh"
check_file "$PROJECT_DIR/hpc/submit_download.sbatch"
check_file "$PROJECT_DIR/hpc/submit_analysis.sbatch"
check_file "$PROJECT_DIR/SpaTran-NS.Rmd"

check_dir  "$PROJECT_DIR/data/spatial_datasets"
check_dir  "$PROJECT_DIR/results"
check_dir  "$PROJECT_DIR/logs"

# ---- 6. SBATCH script content -----------------------------------------------
echo ""
echo "--- 6. SBATCH script content ---"

for sbatch in "$PROJECT_DIR/hpc/submit_download.sbatch" \
              "$PROJECT_DIR/hpc/submit_analysis.sbatch"; do
  bname=$(basename "$sbatch")
  errs=0

  grep -q "account=def-nathoo"          "$sbatch" || { fail "$bname: missing --account=def-nathoo"; errs=$((errs+1)); }
  grep -q "partition="                  "$sbatch" && { fail "$bname: has --partition (remove it)";   errs=$((errs+1)); }
  grep -q "r-bundle-bioconductor"       "$sbatch" || { fail "$bname: missing r-bundle-bioconductor module"; errs=$((errs+1)); }

  [[ $errs -eq 0 ]] && ok "SBATCH ok: $bname"
done

# ---- 7. Downloaded datasets -------------------------------------------------
echo ""
echo "--- 7. Downloaded Visium datasets ---"

DATA_DIR="$PROJECT_DIR/data/spatial_datasets"
if [[ -d "$DATA_DIR" ]]; then
  N=$(find "$DATA_DIR" -maxdepth 1 -name "*.rds" \
      ! -name "download_log*" ! -name "manifest*" 2>/dev/null | wc -l)
  if [[ $N -gt 0 ]]; then
    ok "$N dataset .rds files found"
    find "$DATA_DIR" -maxdepth 1 -name "*.rds" \
         ! -name "download_log*" ! -name "manifest*" \
      | sort | while read f; do
          sz=$(du -h "$f" | cut -f1)
          info "  $(basename $f)  ($sz)"
        done
  else
    fail "No dataset files yet  →  sbatch hpc/submit_download.sbatch"
  fi
else
  fail "Data directory missing  →  mkdir -p $DATA_DIR"
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

#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/cov_submit_test.XXXXXX")
run_dir="$test_dir/run with spaces"
mkdir "$run_dir"
# Generated scheduler-test fixtures, not real prepared analysis data.
touch "$run_dir/manifest.rds"
printf '2005\n' > "$run_dir/task_count.txt"
export CC_TEST_LOG="$test_dir/submissions.log"
export PATH="$script_dir/test-bin:$PATH"
CC_MAX_CONCURRENT=7 bash "$script_dir/submit.sh" "$run_dir"
log=$(<"$CC_TEST_LOG")
[[ "$log" == *'--array=1-1000%7'*'CC_TASK_OFFSET=0'* ]]
[[ "$log" == *'--array=1-1000%7'*'CC_TASK_OFFSET=1000'* ]]
[[ "$log" == *'--array=1-5%7'*'CC_TASK_OFFSET=2000'* ]]
[[ "$log" == *'--dependency=afterok:9001:9002:9003'* ]]
printf '0\n' > "$run_dir/task_count.txt"
if bash "$script_dir/submit.sh" "$run_dir"; then
  echo "Invalid task count was accepted" >&2; exit 1
fi
[[ "$(<"$CC_TEST_LOG")" == "$log" ]]
printf '2005\n' > "$run_dir/task_count.txt"
export CC_TEST_LOG="$test_dir/failed_submissions.log"
if CC_TEST_FAIL_CALL=2 bash "$script_dir/submit.sh" "$run_dir"; then
  echo "Scheduler failure was ignored" >&2; exit 1
fi
[[ "$(<"$CC_TEST_LOG")" != *'--dependency='* ]]
echo "Slurm wrapper tests passed (no real jobs submitted). Fixtures: $test_dir"

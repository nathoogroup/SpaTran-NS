#!/usr/bin/env bash
# Submit prepared, gene-chunk arrays and an afterok combine job. No R required.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: bash hpc/covariance_comparison/submit.sh PREPARED_RUN_DIR" >&2
  exit 2
fi
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(cd "$script_dir/../.." && pwd -P)
run_dir=$(cd "$1" && pwd -P)
if [[ "$run_dir" == *','* || "$run_dir" == *$'\n'* ]]; then
  echo "RUN_DIR must not contain commas or newlines (Slurm --export syntax)." >&2
  exit 2
fi
[[ -f "$run_dir/manifest.rds" && -f "$run_dir/task_count.txt" ]] || {
  echo "Preparation is incomplete: $run_dir" >&2; exit 2;
}
read -r task_count < "$run_dir/task_count.txt"
batch_size=${CC_ARRAY_BATCH_SIZE:-1000}
concurrent=${CC_MAX_CONCURRENT:-32}
for value in "$task_count" "$batch_size" "$concurrent"; do
  if [[ ! "$value" =~ ^[1-9][0-9]*$ || ${#value} -gt 7 ]]; then
    echo "Task count, batch size and concurrency must be positive integers <= 9999999." >&2
    exit 2
  fi
done
command -v sbatch >/dev/null || { echo "sbatch is not available." >&2; exit 2; }
mkdir -p "$repo_dir/logs"
job_ids=()
failure() {
  echo "Submission stopped. Already submitted array jobs: ${job_ids[*]:-none}" >&2
  echo "Inspect those jobs before resubmitting; no jobs were cancelled." >&2
}
trap failure ERR
offset=0
while (( offset < task_count )); do
  count=$((task_count-offset))
  if (( count > batch_size )); then count=$batch_size; fi
  job=$(sbatch --parsable --chdir="$repo_dir" \
    --array="1-${count}%${concurrent}" \
    --export="ALL,RUN_DIR=$run_dir,CC_TASK_OFFSET=$offset" "$script_dir/array.sbatch")
  job=${job%%;*}
  [[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch job ID: $job" >&2; false; }
  job_ids+=("$job")
  echo "Array $job: global tasks $((offset+1))-$((offset+count))"
  offset=$((offset+count))
done
dependencies=$(IFS=:; echo "${job_ids[*]}")
combine=$(sbatch --parsable --chdir="$repo_dir" --dependency="afterok:$dependencies" \
  --export="ALL,RUN_DIR=$run_dir" "$script_dir/combine.sbatch")
echo "Combine job ${combine%%;*}; results: $run_dir/comparisons.csv"
echo "Concurrency cap is $concurrent PER array (${#job_ids[@]} arrays)."

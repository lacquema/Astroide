#!/bin/bash

set -euo pipefail

usage() {
  cat <<'EOF'
Benchmark one Astroide sub-simulation with several OpenMP thread counts.

Usage:
  benchmark_threads.sh --case-dir PATH --binary PATH [options]

Required arguments:
  --case-dir PATH      Directory containing files_XX.in, options_XX.in, etc.
  --binary PATH        Executable to benchmark.

Optional arguments:
  --index N            Sub-simulation index to benchmark (default: 1).
  --threads "LIST"     Space-separated or comma-separated thread counts.
                       Default: "1 4 8".
  --tstop VALUE        Temporary final time used in the benchmark copy.
                       Default: 2000.
  --stacksize VALUE    Value exported as STACKSIZE. Default: 1000000.
  --temp-parent PATH   Parent directory used for temporary benchmark runs.
                       Default: /tmp.
  --keep-temp          Keep temporary benchmark directories after the run.
  -h, --help           Show this help message.

Notes:
  - The script never edits the original simulation directory.
  - It copies the selected sub-simulation to a temporary directory, rewrites
    the output path, and shortens the run by replacing t_stop with --tstop.
  - Output and dump frequencies are set to t_stop in the benchmark copy to
    reduce I/O noise during the comparison.

Examples:
  benchmark_threads.sh \
    --case-dir /path/to/test_dyna3 \
    --binary /path/to/swift_rmvs3_par

  benchmark_threads.sh \
    --case-dir /path/to/test_dyna2 \
    --binary /path/to/swift_hjs_par \
    --threads "1 2 4 8" \
    --tstop 5000
EOF
}

error() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

case_dir=''
binary=''
index=1
threads='1 4 8'
tstop='2000'
stacksize='1000000'
temp_parent='/tmp'
keep_temp='false'

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir)
      case_dir=${2:-}
      shift 2
      ;;
    --binary)
      binary=${2:-}
      shift 2
      ;;
    --index)
      index=${2:-}
      shift 2
      ;;
    --threads)
      threads=${2:-}
      shift 2
      ;;
    --tstop)
      tstop=${2:-}
      shift 2
      ;;
    --stacksize)
      stacksize=${2:-}
      shift 2
      ;;
    --temp-parent)
      temp_parent=${2:-}
      shift 2
      ;;
    --keep-temp)
      keep_temp='true'
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      error "unknown option: $1"
      ;;
  esac
done

[[ -n "$case_dir" ]] || error '--case-dir is required'
[[ -n "$binary" ]] || error '--binary is required'
[[ -d "$case_dir" ]] || error "case directory not found: $case_dir"
[[ -f "$binary" ]] || error "binary not found: $binary"
[[ -d "$temp_parent" ]] || error "temporary parent directory not found: $temp_parent"
[[ "$index" =~ ^[0-9]+$ ]] || error '--index must be an integer'

printf -v suffix '%02d' "$index"
files_name="files_${suffix}.in"

case_dir=$(cd "$case_dir" && pwd -P)
binary=$(cd "$(dirname "$binary")" && pwd -P)/$(basename "$binary")

[[ -f "$case_dir/$files_name" ]] || error "missing input file: $case_dir/$files_name"

# Keep this compatible with the Bash 3.2 shipped by macOS.
generic_name=$(sed -n '1p' "$case_dir/$files_name")
options_name=$(sed -n '2p' "$case_dir/$files_name")
bodies_name=$(sed -n '3p' "$case_dir/$files_name")
debris_name=$(sed -n '4p' "$case_dir/$files_name")

[[ -n "$generic_name" && -n "$options_name" && -n "$bodies_name" && -n "$debris_name" ]] \
  || error "unexpected format in $case_dir/$files_name"

for required_file in "$generic_name" "$options_name" "$bodies_name" "$debris_name"; do
  [[ -f "$case_dir/$required_file" ]] || error "missing case file: $case_dir/$required_file"
done

thread_list=${threads//,/ }
bench_roots=()
cleanup() {
  if [[ "$keep_temp" == 'false' ]] && [[ ${#bench_roots[@]} -gt 0 ]]; then
    rm -rf "${bench_roots[@]}"
  fi
}
trap cleanup EXIT

printf 'Case: %s\n' "$case_dir"
printf 'Binary: %s\n' "$binary"
printf 'Sub-simulation: %s\n' "$suffix"
printf 'Thread counts: %s\n' "$thread_list"
printf 'Benchmark t_stop: %s\n' "$tstop"

for thread_count in $thread_list; do
  [[ "$thread_count" =~ ^[0-9]+$ ]] || error "invalid thread count: $thread_count"

  bench_root=$(mktemp -d "$temp_parent/astroide_bench_${thread_count}.XXXXXX")
  bench_roots+=("$bench_root")
  work_dir="$bench_root/case"
  run_rel="case/run_${suffix}"

  mkdir -p "$work_dir"
  cp "$case_dir/$generic_name" "$work_dir/"
  cp "$case_dir/$bodies_name" "$work_dir/"
  cp "$case_dir/$debris_name" "$work_dir/"
  cp "$case_dir/$files_name" "$work_dir/"

  awk -v benchmark_tstop="$tstop" -v root="$bench_root" -v rel="$run_rel" '
    NR == 1 {
      printf "%s %s %s\n", $1, benchmark_tstop, $3
      next
    }
    NR == 2 {
      printf "%s %s\n", benchmark_tstop, benchmark_tstop
      next
    }
    NR == 5 {
      print root
      next
    }
    NR == 6 {
      print rel
      next
    }
    {
      print
    }
  ' "$case_dir/$options_name" > "$work_dir/$options_name"

  pushd "$work_dir" >/dev/null
  ulimit -Ss "$(ulimit -Hs)"
  export STACKSIZE="$stacksize"
  export OMP_NUM_THREADS="$thread_count"

  if /usr/bin/time -p -o "$bench_root/time.txt" \
      "$binary" < "./$files_name" \
      >"$bench_root/stdout.txt" 2>"$bench_root/stderr.txt"; then
    printf '  %s threads: ' "$thread_count"
    awk '
      /^real / { real_time = $2 }
      /^user / { user_time = $2 }
      /^sys /  { sys_time = $2 }
      END {
        printf "real=%ss user=%ss sys=%ss\n", real_time, user_time, sys_time
      }
    ' "$bench_root/time.txt"
  else
    status=$?
    printf '  %s threads: failed (exit=%s)\n' "$thread_count" "$status" >&2
    printf '    Temporary directory kept in %s\n' "$bench_root" >&2
    tail -n 20 "$bench_root/stderr.txt" >&2 || true
    keep_temp='true'
    exit "$status"
  fi

  popd >/dev/null
done

if [[ "$keep_temp" == 'true' ]]; then
  printf 'Temporary benchmark directories kept under %s\n' "$temp_parent"
fi
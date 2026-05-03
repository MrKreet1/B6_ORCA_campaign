#!/usr/bin/env bash
set -euo pipefail

# Последовательный запуск ORCA input-файлов.
# На VPS 8 CPU, если каждый расчёт использует %pal nprocs 8, НЕ запускай несколько задач параллельно.
# Использование:
#   bash scripts/run_all.sh calculations/stage1
#   bash scripts/run_all.sh calculations/final

ROOT_REL="${1:-calculations/stage1}"
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ROOT="$PROJECT_DIR/$ROOT_REL"
ORCA_CMD="${ORCA_CMD:-orca}"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/run_all_$(date +%Y%m%d_%H%M%S).log"

if ! command -v "$ORCA_CMD" >/dev/null 2>&1 && [[ ! -x "$ORCA_CMD" ]]; then
  echo "ERROR: ORCA command not found: $ORCA_CMD" | tee -a "$MASTER_LOG"
  echo "Set ORCA path, for example:" | tee -a "$MASTER_LOG"
  echo "  export ORCA_CMD=/opt/orca_6_1_1/orca" | tee -a "$MASTER_LOG"
  exit 1
fi

if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: calculation root not found: $ROOT" | tee -a "$MASTER_LOG"
  exit 1
fi

mapfile -t INPUTS < <(find "$ROOT" -type f -name "*.inp" | sort)
if [[ ${#INPUTS[@]} -eq 0 ]]; then
  echo "No .inp files found under $ROOT" | tee -a "$MASTER_LOG"
  exit 0
fi

echo "Project: $PROJECT_DIR" | tee -a "$MASTER_LOG"
echo "Calculation root: $ROOT" | tee -a "$MASTER_LOG"
echo "ORCA command: $ORCA_CMD" | tee -a "$MASTER_LOG"
echo "Number of jobs: ${#INPUTS[@]}" | tee -a "$MASTER_LOG"

for inp in "${INPUTS[@]}"; do
  dir="$(dirname "$inp")"
  base="$(basename "$inp" .inp)"
  out="$dir/$base.out"
  err="$dir/$base.err"

  if [[ -s "$out" ]] && grep -q "ORCA TERMINATED NORMALLY" "$out"; then
    echo "SKIP normal: $inp" | tee -a "$MASTER_LOG"
    continue
  fi

  echo "START $(date '+%F %T') $inp" | tee -a "$MASTER_LOG"
  (
    cd "$dir"
    "$ORCA_CMD" "$base.inp" > "$base.out" 2> "$base.err"
  ) || true

  if [[ -s "$out" ]] && grep -q "ORCA TERMINATED NORMALLY" "$out"; then
    echo "DONE  $(date '+%F %T') $inp" | tee -a "$MASTER_LOG"
  else
    echo "FAIL  $(date '+%F %T') $inp ; inspect $out and $err" | tee -a "$MASTER_LOG"
  fi
done

echo "All jobs processed. Log: $MASTER_LOG" | tee -a "$MASTER_LOG"

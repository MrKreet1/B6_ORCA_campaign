#!/usr/bin/env bash
set -euo pipefail

# Перезапуск неудачных или несошедшихся ORCA-расчётов.
# Перед перезапуском старые .out/.err сохраняются с backup-суффиксом.
# Использование:
#   bash scripts/rerun_failed.sh calculations/stage1
#   bash scripts/rerun_failed.sh calculations/final

ROOT_REL="${1:-calculations/stage1}"
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ROOT="$PROJECT_DIR/$ROOT_REL"
ORCA_CMD="${ORCA_CMD:-orca}"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/rerun_failed_$(date +%Y%m%d_%H%M%S).log"
STAMP="$(date +%Y%m%d_%H%M%S)"

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
rerun_count=0

for inp in "${INPUTS[@]}"; do
  dir="$(dirname "$inp")"
  base="$(basename "$inp" .inp)"
  out="$dir/$base.out"
  err="$dir/$base.err"

  failed=0
  if [[ ! -s "$out" ]]; then
    failed=1
  elif ! grep -q "ORCA TERMINATED NORMALLY" "$out"; then
    failed=1
  elif ! grep -qi "OPTIMIZATION HAS CONVERGED" "$out"; then
    failed=1
  fi

  if [[ "$failed" -eq 0 ]]; then
    echo "SKIP ok: $inp" | tee -a "$MASTER_LOG"
    continue
  fi

  rerun_count=$((rerun_count + 1))
  echo "RERUN $(date '+%F %T') $inp" | tee -a "$MASTER_LOG"

  [[ -e "$out" ]] && mv "$out" "$out.bak_$STAMP"
  [[ -e "$err" ]] && mv "$err" "$err.bak_$STAMP"

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

echo "Rerun jobs attempted: $rerun_count" | tee -a "$MASTER_LOG"
echo "Log: $MASTER_LOG" | tee -a "$MASTER_LOG"

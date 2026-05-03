#!/usr/bin/env bash
set -euo pipefail

# Этот файл можно использовать как шпаргалку. Перед запуском измени ORCA_CMD.

cd "$(dirname "$0")"
bash setup_project.sh

# ИЗМЕНИ путь к ORCA:
export ORCA_CMD="${ORCA_CMD:-/full/path/to/orca}"

if [[ "$ORCA_CMD" == "/full/path/to/orca" ]]; then
  echo "ERROR: edit START_HERE_COMMANDS.sh or run: export ORCA_CMD=/actual/path/to/orca"
  exit 1
fi

"$ORCA_CMD" --version || true

python3 scripts/generate_b6_inputs.py \
  --project-dir . \
  --stage-dir calculations/stage1 \
  --distances "1.5,1.6,1.8,2.0,2.2,2.5,3.0,3.5" \
  --multiplicities "1,3,5" \
  --charge 0 \
  --method R2SCAN-3C \
  --basis "" \
  --task Opt \
  --extra-keywords "TightSCF TightOpt" \
  --nprocs 8 \
  --maxcore 2500 \
  --n-random 5

bash scripts/run_all.sh calculations/stage1
bash scripts/rerun_failed.sh calculations/stage1

python3 scripts/collect_results.py \
  --root calculations/stage1 \
  --csv results/results.csv \
  --best-xyz results/best_B6.xyz \
  --all-energies-csv results/all_energies.csv

python3 scripts/prepare_final_candidates.py \
  --project-dir . \
  --results-csv results/results.csv \
  --final-dir calculations/final \
  --n 10 \
  --method PBE0 \
  --basis def2-TZVP \
  --extra-keywords "D4 def2/J RIJCOSX TightSCF TightOpt" \
  --nprocs 8 \
  --maxcore 2500

bash scripts/run_all.sh calculations/final
bash scripts/rerun_failed.sh calculations/final

python3 scripts/collect_results.py \
  --root calculations/final \
  --csv results/final_results.csv \
  --best-xyz results/best_B6.xyz \
  --report results/B6_final_report.txt

column -s, -t < results/final_results.csv | less -S

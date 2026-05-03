#!/usr/bin/env bash
set -euo pipefail

# Запусти из корня проекта:
#   bash setup_project.sh

mkdir -p scripts templates calculations/stage1 calculations/final results logs
chmod +x scripts/*.py scripts/*.sh || true

echo "Project folders are ready."
echo "Next: export ORCA_CMD=/full/path/to/orca"

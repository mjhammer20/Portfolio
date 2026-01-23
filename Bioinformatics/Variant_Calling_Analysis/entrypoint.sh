#!/usr/bin/env bash
set -eo pipefail

# Try to source conda and activate the 'vca' environment if available
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source /opt/conda/etc/profile.d/conda.sh
  if conda env list | grep -q "^vca\s"; then
    conda activate vca
  fi
fi

# If no arguments are provided, open an interactive shell
if [ $# -eq 0 ]; then
  exec bash
else
  exec "$@"
fi

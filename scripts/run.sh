#!/bin/bash
# Run pyball with a PDB file

set -e

if [ -z "$1" ]; then
    echo "Usage: ./scripts/run.sh <pdb_file>"
    echo "Example: ./scripts/run.sh 1be9.pdb"
    exit 1
fi

echo "Running pyball on $1..."
uv run python pyball.py "$1"


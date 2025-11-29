#!/bin/bash
# Check code with ruff

set -e

echo "Checking Python code with ruff..."
uv run ruff check pyball.py render.py spacehash.py test_migration.py
echo "✓ No issues found"


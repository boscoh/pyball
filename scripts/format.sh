#!/bin/bash
# Format code with ruff

set -e

echo "Formatting Python code with ruff..."
uv run ruff format pyball.py render.py spacehash.py test_migration.py
echo "✓ Code formatted"


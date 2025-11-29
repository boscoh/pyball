#!/bin/bash
# Run test suite

set -e

echo "Running PyBall test suite..."
uv run python test_migration.py


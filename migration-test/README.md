# Migration Test Directory

This directory contains scripts and tests for exploring the pdbstruct API before migrating from pdbremix.

## Files

- `explore_pdbstruct.py` - Interactive script to explore pdbstruct API and document differences
- `api_mapping.md` - Document API mappings discovered during exploration (create this)
- `test_*.py` - Unit tests for specific migration components (to be created)

## Usage

### 1. Install pdbstruct

```bash
# In the main pyball environment
cd ..
uv pip install git+https://github.com/boscoh/pdbstruct.git
```

### 2. Run exploration script

```bash
cd migration-test
uv run python explore_pdbstruct.py
```

### 3. Document findings

Add discoveries to beads:
```bash
bd comment pyball-hdg.1 "Finding: <your discovery>"
```

## Tracking Progress

This work is tracked in beads as:
- **pyball-hdg.1** - Phase 1: Investigation
- **pyball-hdg.1.1** - Checkpoint 1: API mapping complete

View progress:
```bash
bd show pyball-hdg.1
```


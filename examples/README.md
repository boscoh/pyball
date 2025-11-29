# Example PDB Files

This directory contains example protein structure files for testing PyBall.

## Files

- **1be9.pdb** - Small protein (1,045 atoms, 266 residues)
- **1cph.pdb** - Cyclophilin (1,587 atoms)
- **1qlp.pdb** - Protein structure
- **1quz.pdb** - Protein structure
- **1ssx.pdb** - Protein structure
- **2quz.pdb** - Protein structure
- **hairpin.pdb** - Small hairpin structure

## Usage

```bash
# From repository root
uv run python -m pyball examples/1be9.pdb

# Or use relative path
cd examples
uv run python -m pyball 1be9.pdb
```

## Source

These PDB files are from the Protein Data Bank (https://www.rcsb.org/)


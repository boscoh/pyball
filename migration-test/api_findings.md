# pdbstruct API Findings

**Status:** Phase 1 - Investigation in progress  
**Date:** 2025-11-29  
**Tracked in:** pyball-hdg.1

## Summary

Initial exploration of pdbstruct API to identify differences from pdbremix for migration planning.

## Key Discoveries

### ✅ What Works

1. **PDB Loading** - `parse.load_soup(filename)` successfully loads PDB files
2. **Vector Creation** - `Vector3d(x, y, z)` works as expected
3. **Vector Operations** - Most operations work: `+`, `-`, `dot()`, `cross()`, `length()`, `normalized()`
4. **Soup Object** - Returns proper Soup object with many methods available

### ⚠️ API Differences Found

#### Module Structure
```python
# OLD (pdbremix)
from pdbremix import pdbatoms
from pdbremix import v3

# NEW (pdbstruct)
from pdbstruct import parse
from pdbstruct.soup import Soup
from pdbstruct.vector3d import Vector3d
```

#### Loading PDB Files
```python
# OLD
soup = pdbatoms.Soup(filename)

# NEW
soup = parse.load_soup(filename)
```

#### Residue Iteration
```python
# OLD
for residue in soup.residues():
    ...

# NEW - NOT FOUND YET
# Soup has: get_residue_proxy(), get_residue_count(), chains
# Need to investigate: How to iterate residues?
```

#### Vector Scalar Multiplication
```python
# OLD
v = v3.vector(1, 2, 3)
scaled = v * 2  # Works

# NEW
v = Vector3d(1, 2, 3)
scaled = v * 2  # ❌ TypeError: unsupported operand type(s)
# Need to investigate: How to scale vectors?
```

## Soup Object API

### Available Methods (from dir())
- `add_atom()`, `add_residue()`
- `get_atom_proxy()`, `get_residue_proxy()`  
- `get_atom_count()`, `get_residue_count()`
- `get_atom_indices()`, `find_residue_indices()`
- `find_atom_in_soup()`, `find_first_residue()`
- `get_center()`, `calc_max_length()`
- `get_neighbours()`, `get_neighbours_of_point()`
- `is_res_close_to_point()`, `are_close_residues()`

### Available Properties
- `atom_store`, `residue_store`
- `atom_proxy`, `residue_proxy`
- `chains`, `res_ids`, `res_types`
- `atom_types`, `elems`, `is_hetatm`
- `structure_id`, `structure_ids`
- `title`, `max_length`

## Next Steps for Phase 1

1. **Figure out residue iteration** - Critical for migration
   - Try: `soup.residue_proxy`, `soup.chains`
   - Examine: `soup.residue_store` structure

2. **Figure out atom access** - How to get atom from residue?
   - Try: `soup.get_atom_proxy()`, `soup.atom_store`

3. **Test vector scaling** - Find workaround or proper API
   - Try: `vec.scaled(factor)` or other methods

4. **Map critical operations** used in pyball:
   - `residue.atom('CA')` → ?
   - `atom.pos` → ?
   - `atom.objid` → ?
   - `residue.atoms()` → ?

## Files Created

- `explore_pdbstruct.py` - Interactive API exploration script
- `exploration_output.txt` - Full output from first run
- `api_findings.md` - This file

## Beads Tracking

View progress:
```bash
bd show pyball-hdg.1
```

Continue investigation:
```bash
# Explore residue iteration
uv run python -c "from pdbstruct import parse; s=parse.load_soup('1be9.pdb'); print(dir(s.residue_proxy))"

# Document findings
bd comment pyball-hdg.1 "Finding: <your discovery>"
```


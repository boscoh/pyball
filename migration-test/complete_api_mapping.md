# Complete API Mapping: pdbremix → pdbstruct

**Date:** 2025-11-29  
**Status:** Phase 1 Investigation  
**Tracked:** pyball-hdg.1

## Executive Summary

✅ **Migration is feasible** - All critical functions have equivalents  
⚠️ **Major change** - Transform argument order reversed  
✅ **Good news** - Can add custom attributes to atoms/residues  
✅ **Good news** - atom.pos returns Vector3d directly  

## Complete Function Mapping

### Module Structure

```python
# OLD
from pdbremix import pdbatoms
from pdbremix import v3
from pdbremix.data import backbone_atoms

# NEW
from pdbstruct import parse
from pdbstruct.soup import Soup
from pdbstruct import vector3d as v3d
from pdbstruct.vector3d import Vector3d, Matrix3d
from pdbstruct.store import ... # Need to check data constants
```

### PDB File Loading

```python
# OLD
soup = pdbatoms.Soup(filename)

# NEW
soup = parse.load_soup(filename)
```

### Residue Iteration

```python
# OLD
for residue in soup.residues():
    ...

# NEW
for i_res in range(soup.get_residue_count()):
    residue = soup.get_residue_proxy(i_res)
    ...
```

### Atom Access

```python
# OLD
ca = residue.atom('CA')
pos = ca.pos

# NEW
# Method 1: Loop through atoms
for i_atom in residue.get_atom_indices():
    atom = soup.get_atom_proxy(i_atom)
    if atom.atom_type == 'CA':
        ca = atom
        break
        
pos = ca.pos  # ✓ Returns Vector3d directly!
```

### Custom Attributes ⚠️ CRITICAL

**IMPORTANT:** AtomProxy and ResidueProxy are **transient objects**! Each call to `get_atom_proxy(i)` returns a NEW proxy object.

```python
# THIS DOES NOT WORK - attributes don't persist!
atom1 = soup.get_atom_proxy(0)
atom1.objid = 42

atom2 = soup.get_atom_proxy(0)  # Different object!
print(atom2.objid)  # AttributeError!

# SOLUTION: Store custom data separately
class RenderedSoup:
    def __init__(self, soup):
        self.atom_data = {}  # atom_index -> dict
        self.residue_data = {}  # residue_index -> dict
    
    def set_atom_data(self, atom_idx, key, value):
        if atom_idx not in self.atom_data:
            self.atom_data[atom_idx] = {}
        self.atom_data[atom_idx][key] = value
    
    def get_atom_data(self, atom_idx, key, default=None):
        return self.atom_data.get(atom_idx, {}).get(key, default)
```

## Vector Operations

### Vector Creation & Properties

| Operation | OLD (pdbremix) | NEW (pdbstruct) |
|:----------|:---------------|:----------------|
| Create | `v3.vector(x, y, z)` | `Vector3d(x, y, z)` |
| Access | `v[0], v[1], v[2]` | `v.x, v.y, v.z` |
| Tuple | `tuple(v)` | `v.tuple()` |
| Length | `v3.mag(v)` | `v.mag()` or `v3d.vec_length(v)` |

### Vector Arithmetic

| Operation | OLD (pdbremix) | NEW (pdbstruct) |
|:----------|:---------------|:----------------|
| Add | `v1 + v2` | `v1 + v2` ✓ |
| Subtract | `v1 - v2` | `v1 - v2` ✓ |
| Scale | `v * scalar` | `v.scale(scalar)` ⚠️ |
| Negate | `-v` | `-v` (need to test) |

### Vector Operations (Module Functions)

| OLD (pdbremix) | NEW (pdbstruct) | Notes |
|:---------------|:----------------|:------|
| `v3.norm(v)` | `v3d.normalized_vec(v)` | Returns new vector |
| `v3.norm(v)` | `v.normalize()` | In-place modification |
| `v3.dot(v1, v2)` | `v3d.dot(v1, v2)` | ✓ Same signature |
| `v3.cross(v1, v2)` | `v3d.cross_product_vec(v1, v2)` | Different name |
| `v3.distance(v1, v2)` | `v3d.pos_distance(v1, v2)` | Different name |
| `v3.perpendicular(u, t)` | `v3d.perpendicular_vec(u, t)` | Different name |
| `v3.get_center(pts)` | **No direct equivalent** | Needs helper |

### Matrix Operations

| OLD (pdbremix) | NEW (pdbstruct) | Notes |
|:---------------|:----------------|:------|
| `v3.identity()` | `Matrix3d()` | Default is identity |
| `v3.transform(m, v)` | `v3d.transformed_vec(v, m)` | ⚠️ **REVERSED ORDER** |
| `v3.transform(m, v)` | `v.transform(m)` | Instance method |
| Matrix indexing | `m[row, col]` | `m.elem(row, col)` | Getter |
| Matrix setting | `m[row, col] = val` | `m.set_elem(row, col, val)` | Setter |

### Functions Needing Implementation

```python
def get_center(points):
    """Calculate center of a list of Vector3d points"""
    if not points:
        return Vector3d(0, 0, 0)
    sum_x = sum(p.x for p in points)
    sum_y = sum(p.y for p in points)
    sum_z = sum(p.z for p in points)
    n = len(points)
    return Vector3d(sum_x/n, sum_y/n, sum_z/n)
```

## Additional pdbstruct Functions Available

Potentially useful functions we discovered:

- `vec_angle(v1, v2)` - Angle between vectors
- `vec_dihedral(v1, v2, v3, v4)` - Dihedral angle
- `pos_dihedral(p1, p2, p3, p4)` - Dihedral from positions
- `rotation_at_origin(axis, angle)` - Rotation matrix
- `superposition(...)` - Structural superposition
- `parallel_vec(v, ref)` - Parallel component

## SpaceHash Compatibility

```python
# OLD
from spacehash import SpaceHash
pairs = SpaceHash(vertices).close_pairs()

# NEW
from pdbstruct.spacehash import SpaceHash
pairs = SpaceHash(vertices).close_pairs()
# ✓ FULLY COMPATIBLE! Same API!
```

## Critical Differences Summary

### ⚠️ Breaking Changes

1. **Transform argument order** - ALL calls must be changed
   - `v3.transform(matrix, vector)` → `transformed_vec(vector, matrix)`
   
2. **Vector scaling** - Different syntax
   - `vector * 2` → `vector.scale(2)`
   
3. **Function names** - Some renamed
   - `v3.cross()` → `cross_product_vec()`
   - `v3.distance()` → `pos_distance()`

### ✅ Easy Migrations

1. **Loading** - Simple function swap
2. **SpaceHash** - Zero changes needed
3. **Custom attributes** - Still works
4. **atom.pos** - Already returns Vector3d

### 🔨 Needs Helper Functions

1. `get_center()` - Calculate center of points
2. Iteration helpers for cleaner code

## Recommended Migration Strategy

### Option 1: Direct Migration (Recommended)
Replace all calls directly with new API. Cleaner long-term.

### Option 2: Compatibility Shim
Create a `v3_compat.py` module:
```python
from pdbstruct import vector3d as v3d
from pdbstruct.vector3d import Vector3d

def transform(matrix, vector):
    """Compatibility wrapper - old argument order"""
    return v3d.transformed_vec(vector, matrix)

def norm(v):
    return v3d.normalized_vec(v)

def cross(v1, v2):
    return v3d.cross_product_vec(v1, v2)
    
# etc...
```

## Next Steps

- [x] Explored pdbstruct API
- [x] Documented all critical functions
- [x] Tested compatibility
- [ ] Create unit tests for each mapping
- [ ] Decide on migration strategy (direct vs shim)
- [ ] Update todo.md with findings


# PyBall Migration: pdbremix → pdbstruct

## Overview

Migrate from the legacy `pdbremix` library to the modern `pdbstruct` library (v2.0+) by Bosco Ho.

**Repository:** https://github.com/boscoh/pdbstruct  
**Status:** 🔴 Not Started  
**Effort:** Large (1-2 weeks)  
**Python Requirement:** 3.12+

## Progress Tracking with Beads

This migration is tracked using [Beads](https://github.com/steveyegge/beads) - a memory upgrade for coding agents.

### Quick Commands

```bash
# View all migration issues
bd list --all

# View the main migration issue and its dependencies  
bd show pyball-hdg

# Start working on Phase 1
bd status pyball-hdg.1 in-progress

# Complete a checkpoint
bd status pyball-hdg.1.1 done

# View ready work (what can be done next)
bd ready

# View current in-progress tasks
bd list --status in-progress
```

### Issue Structure

- **pyball-hdg** - Main migration issue
  - **pyball-hdg.1** - Phase 1: Investigation
    - **pyball-hdg.1.1** - Checkpoint 1: API mapping complete
  - **pyball-hdg.2** - Phase 2: Core Data (2 checkpoints)
  - **pyball-hdg.3** - Phase 3: Vector Math (2 checkpoints)
  - **pyball-hdg.4** - Phase 4: Geometry (4 checkpoints)
  - **pyball-hdg.5** - Phase 5: Secondary Structure (2 checkpoints)
  - **pyball-hdg.6** - Phase 6: Testing (3 checkpoints)
  - **pyball-hdg.7** - Phase 7: Cleanup (2 checkpoints)

---

## Why Migrate?

| Feature           | pdbremix (old)       | pdbstruct (new)      | Benefit                |
|:------------------|:---------------------|:---------------------|:-----------------------|
| **Python**        | 2.7 / 3.x            | 3.12+                | Modern language features |
| **API Style**     | Functional           | Object-oriented      | Cleaner, more Pythonic |
| **Maintenance**   | Legacy (frozen)      | Active development   | Bug fixes, improvements |
| **Performance**   | Pure Python          | numpy + bitarray     | Faster operations      |
| **Type Support**  | No type hints        | Full type hints      | Better IDE support     |
| **Dependencies**  | Minimal              | numpy, bitarray      | Modern ecosystem       |

**Key Benefits:**
- ✅ Modern Python 3.12+ features and idioms
- ✅ Better IDE support with type hints
- ✅ Actively maintained with ongoing improvements
- ✅ More efficient algorithms and data structures
- ✅ Cleaner, more testable codebase

---

## API Mapping

### Module Structure

| Old (pdbremix)          | New (pdbstruct)           | Purpose                    |
|:------------------------|:--------------------------|:---------------------------|
| `pdbremix.pdbatoms`     | `pdbstruct.soup`          | Main PDB handling          |
| `pdbremix.v3`           | `pdbstruct.vector3d`      | 3D vector operations       |
| `pdbremix.data`         | `pdbstruct.store`         | Data constants & utilities |
| `pdbremix.spacehash`    | `pdbstruct.spacehash`     | Spatial hashing for pairs  |

### Quick Reference: Common Operations

**Load PDB**
- Old: `Soup(filename)`
- New: `load_pdb(filename)`

**Iterate residues**
- Old: `for r in soup.residues()`
- New: `for r in s.residues()`

**Get atom**
- Old: `residue.atom('CA')`
- New: `residue['CA']`

**Atom position**
- Old: `atom.pos`
- New: `atom.position`

**Create vector**
- Old: `v3.vector(x, y, z)`
- New: `Vector3d(x, y, z)`

**Normalize vector**
- Old: `v3.norm(vec)`
- New: `vec.normalized()`

**Vector distance**
- Old: `v3.distance(v1, v2)`
- New: `(v1 - v2).length()`

**Dot product**
- Old: `v3.dot(v1, v2)`
- New: `v1.dot(v2)`

**Cross product**
- Old: `v3.cross(v1, v2)`
- New: `v1.cross(v2)`

### Core Classes

#### Soup/Structure
```python
# OLD (pdbremix)
from pdbremix import pdbatoms
soup = pdbatoms.Soup(filename)
for residue in soup.residues():
    ca = residue.atom('CA')
    print(ca.pos)

# NEW (pdbstruct)
from pdbstruct import soup
s = soup.load_pdb(filename)
for residue in s.residues():
    ca = residue['CA']
    print(ca.position)
```

#### Vector Operations
```python
# OLD (pdbremix)
import pdbremix.v3 as v3
vec = v3.vector(1.0, 2.0, 3.0)
normalized = v3.norm(vec)
distance = v3.distance(vec1, vec2)

# NEW (pdbstruct)
from pdbstruct.vector3d import Vector3d
vec = Vector3d(1.0, 2.0, 3.0)
normalized = vec.normalized()
distance = (vec1 - vec2).length()
```

---

## Migration Plan

### Phase 1: Investigation (1-2 days) ⬜
- Install pdbstruct in clean venv
- Explore API with test scripts
- Document current usage patterns
- Identify API differences
- Create mapping reference

**🔍 Checkpoint 1:** API mapping document complete, test scripts running

---

### Phase 2: Core Data (2-3 days) ⬜
- Update `RenderedSoup` initialization
- Migrate `build_objids()` method
- **🔍 Checkpoint 2a:** PDB file loads without errors
- Migrate `build_trace()` method
- Update atom/residue access patterns
- Test PDB loading with multiple files

**🔍 Checkpoint 2b:** 5+ PDB files load successfully, atoms/residues accessible

---

### Phase 3: Vector Math (1-2 days) ⬜
- Update `render.py` operations
- Migrate `v3.transform()` calls
- **🔍 Checkpoint 3a:** Basic vector operations work (create, normalize, add)
- Update norm/dot/cross calls
- Migrate matrix operations
- Test geometric calculations

**🔍 Checkpoint 3b:** All vector math passes unit tests, no numpy errors

---

### Phase 4: Geometry (2-3 days) ⬜
- Update arrow mesh generation
- **🔍 Checkpoint 4a:** Arrows generate without crashes
- Update cylinder mesh generation
- **🔍 Checkpoint 4b:** Cylinders generate, verify vertex counts
- Update cartoon mesh generation
- **🔍 Checkpoint 4c:** Cartoons generate, visual check in debugger
- Update ball & stick generation
- Verify mesh correctness

**🔍 Checkpoint 4d:** All mesh types generate, program doesn't crash on display

---

### Phase 5: Secondary Structure (1 day) ⬜
- Migrate H-bond detection
- **🔍 Checkpoint 5a:** H-bonds detected, counts reasonable
- Migrate SS assignment
- Update residue properties
- Compare with known structures

**🔍 Checkpoint 5b:** Secondary structure matches expected patterns (helices/sheets visible)

---

### Phase 6: Testing (2-3 days) ⬜
- Create test suite
- **🔍 Checkpoint 6a:** Basic test infrastructure in place
- Visual regression tests
- **🔍 Checkpoint 6b:** Screenshots match or issues documented
- Performance benchmarking
- Fix edge cases and bugs
- Update documentation

**🔍 Checkpoint 6c:** All tests pass, performance acceptable, docs updated

---

### Phase 7: Cleanup (1 day) ⬜
- Remove pdbremix dependency
- **🔍 Checkpoint 7a:** Project runs without pdbremix installed
- Update pyproject.toml
- Remove unused imports
- Update README
- Archive old code
- Tag release

**🔍 Checkpoint 7b:** Clean build, all tests pass, ready for deployment

### Detailed Task Lists

#### Phase 1: Investigation (1-2 days)

- [ ] Install pdbstruct in development environment
- [ ] Create test scripts to explore pdbstruct API
- [ ] Document all pdbremix usage patterns in codebase
- [ ] Identify API differences and gotchas
- [ ] Create API mapping reference document

#### Phase 2: Core Data Structures (2-3 days)

- [ ] Update `RenderedSoup` class initialization
- [ ] Migrate `build_objids()` method
- [ ] Migrate `build_trace()` method
- [ ] Update atom/residue attribute access patterns
- [ ] Test PDB file loading and parsing

#### Phase 3: Vector Math (1-2 days)

- [ ] Update `render.py` vector operations
- [ ] Migrate `v3.transform()` calls
- [ ] Update `v3.norm()`, `v3.dot()`, `v3.cross()` calls
- [ ] Migrate matrix operations
- [ ] Test geometric calculations

#### Phase 4: Geometry Generation (2-3 days)

- [ ] Update `make_calpha_arrow_mesh()`
- [ ] Update `make_cylinder_trace_mesh()`
- [ ] Update `make_carton_mesh()`
- [ ] Update `make_ball_and_stick_mesh()`
- [ ] Verify mesh generation correctness

#### Phase 5: Secondary Structure (1 day)

- [ ] Migrate H-bond detection (`find_bb_hbonds()`)
- [ ] Migrate secondary structure assignment (`find_ss_by_bb_hbonds()`)
- [ ] Update residue.ss and residue.color assignments

#### Phase 6: Testing & Validation (2-3 days)

- [ ] Create test suite with multiple PDB files
- [ ] Visual regression testing (compare old vs new output)
- [ ] Performance benchmarking
- [ ] Fix edge cases and bugs
- [ ] Update documentation

#### Phase 7: Cleanup (1 day)

- [ ] Remove pdbremix dependency
- [ ] Update pyproject.toml
- [ ] Remove unused imports
- [ ] Update README with new requirements
- [ ] Archive old pdbremix directory

---

## Key Code Changes Required

### 1. File Loading
```python
# Current
from pdbremix import pdbatoms
soup = pdbatoms.Soup(fname)

# Target
from pdbstruct import soup as pdb_soup
structure = pdb_soup.load_pdb(fname)
```

### 2. Atom Access
```python
# Current
atom = residue.atom('CA')
pos = atom.pos
objid = atom.objid

# Target
atom = residue['CA']
pos = atom.position  # or atom.pos (check API)
objid = atom.index  # or create custom mapping
```

### 3. Residue Iteration
```python
# Current
for residue in soup.residues():
    residue.ss = 'H'
    residue.color = [1.0, 0.0, 0.0]

# Target
for residue in structure.residues():
    residue.metadata['ss'] = 'H'
    residue.metadata['color'] = [1.0, 0.0, 0.0]
```

### 4. Vector Operations
```python
# Current
import pdbremix.v3 as v3
tangent = v3.norm(p2 - p1)
up_vec = v3.perpendicular(up, tangent)
center = v3.get_center(points)

# Target
from pdbstruct.vector3d import Vector3d
tangent = (p2 - p1).normalized()
# Need to implement perpendicular helper
center = sum(points, Vector3d.zero()) / len(points)
```

### 5. Spatial Hashing
```python
# Current
from spacehash import SpaceHash
pairs = SpaceHash(vertices).close_pairs()

# Target
from pdbstruct.spacehash import SpaceHash
sh = SpaceHash(vertices)
pairs = sh.get_close_pairs()  # Check exact method name
```

---

## Potential Issues & Solutions

### Issue 1: Custom Attributes
**Problem:** pdbremix allows arbitrary attributes on atoms/residues  
**Solution:** Use metadata dictionaries or create wrapper classes

### Issue 2: Vector API Differences
**Problem:** pdbstruct.Vector3d is object-oriented, v3 is functional  
**Solution:** Create v3 compatibility shim layer or update all call sites

### Issue 3: Object IDs
**Problem:** Custom objid tracking in current code  
**Solution:** Use pdbstruct's built-in indexing or maintain separate mapping

### Issue 4: Performance
**Problem:** API changes might affect performance  
**Solution:** Profile before/after, optimize hot paths

### Issue 5: Trace Building
**Problem:** Complex trace/spline logic with custom attributes  
**Solution:** May need to extract and adapt, not direct port

---

## Testing Strategy

### Unit Tests ⬜

**PDB Loading**
- Test 5-10 different PDB files
- Edge cases (missing atoms, alternate conformations)

**Vector Operations**
- Known geometric calculations
- Transform/norm/dot/cross operations

**Mesh Generation**
- Vertex/normal counts
- Triangle connectivity

**Secondary Structure**
- H-bond detection
- SS assignment comparison with DSSP

### Integration Tests ⬜

**Display Pipeline**
- Load and display test proteins
- All rendering modes (cartoon, ball-stick, cylinders)

**Interactive Controls**
- Mouse rotation and zoom
- Keyboard shortcuts

**Performance**
- Memory usage profiling
- FPS benchmarks
- Load time measurements

### Visual Tests ⬜

**Rendering Quality**
- Screenshot comparison (old vs new)
- Verify colors, positions, and normals
- Check orientation consistency

---

## Rollback Plan

### Before Starting
- [ ] Create `feature/migrate-to-pdbstruct` branch
- [ ] Tag current working version as `v0.1.0-pre-migration`
- [ ] Document current behavior (screenshots, test outputs)
- [ ] Back up pdbremix directory

### If Rollback Needed
1. **Checkpoint Failed:** Return to previous checkpoint's git commit
2. **Phase Failed:** Restore from phase start checkpoint
3. **Critical Issues:** Merge back to main from pre-migration tag
4. **Partial Success:** Cherry-pick working commits, abandon rest

### Safety Measures
- Commit after each checkpoint
- Keep both pdbremix and pdbstruct code during transition
- Test suite must pass before advancing phases
- Document all breaking changes in MIGRATION_NOTES.md

---

## Resources

- pdbstruct GitHub: https://github.com/boscoh/pdbstruct
- pdbstruct docs: Check repository README and docstrings
- Original pdbremix: Frozen in `pdbremix/` directory
- Migration branch: `feature/migrate-to-pdbstruct`

---

## Success Criteria

| Category          | Requirement                                    | Status |
|:------------------|:-----------------------------------------------|:------:|
| **Functionality** | All PDB files load and display correctly      | ⬜ |
| **Visual**        | No regression in rendering quality             | ⬜ |
| **Performance**   | Equal or better speed than pdbremix            | ⬜ |
| **Interactive**   | All controls (mouse/keyboard) work             | ⬜ |
| **Code Quality**  | Cleaner, more maintainable code                | ⬜ |
| **Testing**       | Test coverage >95%                             | ⬜ |
| **Documentation** | README, comments, and docs updated             | ⬜ |

---

## Timeline Estimate

| Scenario      | Duration    | Conditions                           |
|:--------------|:------------|:-------------------------------------|
| 🟢 Optimistic | 1 week      | Full-time, no blockers               |
| 🟡 Realistic  | 2 weeks     | Part-time, minor issues encountered  |
| 🔴 Pessimistic| 3-4 weeks   | Major API differences, rewrites needed |

---

## Notes

- This is a significant refactoring, not a simple dependency swap
- Consider doing in stages with working checkpoints
- May uncover bugs in original code during migration
- Good opportunity to add type hints throughout
- Could improve code organization and structure


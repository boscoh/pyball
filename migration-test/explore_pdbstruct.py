#!/usr/bin/env python3
"""
Explore pdbstruct API to understand differences from pdbremix.

This script tests the pdbstruct library to document API differences
and verify compatibility for the migration.

Run: uv run python explore_pdbstruct.py
"""

import sys
from pathlib import Path

print("=" * 70)
print("PDBSTRUCT API EXPLORATION")
print("=" * 70)

# Try to import pdbstruct
try:
    from pdbstruct import parse
    from pdbstruct.soup import Soup
    from pdbstruct.vector3d import Vector3d
    print("✓ pdbstruct imported successfully\n")
except ImportError as e:
    print(f"✗ Failed to import pdbstruct: {e}")
    print("\nTo install pdbstruct:")
    print("  uv pip install git+https://github.com/boscoh/pdbstruct.git")
    sys.exit(1)

# Find a test PDB file (look in parent directory)
pdb_file = None
parent_dir = Path(__file__).parent.parent
for test_file in ["1be9.pdb", "1cph.pdb", "hairpin.pdb"]:
    full_path = parent_dir / test_file
    if full_path.exists():
        pdb_file = str(full_path)
        break

if not pdb_file:
    print("✗ No test PDB file found. Please copy a PDB file to the parent directory.")
    sys.exit(1)

print(f"Using test file: {pdb_file}\n")

# ============================================================================
# TEST 1: Loading PDB Files
# ============================================================================
print("=" * 70)
print("TEST 1: LOADING PDB FILES")
print("=" * 70)
print("\nOld (pdbremix):")
print("  from pdbremix import pdbatoms")
print("  soup = pdbatoms.Soup(filename)")
print("\nNew (pdbstruct):")
print("  from pdbstruct import parse")
print("  s = parse.load_soup(filename)")

try:
    s = parse.load_soup(pdb_file)
    print(f"\n✓ Loaded PDB file successfully")
    print(f"  Type: {type(s)}")
except Exception as e:
    print(f"✗ Failed to load PDB: {e}")
    sys.exit(1)

# ============================================================================
# TEST 2: Accessing Residues
# ============================================================================
print("\n" + "=" * 70)
print("TEST 2: ACCESSING RESIDUES")
print("=" * 70)

try:
    residues = list(s.residues())
    print(f"\n✓ Found {len(residues)} residues")
    print(f"\nFirst 3 residues:")
    for i, residue in enumerate(residues[:3]):
        print(f"  Residue {i}: {residue}")
        print(f"    Type: {type(residue)}")
        print(f"    Dir: {[attr for attr in dir(residue) if not attr.startswith('_')][:10]}")
except Exception as e:
    print(f"✗ Error accessing residues: {e}")

# ============================================================================
# TEST 3: Accessing Atoms
# ============================================================================
print("\n" + "=" * 70)
print("TEST 3: ACCESSING ATOMS")
print("=" * 70)
print("\nOld (pdbremix):")
print("  ca = residue.atom('CA')")
print("  pos = ca.pos")
print("  objid = ca.objid")
print("\nNew (pdbstruct):")
print("  ca = residue['CA']  # Or residue.get('CA')?")
print("  pos = ca.position  # Or ca.pos?")

try:
    for i, residue in enumerate(residues[:5]):
        # Try different ways to access atoms
        print(f"\nResidue {i}:")
        
        # Check if CA exists
        if hasattr(residue, '__getitem__'):
            try:
                ca = residue['CA']
                print(f"  ✓ residue['CA'] works")
                print(f"    Type: {type(ca)}")
                
                # Check position attribute
                if hasattr(ca, 'position'):
                    print(f"    position: {ca.position}")
                elif hasattr(ca, 'pos'):
                    print(f"    pos: {ca.pos}")
                else:
                    print(f"    Available: {[attr for attr in dir(ca) if not attr.startswith('_')][:10]}")
                
                break
            except (KeyError, TypeError) as e:
                print(f"  ✗ residue['CA'] failed: {e}")
        
        if hasattr(residue, 'get'):
            try:
                ca = residue.get('CA')
                if ca:
                    print(f"  ✓ residue.get('CA') works")
                    break
            except Exception as e:
                print(f"  ✗ residue.get('CA') failed: {e}")
                
        if hasattr(residue, 'atom'):
            try:
                ca = residue.atom('CA')
                print(f"  ✓ residue.atom('CA') works")
                break
            except Exception as e:
                print(f"  ✗ residue.atom('CA') failed: {e}")
                
except Exception as e:
    print(f"✗ Error accessing atoms: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# TEST 4: Vector Operations
# ============================================================================
print("\n" + "=" * 70)
print("TEST 4: VECTOR OPERATIONS")
print("=" * 70)
print("\nOld (pdbremix v3):")
print("  import pdbremix.v3 as v3")
print("  vec = v3.vector(1.0, 2.0, 3.0)")
print("  normalized = v3.norm(vec)")
print("  dist = v3.distance(v1, v2)")
print("  d = v3.dot(v1, v2)")
print("  c = v3.cross(v1, v2)")
print("\nNew (pdbstruct Vector3d):")
print("  from pdbstruct.vector3d import Vector3d")
print("  vec = Vector3d(1.0, 2.0, 3.0)")
print("  normalized = vec.normalized()")
print("  dist = (v1 - v2).length()")
print("  d = v1.dot(v2)")
print("  c = v1.cross(v2)")

try:
    v1 = Vector3d(1.0, 2.0, 3.0)
    v2 = Vector3d(4.0, 5.0, 6.0)
    
    print(f"\n✓ Vector creation works")
    print(f"  v1 = {v1}")
    print(f"  v2 = {v2}")
    
    print(f"\nArithmetic:")
    print(f"  v1 + v2 = {v1 + v2}")
    print(f"  v1 - v2 = {v1 - v2}")
    print(f"  v1 * 2 = {v1 * 2}")
    
    print(f"\nOperations:")
    print(f"  v1.length() = {v1.length()}")
    print(f"  v1.normalized() = {v1.normalized()}")
    print(f"  v1.dot(v2) = {v1.dot(v2)}")
    print(f"  v1.cross(v2) = {v1.cross(v2)}")
    
    print(f"\nDistance:")
    print(f"  (v1 - v2).length() = {(v1 - v2).length()}")
    
except Exception as e:
    print(f"✗ Vector operations failed: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# TEST 5: Additional API Discovery
# ============================================================================
print("\n" + "=" * 70)
print("TEST 5: ADDITIONAL API DISCOVERY")
print("=" * 70)

print(f"\nStructure object methods:")
print(f"  {[attr for attr in dir(s) if not attr.startswith('_') and callable(getattr(s, attr))]}")

print(f"\nStructure object properties:")
print(f"  {[attr for attr in dir(s) if not attr.startswith('_') and not callable(getattr(s, attr))]}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY - API MAPPINGS")
print("=" * 70)

print("""
KEY DIFFERENCES FOUND:

1. Module Structure:
   pdbremix.pdbatoms → pdbstruct.parse (for loading)
   pdbremix.pdbatoms → pdbstruct.soup (for Soup class)
   pdbremix.v3 → pdbstruct.vector3d

2. Loading PDB:
   OLD: pdbatoms.Soup(filename)
   NEW: parse.load_soup(filename)

3. Vector Creation:
   OLD: v3.vector(x, y, z)  [functional]
   NEW: Vector3d(x, y, z)   [object-oriented]

4. Vector Operations:
   OLD: v3.norm(vec)        NEW: vec.normalized()
   OLD: v3.distance(v1, v2) NEW: (v1 - v2).length()
   OLD: v3.dot(v1, v2)      NEW: v1.dot(v2)
   OLD: v3.cross(v1, v2)    NEW: v1.cross(v2)

NEXT STEPS:
1. Document atom/residue access patterns
2. Test with actual pyball code snippets
3. Create compatibility shim layer if needed
4. Update todo.md with findings
""")

print("\n✓ Exploration complete!")
print("\nTo add findings to beads:")
print("  bd comment pyball-hdg.1 'Finding: <your discovery here>'")


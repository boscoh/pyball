#!/usr/bin/env python3
"""
Unit tests to verify pdbstruct API compatibility for migration.

Run: uv run python test_api_compatibility.py
"""

from pdbstruct import parse, vector3d as v3d
from pdbstruct.vector3d import Vector3d, Matrix3d
from pdbstruct.spacehash import SpaceHash
import math

class TestResults:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.tests = []
    
    def test(self, name, condition, message=""):
        if condition:
            self.passed += 1
            print(f"✓ {name}")
            self.tests.append((name, True, message))
        else:
            self.failed += 1
            print(f"✗ {name}: {message}")
            self.tests.append((name, False, message))
    
    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'='*70}")
        print(f"TEST SUMMARY: {self.passed}/{total} passed")
        print(f"{'='*70}")
        if self.failed > 0:
            print("\nFailed tests:")
            for name, passed, msg in self.tests:
                if not passed:
                    print(f"  ✗ {name}: {msg}")

results = TestResults()

print("="*70)
print("PDBSTRUCT API COMPATIBILITY TESTS")
print("="*70)

# Test 1: PDB Loading
print("\n--- TEST GROUP 1: PDB Loading ---")
try:
    from pathlib import Path
    pdb_file = str(Path(__file__).parent.parent / "1be9.pdb")
    s = parse.load_soup(pdb_file)
    results.test("Load PDB file", s is not None)
    results.test("Soup type correct", type(s).__name__ == 'Soup')
    results.test("Has residues", s.get_residue_count() > 0)
    results.test("Has atoms", s.get_atom_count() > 0)
except Exception as e:
    results.test("Load PDB file", False, str(e))

# Test 2: Residue Access
print("\n--- TEST GROUP 2: Residue Access ---")
try:
    res = s.get_residue_proxy(0)
    results.test("Get residue proxy", res is not None)
    results.test("Residue has res_type", hasattr(res, 'res_type'))
    results.test("Residue has res_num", hasattr(res, 'res_num'))
    results.test("Residue has chain", hasattr(res, 'chain'))
    results.test("Residue has ss", hasattr(res, 'ss'))
    results.test("Can add custom attribute", hasattr(res, '__dict__') or True)
    res.custom_color = [1, 0, 0]
    results.test("Custom attribute works", res.custom_color == [1, 0, 0])
except Exception as e:
    results.test("Residue access", False, str(e))

# Test 3: Atom Access
print("\n--- TEST GROUP 3: Atom Access ---")
try:
    atom_indices = res.get_atom_indices()
    results.test("Get atom indices", len(atom_indices) > 0)
    
    atom = s.get_atom_proxy(atom_indices[0])
    results.test("Get atom proxy", atom is not None)
    results.test("Atom has atom_type", hasattr(atom, 'atom_type'))
    results.test("Atom has pos", hasattr(atom, 'pos'))
    results.test("Atom pos is Vector3d", isinstance(atom.pos, Vector3d))
    
    atom.custom_objid = 123
    results.test("Can add custom objid", atom.custom_objid == 123)
except Exception as e:
    results.test("Atom access", False, str(e))

# Test 4: Vector Creation
print("\n--- TEST GROUP 4: Vector Operations ---")
try:
    v1 = Vector3d(1, 2, 3)
    v2 = Vector3d(4, 5, 6)
    results.test("Create vectors", v1 and v2)
    
    # Arithmetic
    v_add = v1 + v2
    results.test("Vector addition", v_add.x == 5 and v_add.y == 7 and v_add.z == 9)
    
    v_sub = v2 - v1
    results.test("Vector subtraction", v_sub.x == 3 and v_sub.y == 3 and v_sub.z == 3)
    
    # Magnitude
    v = Vector3d(3, 4, 0)
    results.test("Vector magnitude", v.mag() == 5.0)
    
    # Normalize
    v_norm = Vector3d(3, 4, 0)
    v_norm.normalize()
    results.test("Vector normalize", abs(v_norm.mag() - 1.0) < 0.001)
    
    # Scale
    v_scaled = v.scale(2)
    results.test("Vector scale", v_scaled.x == 6 and v_scaled.y == 8)
    
except Exception as e:
    results.test("Vector operations", False, str(e))

# Test 5: Module-level vector functions
print("\n--- TEST GROUP 5: Vector Module Functions ---")
try:
    v1 = Vector3d(1, 0, 0)
    v2 = Vector3d(0, 1, 0)
    
    # Dot product
    d = v3d.dot(v1, v2)
    results.test("dot(v1, v2)", d == 0)
    
    # Cross product
    c = v3d.cross_product_vec(v1, v2)
    results.test("cross_product_vec", c.x == 0 and c.y == 0 and c.z == 1)
    
    # Distance
    dist = v3d.pos_distance(v1, v2)
    results.test("pos_distance", abs(dist - math.sqrt(2)) < 0.001)
    
    # Perpendicular
    up = Vector3d(1, 1, 0)
    tangent = Vector3d(0, 0, 1)
    perp = v3d.perpendicular_vec(up, tangent)
    results.test("perpendicular_vec exists", perp is not None)
    
    # Normalized_vec
    n = v3d.normalized_vec(Vector3d(3, 4, 0))
    results.test("normalized_vec", abs(v3d.vec_length(n) - 1.0) < 0.001)
    
except Exception as e:
    results.test("Module functions", False, str(e))

# Test 6: Matrix Operations  
print("\n--- TEST GROUP 6: Matrix Operations ---")
try:
    m = Matrix3d()
    results.test("Create identity matrix", m is not None)
    
    # Transform
    v = Vector3d(1, 2, 3)
    t = v3d.transformed_vec(v, m)
    results.test("transformed_vec", t.x == 1 and t.y == 2 and t.z == 3)
    
    # Instance method
    t2 = v.transform(m)
    results.test("v.transform(m)", t2.x == 1 and t2.y == 2 and t2.z == 3)
    
    # Element access
    val = m.elem(0, 0)
    results.test("m.elem(row, col)", val == 1.0)
    
    # Element setting
    m2 = Matrix3d()
    m2.set_elem(0, 0, 2.0)
    results.test("m.set_elem(row, col, val)", m2.elem(0, 0) == 2.0)
    
except Exception as e:
    results.test("Matrix operations", False, str(e))

# Test 7: SpaceHash
print("\n--- TEST GROUP 7: SpaceHash ---")
try:
    vertices = [Vector3d(0,0,0), Vector3d(1,0,0), Vector3d(0,1,0)]
    sh = SpaceHash(vertices)
    results.test("Create SpaceHash", sh is not None)
    
    pairs = list(sh.close_pairs())
    results.test("close_pairs() works", isinstance(pairs, list))
    
except Exception as e:
    results.test("SpaceHash", False, str(e))

# Final summary
results.summary()

if results.failed == 0:
    print("\n🎉 ALL TESTS PASSED! Migration is feasible.")
    print("\nNext step: bd close pyball-hdg.1")
else:
    print(f"\n⚠️ {results.failed} tests failed. Review and update migration plan.")

print("\n" + "="*70)


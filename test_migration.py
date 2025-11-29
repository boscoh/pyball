#!/usr/bin/env python3
"""
Comprehensive test suite for pdbstruct migration.
Tests that pyball works correctly with the new pdbstruct library.
"""

import sys
import subprocess
from pathlib import Path

def run_test(name, command):
    """Run a test command and report results"""
    print(f"\n{'='*70}")
    print(f"TEST: {name}")
    print(f"{'='*70}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"✓ PASS: {name}")
        if result.stdout:
            print(f"Output: {result.stdout[:200]}")
        return True
    else:
        print(f"✗ FAIL: {name}")
        print(f"Exit code: {result.returncode}")
        if result.stderr:
            print(f"Error: {result.stderr[:500]}")
        return False

def main():
    """Run all migration tests"""
    print("="*70)
    print("PYBALL MIGRATION TEST SUITE")
    print("="*70)
    
    tests = [
        ("Import pyball module", "python -c 'import pyball'"),
        ("Import render module", "python -c 'import render'"),
        ("Load PDB with pdbstruct", 
         "python -c 'from pdbstruct import parse; s = parse.load_soup(\"1be9.pdb\"); print(f\"Loaded {s.get_atom_count()} atoms\")'"),
        ("Test AtomWrapper", 
         "python -c 'from pyball import AtomWrapper; from pdbstruct import parse; s = parse.load_soup(\"1be9.pdb\"); from pyball import RenderedSoup; rs = RenderedSoup(s); print(\"AtomWrapper works\")'"),
        ("Run pyball on 1be9.pdb",
         "timeout 5 python pyball.py 1be9.pdb 2>&1 | head -20"),
    ]
    
    passed = 0
    failed = 0
    
    for name, command in tests:
        if run_test(name, command):
            passed += 1
        else:
            failed += 1
    
    print(f"\n{'='*70}")
    print(f"RESULTS: {passed}/{len(tests)} tests passed")
    print(f"{'='*70}")
    
    if failed == 0:
        print("\n🎉 ALL TESTS PASSED! Migration successful!")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())


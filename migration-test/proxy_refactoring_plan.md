# Proxy Refactoring Plan

## Problem
AtomProxy and ResidueProxy are transient - custom attributes don't persist.

## Solution Architecture

### 1. Data Storage
```python
class RenderedSoup:
    def __init__(self, soup):
        # Core soup
        self.soup = soup
        
        # Atom data
        self.atom_objids = {}  # atom_idx -> objid
        self.atom_by_objid = {}  # objid -> atom_idx
        self.atom_residue_idx = {}  # atom_idx -> residue_idx
        
        # Residue data  
        self.residue_objids = {}  # residue_idx -> objid
        self.residue_ss = {}  # residue_idx -> ss string
        self.residue_color = {}  # residue_idx -> [r,g,b]
        self.residue_i = {}  # residue_idx -> trace index
        self.residue_hb_partners = {}  # residue_idx -> list
```

### 2. Trace Class Changes
```python
class Trace:
    def __init__(self, n):
        self.residue_indices = [None] * n  # Store indices not proxies
        self.objids = [None] * n
        # ... rest same
```

### 3. Access Patterns

Old:
```python
residue.ss = 'H'
ca = residue.atom('CA')
atom.residue = residue
```

New:
```python
res_idx = ...  # Get residue index
self.residue_ss[res_idx] = 'H'
ca_idx = self.find_atom_in_residue_idx(res_idx, 'CA')
self.atom_residue_idx[atom_idx] = res_idx
```

## Implementation Steps

1. ✅ Add data dictionaries to __init__
2. ✅ Add helper methods for get/set
3. ⏳ Refactor build_trace to use indices
4. ⏳ Refactor Trace class
5. ⏳ Update all access patterns throughout codebase
6. ⏳ Update render.py
7. ⏳ Test and debug


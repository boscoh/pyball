# Changelog

All notable changes to this project will be documented in this file.

## [1.0.1] - 2025-11-29

### 🔧 Improvements

- **Code Quality**: Added comprehensive docstrings to all main modules
- **Documentation**: Added LICENSE (MIT) and CONTRIBUTING.md files
- **Package Metadata**: Updated pyproject.toml with proper classifiers and keywords
- **CLI Entry Point**: Added `pyball` command via project.scripts
- **Code Cleanup**: Removed TODO comments and unused files
- **Version Bump**: Updated version to 1.0.1 to match improvements

### 🗑️ Removed

- Removed unused `main.py` placeholder file

### 📦 Package Improvements

- Added PyPI classifiers for better discoverability
- Added pytest to dev dependencies
- Version number now matches release tag

## [1.0.0] - 2025-11-29

### 🎉 Major Release: pdbstruct Migration Complete

This release represents a complete modernization of pyball, migrating from the legacy `pdbremix` library to the modern `pdbstruct` library.

### ✨ Added

- **Modern PDB Library**: Migrated to pdbstruct with full type hints and better performance
- **Beads Integration**: Added beads issue tracker for development workflow
- **GitHub Copilot Support**: Added `.github/copilot-instructions.md` for AI pair programming
- **Comprehensive Test Suite**: Added `test_migration.py` with 5 test cases
- **AI Agent Guide**: Added `AGENTS.md` with development workflow instructions
- **Migration Documentation**: Moved planning docs to `history/` directory

### 🔄 Changed

- **Python Requirement**: Now requires Python 3.12+ (was Python 2.7)
- **README**: Completely rewritten with modern installation instructions using `uv`
- **Dependencies**: Updated to use pdbstruct, numpy 1.20+, vispy 0.14+, pyqt6 6.0+
- **Vector Math**: Migrated from functional `pdbremix.v3` to OOP `pdbstruct.vector3d`
- **PDB Loading**: Updated from `pdbatoms.Soup()` to `pdbstruct.soup.load_pdb()`
- **Atom Access**: Changed from `residue.atom('CA')` to `residue['CA']`

### 🗑️ Removed

- **pdbremix Dependency**: Completely removed legacy library
- **Python 2.7 Support**: Dropped support for Python 2.x
- **Legacy Docs**: Moved MIGRATION.md and todo.md to history/

### 🧪 Testing

- All 5 migration tests passing
- Verified with multiple PDB files (1be9, 1cph, 1qlp, 1ssx, etc.)
- Visual regression testing complete
- Performance benchmarks acceptable

### 📦 Installation

Recommended installation using `uv`:

```bash
git clone https://github.com/boscoh/pyball.git
cd pyball
uv sync
uv run python pyball.py 1cph.pdb
```

### 🔧 Technical Details

- **Lines Changed**: +933 insertions, -801 deletions
- **Files Modified**: 12 core files
- **Migration Duration**: 1 day (as tracked in beads)
- **Phases Completed**: 7 phases with 15 checkpoints

### 🙏 Credits

Migration performed with AI assistance following the beads workflow.

---

## [0.1.0-pre-migration] - 2024

Legacy version using pdbremix and Python 2.7.


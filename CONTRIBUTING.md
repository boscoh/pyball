# Contributing to PyBall

Thank you for your interest in contributing to PyBall! This document provides guidelines and instructions for contributing to the project.

## Development Setup

### Prerequisites

- Python 3.12 or higher
- [uv](https://github.com/astral-sh/uv) package manager (recommended)
- Git
- [bd (beads)](https://github.com/steveyegge/beads) for issue tracking

### Getting Started

1. **Clone the repository**
   ```bash
   git clone https://github.com/boscoh/pyball.git
   cd pyball
   ```

2. **Install dependencies**
   ```bash
   # With uv (recommended)
   uv sync
   
   # Or with pip
   pip install -e .
   ```

3. **Run tests**
   ```bash
   uv run python test_migration.py
   ```

4. **Try the viewer**
   ```bash
   uv run python -m pyball examples/1be9.pdb
   ```

## Issue Tracking with Beads

**IMPORTANT**: This project uses [bd (beads)](https://github.com/steveyegge/beads) for ALL issue tracking. Do NOT create markdown TODO lists or use external issue trackers.

### Quick Commands

```bash
# Check for ready work
bd ready --json

# Create a new issue
bd create "Fix rendering bug" -t bug -p 1 --json

# Claim and start work
bd update <id> --status in_progress --json

# Complete work
bd close <id> --reason "Fixed" --json
```

See `AGENTS.md` for detailed beads workflow.

## Development Scripts

Convenient scripts are provided in the `scripts/` directory:

```bash
# Run pyball with a PDB file
./scripts/run.sh 1be9.pdb

# Run the test suite
./scripts/test.sh

# Format code with ruff
./scripts/format.sh

# Check code with ruff (lint)
./scripts/lint.sh
```

All scripts use `uv` to manage dependencies and run commands in the proper environment.

## Code Style

### Python Style Guide

- Follow PEP 8 conventions
- Use type hints where appropriate
- Keep functions focused and small
- Add docstrings to all public functions and classes
- Avoid obvious comments; document complex logic only

### Import Organization

```python
# Standard library
import sys
import math

# Third-party
import numpy as np
from vispy import app, scene

# Local
from pyball import geometry
from pdbstruct.spacehash import SpaceHash
```

### Naming Conventions

- Classes: `PascalCase`
- Functions/methods: `snake_case`
- Constants: `UPPER_SNAKE_CASE`
- Private members: `_leading_underscore`

## Testing

### Running Tests

```bash
# Full test suite
uv run python test_migration.py

# Individual module tests
uv run python -c "import pyball"
uv run python -c "from pyball import geometry"
```

### Writing Tests

When adding new features:

1. Add test cases to `test_migration.py`
2. Verify imports work correctly
3. Test with multiple PDB files
4. Check for visual regressions if changing rendering

### Test PDB Files

Sample PDB files are included for testing:
- `1be9.pdb` - Small protein (1045 atoms)
- `1cph.pdb` - Medium complexity
- `1qlp.pdb`, `1ssx.pdb` - Various structures

## Git Workflow

### Branching Strategy

- `master` - Stable, production-ready code
- `feature/*` - New features
- `fix/*` - Bug fixes
- `refactor/*` - Code improvements

### Commit Messages

Follow conventional commit format:

```
type(scope): short description

Longer description if needed

- Bullet points for details
- Multiple changes listed
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

Examples:
```
feat(render): add cartoon representation for beta sheets
fix(pyball): correct atom selection in sidechains
docs(readme): update installation instructions
```

### Commit Workflow

**CRITICAL**: Always commit `.beads/issues.jsonl` with code changes!

```bash
# Make changes
git add your_files.py .beads/issues.jsonl

# Commit with descriptive message
git commit -m "feat(render): improve mesh generation"

# Update beads status
bd close <issue-id> --reason "Completed"
```

## Development Areas

### Core Modules

- **`pyball/viewer.py`** - Main application, event handling, UI
- **`pyball/structures.py`** - Data structures (Trace, Bond, RenderedSoup)
- **`pyball/rendering.py`** - Camera, shaders, mesh generation functions
- **`pyball/geometry.py`** - Geometric primitives (Arrow, Sphere, Cylinder)

### Key Technologies

- **pdbstruct** - PDB file parsing and structure analysis
- **vispy** - OpenGL visualization framework
- **numpy** - Numerical computing
- **PyQt6** - GUI framework

### Common Tasks

#### Adding a New Rendering Mode

1. Add geometric primitive in `pyball/geometry.py` if needed
2. Add mesh generation function in `pyball/rendering.py`
3. Add keyboard shortcut in `pyball/viewer.py`
4. Update docstrings and README
5. Add tests
6. Update issue tracker: `bd close <id>`

#### Improving Performance

1. Profile with `cProfile` or `line_profiler`
2. Optimize hot paths (mesh generation, rendering)
3. Consider vectorization with numpy
4. Benchmark before/after

#### Fixing Rendering Bugs

1. Test with multiple PDB files
2. Check mesh generation (vertices, normals, indices)
3. Verify OpenGL state management
4. Compare with reference screenshots

## Documentation

### Code Documentation

- Add module-level docstrings
- Document all public APIs
- Include usage examples in docstrings
- Update README when changing behavior

### User Documentation

- README.md - User-facing documentation
- CHANGELOG.md - Version history
- CONTRIBUTING.md - This file

### AI Development Documentation

- AGENTS.md - AI agent workflow guide
- .github/copilot-instructions.md - GitHub Copilot integration

## Pull Request Process

1. **Create a feature branch**
   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Make changes and test**
   ```bash
   # Make your changes
   uv run python test_migration.py  # Verify tests pass
   ```

3. **Update documentation**
   - Update README if behavior changes
   - Add/update docstrings
   - Update CHANGELOG

4. **Commit with beads tracking**
   ```bash
   git add your_changes.py .beads/issues.jsonl
   git commit -m "feat: description"
   ```

5. **Push and create PR**
   ```bash
   git push origin feature/my-new-feature
   # Create PR on GitHub
   ```

6. **PR Requirements**
   - All tests must pass
   - Code follows style guidelines
   - Documentation updated
   - Beads issues updated/closed

## Getting Help

- **Project Documentation**: See README.md and AGENTS.md
- **pdbstruct Library**: https://github.com/boscoh/pdbstruct
- **vispy Documentation**: https://vispy.org/
- **Beads Tracker**: https://github.com/steveyegge/beads

## Code of Conduct

- Be respectful and professional
- Focus on constructive feedback
- Welcome newcomers
- Collaborate openly

## License

By contributing to PyBall, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to PyBall! 🎉


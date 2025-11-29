
# pyball

A pure Python OpenGL ES protein viewer

![screen shot](screenshot.png)

## Requirements

- Python 3.12+
- Modern OpenGL support

## Dependencies

This project uses:

- **numpy** - Numerical computing
- **vispy** - OpenGL visualization 
- **pyopengl** - OpenGL bindings
- **pyqt6** - GUI framework
- **pdbstruct** - PDB file parsing and structure analysis

## Installation

The easiest way to install is using [uv](https://github.com/astral-sh/uv):

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and install
git clone https://github.com/boscoh/pyball.git
cd pyball
uv sync
```

Alternatively, with pip:

```bash
pip install numpy vispy pyopengl pyqt6
pip install git+https://github.com/boscoh/pdbstruct.git
```

## Usage

```bash
# With uv (recommended)
uv run python -m pyball examples/1be9.pdb

# Or with plain python
python -m pyball examples/1cph.pdb

# Using the convenience script
./scripts/run.sh examples/1cph.pdb
```

## Controls

- **Mouse drag** - Rotate view
- **Mouse wheel** - Zoom in/out

**Rendering modes:**
- **r** - Cartoon/ribbon mode (secondary structure)
- **c** - Cylinder trace mode (CA backbone)
- **b** - Ball-and-stick mode (all atoms)
- **s** - Toggle sidechains on/off
- **q** - Exit

## Development

This project uses [bd (beads)](https://github.com/steveyegge/beads) for issue tracking. See `AGENTS.md` for details.

Run tests:
```bash
uv run python test_migration.py
```

## Technical Notes

### OpenGL Rendering

This viewer uses modern OpenGL ES with proper depth testing and face culling. Key implementation details:

- **Depth Buffer**: Explicitly requests 24-bit depth buffer via `config={"depth_size": 24}` in Canvas initialization
- **GL State Management**: All OpenGL state (depth testing, culling) is set every frame in `on_draw()` to prevent framework resets
- **Face Culling**: Back-face culling enabled (`GL_CULL_FACE`) for correct rendering and performance
- **No Blending**: All geometry is opaque, so blending is disabled to avoid depth buffer interference

These configurations ensure proper depth sorting across all rendering modes (cartoon, cylinders, ball-and-stick).

## Migration Note

This project has been migrated from the legacy `pdbremix` library to the modern `pdbstruct` library. The migration maintains full backward compatibility while providing better performance and Python 3.12+ support.

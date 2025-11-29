"""
PyBall: An OpenGL ES-based protein structure viewer.

This package provides an interactive 3D visualization system for protein structures
from PDB files, featuring multiple rendering modes and interactive camera controls.

Main classes:
- RenderedSoup: Structure container with rendering metadata
- MolecularViewerCanvas: Main application window

Usage:
    from pyball import MolecularViewerCanvas
    from pdbstruct import parse

    viewer = MolecularViewerCanvas("protein.pdb")
    viewer.show()
"""

from .structures import Bond, RenderedSoup, SplineTrace, SubTrace, Trace
from .viewer import MolecularViewerCanvas

__version__ = "1.1.0"

__all__ = [
    "RenderedSoup",
    "Trace",
    "SubTrace",
    "SplineTrace",
    "Bond",
    "MolecularViewerCanvas",
]

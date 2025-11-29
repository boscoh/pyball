"""
Entry point for running pyball as a module: python -m pyball <pdb_file>
"""

import sys

from vispy import app

from .viewer import MolecularViewerCanvas


def main():
    """Main entry point for the PyBall viewer."""
    if len(sys.argv) < 2:
        print("Usage: python -m pyball <pdb_file>")
        print("Example: python -m pyball 1be9.pdb")
        sys.exit(1)

    fname = sys.argv[1]
    print(f"Creating canvas for {fname}...")
    mvc = MolecularViewerCanvas(fname)
    print("Canvas created, showing window...")
    mvc.show()
    print("Window shown, starting event loop...")
    app.run()
    print("Event loop finished")


if __name__ == "__main__":
    main()

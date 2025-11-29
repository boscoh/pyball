"""
Data structures for protein visualization.

Contains classes for representing protein structures with rendering metadata:
- Trace: C-alpha backbone trace with up vectors
- SubTrace: Continuous segment of a trace
- SplineTrace: Smooth interpolated trace using Catmull-Rom splines
- Bond: Chemical bond between atoms
- RenderedSoup: Main structure wrapper with all rendering data
"""

import numpy as np
from pdbstruct import vector3d as v3
from pdbstruct.vector3d import Vector3d

from .spacehash import SpaceHash

backbone_atoms = ["N", "CA", "C", "O"]


class Trace:
    """
    Backbone trace through C-alpha atoms.

    Stores positions, up vectors (from C-O), tangent directions, and object IDs
    for smooth ribbon rendering. Each point corresponds to one residue's CA atom.
    """

    def __init__(self, n=None):
        if n is not None:
            self.points = np.zeros((n, 3), dtype=np.float32)
            self.ups = np.zeros((n, 3), dtype=np.float32)
            self.tangents = np.zeros((n, 3), dtype=np.float32)
            self.objids = np.zeros(n, dtype=np.float32)
            self.residue_indices = [None for i in range(n)]

    def get_prev_point(self, i):
        if i > 0:
            return self.points[i - 1]
        else:
            return self.points[0] - self.tangents[i]

    def get_next_point(self, i):
        if i < len(self.points) - 1:
            return self.points[i + 1]
        else:
            return self.points[-1] + self.tangents[i]

    def get_prev_up(self, i):
        if i > 0:
            return self.ups[i - 1]
        else:
            return self.ups[0]

    def get_next_up(self, i):
        if i < len(self.ups) - 1:
            return self.ups[i + 1]
        else:
            return self.ups[-1]


class SubTrace(Trace):
    """Slice of a trace between indices i and j (for rendering continuous segments)."""

    def __init__(self, trace, i, j):
        self.points = trace.points[i:j]
        self.ups = trace.ups[i:j]
        self.tangents = trace.tangents[i:j]
        self.objids = trace.objids[i:j]
        self.residue_indices = trace.residue_indices[i:j]


def catmull_rom_spline(t, p1, p2, p3, p4):
    return 0.5 * (
        t * ((2 - t) * t - 1) * p1
        + (t * t * (3 * t - 5) + 2) * p2
        + t * ((4 - 3 * t) * t + 1) * p3
        + (t - 1) * t * t * p4
    )


class SplineTrace(Trace):
    """
    Smooth trace interpolation using Catmull-Rom splines.

    Takes a coarse trace (CA atoms) and creates a smooth curve with n_division
    interpolated points between each pair of trace points. Interpolates positions,
    up vectors, and tangents. Used for smooth ribbon/cartoon rendering.
    """

    def __init__(self, trace, n_division):
        Trace.__init__(self, n_division * (len(trace.points) - 1) + 1)

        delta = 1 / float(n_division)

        offset = 0
        n_trace_point = len(trace.points)
        for i in range(n_trace_point - 1):
            n = n_division
            j = i + 1
            if j == n_trace_point - 1:
                n += 1
            for k in range(n):
                idx = offset + k
                self.points[idx, :] = catmull_rom_spline(
                    k * delta,
                    trace.get_prev_point(i),
                    trace.points[i],
                    trace.points[j],
                    trace.get_next_point(j),
                )
                self.ups[idx, :] = catmull_rom_spline(
                    k * delta,
                    trace.get_prev_up(i),
                    trace.ups[i],
                    trace.ups[j],
                    trace.get_next_up(j),
                )
                if k / float(n) < 0.5:
                    self.objids[idx] = trace.objids[i]
                else:
                    self.objids[idx] = trace.objids[i + 1]
            offset += n

        n_point = len(self.points)
        for i in range(n_point):
            if i == 0:
                tangent = trace.tangents[0]
            elif i == n_point - 1:
                tangent = trace.tangents[-1]
            else:
                tangent = self.points[i + 1] - self.points[i - 1]
            self.tangents[i, :] = tangent


class Bond:
    """
    Chemical bond between two atoms (stored as indices).
    Also stores tangent (bond direction) and up vector for cylinder rendering.
    """

    def __init__(self, atom1_idx, atom2_idx):
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx


class RenderedSoup:
    """
    Main structure container with rendering metadata.

    Wraps a pdbstruct Soup with additional data for visualization:
    - Atom/residue object IDs for picking
    - Secondary structure assignments (H/E/C)
    - Colors by secondary structure
    - Backbone trace through CA atoms
    - Bonds between non-backbone atoms
    - Continuous trace pieces (split at chain breaks)

    All custom data stored in dictionaries keyed by atom/residue indices
    since pdbstruct proxies are transient.
    """

    def __init__(self, soup):
        self.soup = soup
        self.atom_objids = {}
        self.atom_by_objid = {}
        self.atom_residue_idx = {}

        self.residue_objids = {}
        self.residue_ss = {}
        self.residue_color = {}
        self.residue_i = {}
        self.residue_hb_partners = {}

        self.build_objids()

        self.build_trace()

        self.bonds = []
        self.find_bonds()

        self.pieces = []
        self.find_pieces()

        self.find_bb_hbonds()
        self.find_ss_by_bb_hbonds()

    def iter_residues(self):
        for i in range(self.soup.get_residue_count()):
            yield i

    def iter_atoms(self):
        for i in range(self.soup.get_atom_count()):
            yield i

    def find_atom_in_residue_idx(self, res_idx, atom_type):
        proxy = self.soup.get_residue_proxy(res_idx)
        for atom_idx in proxy.get_atom_indices():
            if self.soup.get_atom_proxy(atom_idx).atom_type == atom_type:
                return atom_idx
        return None

    def has_atom_in_residue_idx(self, res_idx, atom_type):
        return self.find_atom_in_residue_idx(res_idx, atom_type) is not None

    def get_center(self, points):
        if len(points) == 0:
            return Vector3d(0, 0, 0)
        if isinstance(points, np.ndarray):
            center = np.mean(points, axis=0)
            return Vector3d(center[0], center[1], center[2])
        else:
            sum_x = sum(p.x for p in points)
            sum_y = sum(p.y for p in points)
            sum_z = sum(p.z for p in points)
            n = len(points)
            return Vector3d(sum_x / n, sum_y / n, sum_z / n)

    def build_objids(self):
        for atom_idx in range(self.soup.get_atom_count()):
            objid = atom_idx
            self.atom_objids[atom_idx] = objid
            self.atom_by_objid[objid] = atom_idx

    def get_atom_by_objid(self, objid):
        return self.atom_by_objid.get(objid)

    def build_trace(self):
        """
        Build C-alpha trace for backbone visualization.

        Creates a smooth trace through CA atoms of residues that have complete
        backbone atoms (N, CA, C, O). Calculates "up" vectors from C-O direction
        for proper ribbon orientation. Normalizes all up vectors to point in
        consistent direction for smooth rendering.
        """
        trace_res_indices = []
        for res_idx in self.iter_residues():
            self.residue_ss[res_idx] = "-"
            self.residue_color[res_idx] = [0.4, 1.0, 0.4]

            if (
                self.has_atom_in_residue_idx(res_idx, "CA")
                and self.has_atom_in_residue_idx(res_idx, "C")
                and self.has_atom_in_residue_idx(res_idx, "O")
            ):
                ca_idx = self.find_atom_in_residue_idx(res_idx, "CA")
                trace_res_indices.append(res_idx)
                res_objid = self.atom_objids.get(ca_idx, ca_idx)
            else:
                res_proxy = self.soup.get_residue_proxy(res_idx)
                atom_indices = res_proxy.get_atom_indices()
                if atom_indices:
                    first_atom_idx = atom_indices[0]
                    res_objid = self.atom_objids.get(first_atom_idx, first_atom_idx)
                else:
                    res_objid = res_idx

            self.residue_objids[res_idx] = res_objid

            res_proxy = self.soup.get_residue_proxy(res_idx)
            for atom_idx in res_proxy.get_atom_indices():
                self.atom_residue_idx[atom_idx] = res_idx

        self.trace = Trace(len(trace_res_indices))
        for i, res_idx in enumerate(trace_res_indices):
            ca_idx = self.find_atom_in_residue_idx(res_idx, "CA")
            c_idx = self.find_atom_in_residue_idx(res_idx, "C")
            o_idx = self.find_atom_in_residue_idx(res_idx, "O")

            self.residue_i[res_idx] = i
            self.trace.residue_indices[i] = res_idx
            self.trace.objids[i] = self.residue_objids[res_idx]
            self.trace.points[i] = self.soup.get_atom_proxy(ca_idx).pos
            self.trace.ups[i] = (
                self.soup.get_atom_proxy(c_idx).pos - self.soup.get_atom_proxy(o_idx).pos
            )

        for i in range(1, len(self.trace.points)):
            if v3.dot(self.trace.ups[i - 1], self.trace.ups[i]) < 0:
                self.trace.ups[i] = -self.trace.ups[i]

        self.center = self.get_center(self.trace.points)
        centered_points = [p - self.center for p in self.trace.points]

        self.scale = 1.0 / max(map(max, centered_points))

    def find_bb_hbonds(self):
        """
        Detect hydrogen bonds between backbone N and O atoms.

        Uses spatial hashing for efficient proximity detection. Identifies potential
        H-bonds when N and O atoms from different residues are within 3.5 Å.
        Stores H-bond partners as trace indices in residue_hb_partners dict.
        """
        print("Find H-Bonds...")
        vertices = []
        atoms = []
        for res_idx in self.trace.residue_indices:
            if res_idx is None:
                continue
            if self.has_atom_in_residue_idx(res_idx, "O"):
                o_atom_idx = self.find_atom_in_residue_idx(res_idx, "O")
                atoms.append(o_atom_idx)
                vertices.append(self.soup.get_atom_proxy(o_atom_idx).pos)
            if self.has_atom_in_residue_idx(res_idx, "N"):
                n_atom_idx = self.find_atom_in_residue_idx(res_idx, "N")
                atoms.append(n_atom_idx)
                vertices.append(self.soup.get_atom_proxy(n_atom_idx).pos)
            self.residue_hb_partners[res_idx] = []
        d = 3.5
        for i, j in SpaceHash(vertices).close_pairs():
            atom1_idx = atoms[i]
            atom2_idx = atoms[j]
            atom1_proxy = self.soup.get_atom_proxy(atom1_idx)
            atom2_proxy = self.soup.get_atom_proxy(atom2_idx)
            if atom1_proxy.atom_type == atom2_proxy.atom_type:
                continue
            if v3.pos_distance(atom1_proxy.pos, atom2_proxy.pos) < d:
                res1_idx = self.atom_residue_idx.get(atom1_idx)
                res2_idx = self.atom_residue_idx.get(atom2_idx)
                if res1_idx is not None and res2_idx is not None:
                    i1 = self.residue_i.get(res1_idx)
                    i2 = self.residue_i.get(res2_idx)
                    if i1 is not None and i2 is not None:
                        self.residue_hb_partners[res1_idx].append(i2)
                        self.residue_hb_partners[res2_idx].append(i1)

    def find_ss_by_bb_hbonds(self):
        """
        Assign secondary structure (SS) based on hydrogen bond patterns.

        Implements simplified DSSP-like algorithm:
        - Alpha helix ('H'): i->i+4 and i+1->i+5 H-bonds
        - 3-10 helix ('H'): i->i+3 and i+1->i+4 H-bonds
        - Beta sheet ('E'): parallel or anti-parallel H-bond pairs with distant residues
        - Coil ('C'): everything else

        Colors residues: helix=red, sheet=blue, coil=gray
        """

        def is_hb(i_res, j_res):
            if not (0 <= i_res <= len(self.trace.residue_indices) - 1):
                return False
            res_idx = self.trace.residue_indices[i_res]
            if res_idx is None:
                return False
            return j_res in self.residue_hb_partners.get(res_idx, [])

        print("Find Secondary Structure...")
        for res_idx in self.trace.residue_indices:
            if res_idx is not None:
                self.residue_ss[res_idx] = "C"

        n_res = len(self.trace.residue_indices)
        for i_res1 in range(n_res):
            if is_hb(i_res1, i_res1 + 4) and is_hb(i_res1 + 1, i_res1 + 5):
                for i_res in range(i_res1 + 1, i_res1 + 5):
                    res_idx = self.trace.residue_indices[i_res]
                    if res_idx is not None:
                        self.residue_ss[res_idx] = "H"

            if is_hb(i_res1, i_res1 + 3) and is_hb(i_res1 + 1, i_res1 + 4):
                for i_res in range(i_res1 + 1, i_res1 + 4):
                    res_idx = self.trace.residue_indices[i_res]
                    if res_idx is not None:
                        self.residue_ss[res_idx] = "H"

            for i_res2 in range(n_res):
                if abs(i_res1 - i_res2) > 5:
                    if is_hb(i_res1, i_res2):
                        beta_residues = []

                        if is_hb(i_res1 - 2, i_res2 - 2):
                            beta_residues.extend(
                                [i_res1 - 2, i_res1 - 1, i_res1, i_res2 - 2, i_res2 - 1, i_res2]
                            )
                        if is_hb(i_res1 + 2, i_res2 + 2):
                            beta_residues.extend(
                                [i_res1 + 2, i_res1 + 1, i_res1, i_res2 + 2, i_res2 + 1, i_res2]
                            )

                        if is_hb(i_res1 - 2, i_res2 + 2):
                            beta_residues.extend(
                                [i_res1 - 2, i_res1 - 1, i_res1, i_res2 + 2, i_res2 + 1, i_res2]
                            )
                        if is_hb(i_res1 + 2, i_res2 - 2):
                            beta_residues.extend(
                                [i_res1 + 2, i_res1 + 1, i_res1, i_res2 - 2, i_res2 - 1, i_res2]
                            )

                        for i_res in beta_residues:
                            res_idx = self.trace.residue_indices[i_res]
                            if res_idx is not None:
                                self.residue_ss[res_idx] = "E"

        color_by_ss = {
            "-": (0.5, 0.5, 0.5),
            "C": (0.5, 0.5, 0.5),
            "H": (0.8, 0.4, 0.4),
            "E": (0.4, 0.4, 0.8),
        }
        for res_idx in self.trace.residue_indices:
            if res_idx is not None:
                ss = self.residue_ss.get(res_idx, "C")
                self.residue_color[res_idx] = color_by_ss[ss]

    def find_pieces(self, cutoff=5.5):
        """
        Split trace into continuous pieces where CA atoms are close.

        Breaks the trace when CA-CA distance exceeds cutoff (default 5.5Å),
        indicating chain breaks or missing residues. Smooths tangent and up
        vectors within each piece for better rendering. Each piece becomes
        a separate continuous ribbon segment.
        """
        self.pieces = []

        i = 0
        n_point = len(self.trace.points)

        for j in range(1, n_point + 1):
            is_new_piece = False
            if j == n_point:
                is_new_piece = True
            else:
                dist = v3.pos_distance(self.trace.points[j - 1], self.trace.points[j])
                if dist > cutoff:
                    is_new_piece = True

            if is_new_piece:
                for k in range(i, j):
                    if k == i:
                        tangent = self.trace.points[i + 1] - self.trace.points[i]
                    elif k == j - 1:
                        tangent = self.trace.points[k] - self.trace.points[k - 1]
                    else:
                        tangent = self.trace.points[k + 1] - self.trace.points[k - 1]
                    norm = np.linalg.norm(tangent)
                    if norm > 0:
                        self.trace.tangents[k] = tangent / norm
                    else:
                        self.trace.tangents[k] = tangent

                ups = []
                for k in range(i, j):
                    up = self.trace.ups[k]
                    if k > i:
                        up = up + self.trace.ups[k - 1]
                    elif k < j - 1:
                        up = up + self.trace.ups[k + 1]
                    tangent = self.trace.tangents[k]
                    perp = up - np.dot(up, tangent) * tangent
                    norm = np.linalg.norm(perp)
                    if norm > 0:
                        ups.append(perp / norm)
                    else:
                        ups.append(perp)
                self.trace.ups[i:j] = ups

                self.pieces.append(SubTrace(self.trace, i, j))

                i = j

    def find_bonds(self):
        all_atom_indices = list(range(self.soup.get_atom_count()))
        backbone_atoms_temp = backbone_atoms.copy()
        backbone_atoms_temp.remove("CA")

        self.draw_to_screen_atoms = []
        for atom_idx in all_atom_indices:
            proxy = self.soup.get_atom_proxy(atom_idx)
            if proxy.atom_type not in backbone_atoms_temp and proxy.elem != "H":
                self.draw_to_screen_atoms.append(atom_idx)

        vertices = [self.soup.get_atom_proxy(idx).pos for idx in self.draw_to_screen_atoms]
        self.bonds = []
        print("Finding bonds...")
        for i, j in SpaceHash(vertices).close_pairs():
            atom1_idx = self.draw_to_screen_atoms[i]
            atom2_idx = self.draw_to_screen_atoms[j]
            atom1_proxy = self.soup.get_atom_proxy(atom1_idx)
            atom2_proxy = self.soup.get_atom_proxy(atom2_idx)

            d = 2
            if atom1_proxy.elem == "H" or atom2_proxy.elem == "H":
                continue
            if v3.pos_distance(atom1_proxy.pos, atom2_proxy.pos) < d:
                if atom1_proxy.alt != " " and atom2_proxy.alt != " ":
                    if atom1_proxy.alt != atom2_proxy.alt:
                        continue
                bond = Bond(atom1_idx, atom2_idx)
                bond.tangent = atom2_proxy.pos - atom1_proxy.pos
                bond.up = v3.cross_product_vec(atom1_proxy.pos, bond.tangent)
                self.bonds.append(bond)

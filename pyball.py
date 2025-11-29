#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PyBall: An OpenGL ES-based protein structure viewer.

This module provides a complete interactive 3D visualization system for
protein structures from PDB files, featuring:
- Ball-and-stick rendering
- Cartoon (ribbon) representation
- Cylinder trace visualization
- C-alpha arrow representation
- Secondary structure detection (helices/sheets)
- Interactive rotation and zoom controls

Usage:
    python pyball.py <pdb_file>

Controls:
    Mouse drag: Rotate view
    Mouse wheel: Zoom
    's': Toggle sidechain display
    'q': Quit application
"""

import itertools
import math
import sys
from ctypes import c_float
from pprint import pprint

import numpy as np
import numpy.linalg as linalg
import OpenGL.GL as gl
from pdbstruct import parse
from pdbstruct import vector3d as v3
from pdbstruct.vector3d import Matrix3d, Vector3d
from vispy import app, gloo, scene
from vispy.gloo import IndexBuffer, Program, VertexBuffer
from vispy.scene import visuals
from vispy.scene.visuals import Text
from vispy.util.transforms import ortho, perspective, rotate, translate

import render
from spacehash import SpaceHash

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
            # last division includes the very last trace point
            if j == n_trace_point - 1:
                n += 1
            for k in range(n):
                l = offset + k
                self.points[l, :] = catmull_rom_spline(
                    k * delta,
                    trace.get_prev_point(i),
                    trace.points[i],
                    trace.points[j],
                    trace.get_next_point(j),
                )
                self.ups[l, :] = catmull_rom_spline(
                    k * delta,
                    trace.get_prev_up(i),
                    trace.ups[i],
                    trace.ups[j],
                    trace.get_next_up(j),
                )
                if k / float(n) < 0.5:
                    self.objids[l] = trace.objids[i]
                else:
                    self.objids[l] = trace.objids[i + 1]
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
        self.atom_objids = {}  # atom_idx -> objid
        self.atom_by_objid = {}  # objid -> atom_idx
        self.atom_residue_idx = {}  # atom_idx -> residue_idx

        self.residue_objids = {}  # residue_idx -> objid
        self.residue_ss = {}  # residue_idx -> ss string
        self.residue_color = {}  # residue_idx -> [r,g,b]
        self.residue_i = {}  # residue_idx -> trace index
        self.residue_hb_partners = {}  # residue_idx -> list

        self.build_objids()

        self.build_trace()

        self.bonds = []
        self.find_bonds()

        self.pieces = []
        self.find_pieces()

        # self.find_ss_by_zhang_skolnick()
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
        # Check if numpy array or Vector3d list
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

        # remove alternate conformation by looking for orphaned atoms
        atoms = [
            atom_idx
            for atom_idx in self.iter_atoms()
            if self.atom_residue_idx.get(atom_idx) is not None
        ]

        # make ups point in the same direction
        for i in range(1, len(self.trace.points)):
            if v3.dot(self.trace.ups[i - 1], self.trace.ups[i]) < 0:
                self.trace.ups[i] = -self.trace.ups[i]

        # find geometrical center of points
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
            # alpha-helix
            if is_hb(i_res1, i_res1 + 4) and is_hb(i_res1 + 1, i_res1 + 5):
                for i_res in range(i_res1 + 1, i_res1 + 5):
                    res_idx = self.trace.residue_indices[i_res]
                    if res_idx is not None:
                        self.residue_ss[res_idx] = "H"

            # 3-10 helix
            if is_hb(i_res1, i_res1 + 3) and is_hb(i_res1 + 1, i_res1 + 4):
                for i_res in range(i_res1 + 1, i_res1 + 4):
                    res_idx = self.trace.residue_indices[i_res]
                    if res_idx is not None:
                        self.residue_ss[res_idx] = "H"

            for i_res2 in range(n_res):
                if abs(i_res1 - i_res2) > 5:
                    if is_hb(i_res1, i_res2):
                        beta_residues = []

                        # parallel beta sheet pairs
                        if is_hb(i_res1 - 2, i_res2 - 2):
                            beta_residues.extend(
                                [i_res1 - 2, i_res1 - 1, i_res1, i_res2 - 2, i_res2 - 1, i_res2]
                            )
                        if is_hb(i_res1 + 2, i_res2 + 2):
                            beta_residues.extend(
                                [i_res1 + 2, i_res1 + 1, i_res1, i_res2 + 2, i_res2 + 1, i_res2]
                            )

                        # anti-parallel beta sheet pairs
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
                    # Normalize numpy array directly
                    norm = np.linalg.norm(tangent)
                    if norm > 0:
                        self.trace.tangents[k] = tangent / norm
                    else:
                        self.trace.tangents[k] = tangent

                ups = []
                # smooth then rotate
                for k in range(i, j):
                    up = self.trace.ups[k]
                    if k > i:
                        up = up + self.trace.ups[k - 1]
                    elif k < j - 1:
                        up = up + self.trace.ups[k + 1]
                    # Compute perpendicular with numpy (up - (up·tangent)*tangent)
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


def identity():
    return np.eye(4, dtype=np.float32)


class Camera:
    """View transformation manager for interactive 3D camera control."""

    def __init__(self):
        self.view = identity()
        self.model = identity()
        self.rotation = identity()
        self.projection = identity()
        self.zoom = 40
        self.center = (0, 0, 0, 0)
        self.view = np.dot(self.view, translate((0, 0, -self.zoom)))
        self.is_fog = True
        self.fog_near = -1
        self.fog_far = 50
        self.fog_color = [0, 0, 0]

    def recalc_projection(self):
        self.projection = perspective(25.0, self.width / float(self.height), 1.0, 50.0 + self.zoom)
        self.fog_near = self.zoom
        self.fog_far = 20.0 + self.zoom

    def resize(self, width, height):
        self.width, self.height = width, height
        self.size = [width, height]
        self.recalc_projection()

    def recalc_model(self):
        self.model = np.dot(self.translation, self.rotation)

    def rotate(self, phi_diff, theta_diff, psi_diff):
        self.rotation = np.dot(self.rotation, rotate(phi_diff, (0, 1, 0)))
        self.rotation = np.dot(self.rotation, rotate(theta_diff, (1, 0, 0)))
        self.rotation = np.dot(self.rotation, rotate(psi_diff, (0, 0, -1)))
        self.recalc_model()

    def rezoom(self, zoom_diff):
        self.zoom = max(10, self.zoom + zoom_diff)
        self.view = translate((0, 0, -self.zoom))
        self.recalc_projection()

    def set_center(self, center):
        self.center = center
        self.translation = translate((-self.center[0], -self.center[1], -self.center[2]))
        self.recalc_model()


class TriangleStore:
    """
    Efficient builder for indexed triangle meshes.

    Accumulates vertices (position, normal, color, objid) and triangle indices
    into flat numpy arrays for OpenGL rendering. Supports triangle strips for
    efficient geometry representation.
    """

    def __init__(self, n_vertex):
        self.data = np.zeros(
            n_vertex,
            [
                ("a_position", np.float32, 3),
                ("a_normal", np.float32, 3),
                ("a_color", np.float32, 3),
                ("a_objid", np.float32, 1),
            ],
        )
        self.i_vertex = 0
        self.n_vertex = n_vertex
        self.indices = []

    def add_vertex(self, vertex, normal, color, objid):
        self.data["a_position"][self.i_vertex, :] = vertex
        self.data["a_normal"][self.i_vertex, :] = normal
        self.data["a_color"][self.i_vertex, :] = color
        self.data["a_objid"][self.i_vertex] = objid
        self.i_vertex += 1

    def vertex_buffer(self):
        return gloo.VertexBuffer(self.data)

    def index_buffer(self):
        return gloo.IndexBuffer(self.indices)

    def setup_next_strip(self, indices):
        """
        Add triangular indices relative to self.i_vertex_in_buffer
        """
        indices = [i + self.i_vertex for i in indices]
        self.indices.extend(indices)


def group(lst, n):
    """
    Returns iterable of n-tuple from a list.Incomplete tuples discarded
    http://code.activestate.com/recipes/303060-group-a-list-into-sequential-n-tuples/
    >>> list(group(range(10), 3))
        [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    return zip(*[itertools.islice(lst, i, None, n) for i in range(n)])


def make_calpha_arrow_mesh(rendered_soup, trace, length=0.7, width=0.35, thickness=0.3):
    """
    Generate arrows at each CA position pointing along backbone direction.
    Useful for quick visualization of chain direction and secondary structure.
    """
    arrow = render.Arrow(length, width, thickness)

    n_point = len(trace.points)
    triangle_store = TriangleStore(n_point * len(arrow.indices))

    for i_point in range(n_point):
        orientate = arrow.get_orientate(trace.tangents[i_point], trace.ups[i_point], 1.0)

        for indices in group(arrow.indices, 3):
            points = [arrow.vertices[i] for i in indices]

            normal = v3.cross_product_vec(points[1] - points[0], points[0] - points[2])
            normal = np.dot(orientate[:3, :3], normal)  # Matrix-vector multiply with numpy

            res_idx = trace.residue_indices[i_point]
            color = (
                rendered_soup.residue_color.get(res_idx, [0.4, 1.0, 0.4])
                if res_idx is not None
                else [0.4, 1.0, 0.4]
            )
            for point in points:
                triangle_store.add_vertex(
                    np.dot(orientate[:3, :3], point) + trace.points[i_point],
                    normal,
                    color,
                    trace.objids[i_point],
                )

    return triangle_store.vertex_buffer()


def make_cylinder_trace_mesh(rendered_soup, pieces, coil_detail=4, radius=0.3):
    """
    Simple tube representation through CA atoms.
    Creates cylindrical segments connecting consecutive CA positions, colored by
    residue secondary structure. Faster to render than cartoons.
    """
    cylinder = render.Cylinder(coil_detail)

    n_point = sum(len(piece.points) for piece in pieces)
    triangle_store = TriangleStore(2 * n_point * cylinder.n_vertex)

    for piece in pieces:
        points = piece.points

        for i_point in range(len(points) - 1):
            tangent = 0.5 * (points[i_point + 1] - points[i_point])

            res_idx1 = piece.residue_indices[i_point]
            res_idx2 = piece.residue_indices[i_point + 1]
            color1 = (
                rendered_soup.residue_color.get(res_idx1, [0.4, 1.0, 0.4])
                if res_idx1 is not None
                else [0.4, 1.0, 0.4]
            )
            color2 = (
                rendered_soup.residue_color.get(res_idx2, [0.4, 1.0, 0.4])
                if res_idx2 is not None
                else [0.4, 1.0, 0.4]
            )

            up = piece.ups[i_point] + piece.ups[i_point + 1]

            orientate = cylinder.get_orientate(tangent, up, radius)
            triangle_store.setup_next_strip(cylinder.indices)
            for point, normal in zip(cylinder.points, cylinder.normals):
                triangle_store.add_vertex(
                    np.dot(orientate[:3, :3], point) + points[i_point],
                    np.dot(orientate[:3, :3], normal),
                    color1,
                    piece.objids[i_point],
                )

            orientate = cylinder.get_orientate(-tangent, up, radius)
            triangle_store.setup_next_strip(cylinder.indices)
            for point, normal in zip(cylinder.points, cylinder.normals):
                triangle_store.add_vertex(
                    np.dot(orientate[:3, :3], point) + points[i_point + 1],
                    np.dot(orientate[:3, :3], normal),
                    color2,
                    piece.objids[i_point + 1],
                )

    return triangle_store.index_buffer(), triangle_store.vertex_buffer()


def make_carton_mesh(
    rendered_soup, pieces, coil_detail=5, spline_detail=3, width=1.6, thickness=0.2
):
    """
    Cartoon/ribbon representation with secondary structure-based geometry.

    Interpolates smooth splines through CA atoms then extrudes profiles:
    - Rectangles (ribbons) for helices and sheets
    - Circles (tubes) for coils
    Segments change geometry at SS boundaries for visual distinction.
    """
    rect = render.RectProfile(width, 0.15)
    circle = render.CircleProfile(coil_detail, 0.3)

    builders = []
    for piece in pieces:
        spline = SplineTrace(piece, 2 * spline_detail)

        n_point = len(piece.points)

        i_point = 0
        j_point = 1
        while i_point < n_point:
            res_idx = piece.residue_indices[i_point]
            ss = rendered_soup.residue_ss.get(res_idx, "C") if res_idx is not None else "C"
            color = (
                rendered_soup.residue_color.get(res_idx, [0.4, 1.0, 0.4])
                if res_idx is not None
                else [0.4, 1.0, 0.4]
            )
            color = [min(1.0, 1.2 * c) for c in color]
            profile = circle if ss == "C" else rect

            while j_point < n_point:
                res_idx_j = piece.residue_indices[j_point]
                ss_j = (
                    rendered_soup.residue_ss.get(res_idx_j, "C") if res_idx_j is not None else "C"
                )
                if ss_j != ss:
                    break
                j_point += 1

            i_spline = 2 * i_point * spline_detail - spline_detail
            if i_spline < 0:
                i_spline = 0
            j_spline = (j_point - 1) * 2 * spline_detail + spline_detail + 1
            if j_spline > len(spline.points) - 1:
                j_spline = len(spline.points) - 1

            sub_spline = SubTrace(spline, i_spline, j_spline)

            builders.append(render.TubeBuilder(sub_spline, profile, color))

            i_point = j_point
            j_point = i_point + 1

    n_vertex = sum(r.n_vertex for r in builders)
    triangle_store = TriangleStore(n_vertex)

    for r in builders:
        r.build_triangles(triangle_store)

    return triangle_store.index_buffer(), triangle_store.vertex_buffer()


def make_ball_and_stick_mesh(rendered_soup, sphere_stack=5, sphere_arc=5, tube_arc=5, radius=0.2):
    """
    Classic ball-and-stick representation showing atoms and bonds.

    Renders non-backbone heavy atoms as spheres, connected by cylindrical bonds.
    Bonds detected via spatial proximity (<2Å between non-H atoms). Each atom/bond
    half colored by its residue's secondary structure.
    """
    sphere = render.Sphere(sphere_stack, sphere_arc)
    cylinder = render.Cylinder(4)

    n_vertex = len(rendered_soup.draw_to_screen_atoms) * sphere.n_vertex
    n_vertex += 2 * len(rendered_soup.bonds) * cylinder.n_vertex
    triangle_store = TriangleStore(n_vertex)

    for atom_idx in rendered_soup.draw_to_screen_atoms:
        atom_proxy = rendered_soup.soup.get_atom_proxy(atom_idx)
        res_idx = rendered_soup.atom_residue_idx.get(atom_idx)
        color = (
            rendered_soup.residue_color.get(res_idx, [0.4, 1.0, 0.4])
            if res_idx is not None
            else [0.4, 1.0, 0.4]
        )
        objid = rendered_soup.atom_objids.get(atom_idx, atom_idx)

        triangle_store.setup_next_strip(sphere.indices)
        orientate = sphere.get_orientate(radius)
        for point in sphere.points:
            triangle_store.add_vertex(
                np.dot(orientate[:3, :3], point) + atom_proxy.pos,
                point,  # same as normal!
                color,
                objid,
            )

    for bond in rendered_soup.bonds:
        atom1_proxy = rendered_soup.soup.get_atom_proxy(bond.atom1_idx)
        atom2_proxy = rendered_soup.soup.get_atom_proxy(bond.atom2_idx)

        res1_idx = rendered_soup.atom_residue_idx.get(bond.atom1_idx)
        res2_idx = rendered_soup.atom_residue_idx.get(bond.atom2_idx)
        color1 = (
            rendered_soup.residue_color.get(res1_idx, [0.4, 1.0, 0.4])
            if res1_idx is not None
            else [0.4, 1.0, 0.4]
        )
        color2 = (
            rendered_soup.residue_color.get(res2_idx, [0.4, 1.0, 0.4])
            if res2_idx is not None
            else [0.4, 1.0, 0.4]
        )
        objid1 = rendered_soup.atom_objids.get(bond.atom1_idx, bond.atom1_idx)
        objid2 = rendered_soup.atom_objids.get(bond.atom2_idx, bond.atom2_idx)

        tangent = bond.tangent.scale(0.5)

        orientate = cylinder.get_orientate(tangent, bond.up, radius)
        triangle_store.setup_next_strip(cylinder.indices)
        for point, normal in zip(cylinder.points, cylinder.normals):
            triangle_store.add_vertex(
                np.dot(orientate[:3, :3], point) + atom1_proxy.pos,
                np.dot(orientate[:3, :3], normal),
                color1,
                objid1,
            )

        orientate = cylinder.get_orientate(-tangent, bond.up, radius)
        triangle_store.setup_next_strip(cylinder.indices)
        for point, normal in zip(cylinder.points, cylinder.normals):
            triangle_store.add_vertex(
                np.dot(orientate[:3, :3], point) + atom2_proxy.pos,
                np.dot(orientate[:3, :3], normal),
                color2,
                objid2,
            )

    return triangle_store.index_buffer(), triangle_store.vertex_buffer()


semilight_vertex = """
uniform mat4 u_model;
uniform mat4 u_normal;
uniform mat4 u_view;
uniform mat4 u_projection;

attribute vec3  a_position;
attribute vec3  a_normal;
attribute vec3  a_color;
attribute float a_objid;

varying vec4 N;
varying vec3 v_color;

void main (void)
{
    gl_Position = u_projection * u_view * u_model * vec4(a_position, 1.0);
    N = normalize(u_normal * vec4(a_normal, 1.0));
    v_color = a_color;
}
"""


semilight_fragment = """
uniform bool u_is_lighting;
uniform vec3 u_light_position;
uniform bool u_is_fog;
uniform float u_fog_near;
uniform float u_fog_far;
uniform vec3 u_fog_color;

const vec4 ambient_color = vec4(.2, .2, .2, 1.);
const vec4 diffuse_intensity = vec4(1., 1., 1., 1.); 

varying vec4 N;
varying vec3 v_color;

void main()
{
    if (u_is_lighting) {
        vec4 color = vec4(v_color, 1.0);
        vec4 L = vec4(normalize(u_light_position.xyz), 1);
        vec4 ambient = color * ambient_color;
        vec4 diffuse = color * diffuse_intensity;
        float d = max(0., dot(N, L));
        color = clamp(ambient + diffuse * d, 0., 1.);
        gl_FragColor = color;
    }

    if (u_is_fog) {
        float depth = gl_FragCoord.z / gl_FragCoord.w;
        float fog_factor = smoothstep(u_fog_near, u_fog_far, depth);
        gl_FragColor = mix(
            gl_FragColor, 
            vec4(u_fog_color, gl_FragColor.w), 
            fog_factor);
    }
}
"""


picking_vertex = """

uniform mat4 u_model;
uniform mat4 u_normal;
uniform mat4 u_view;
uniform mat4 u_projection;

attribute vec3 a_position;
attribute vec3 a_normal;
attribute vec3 a_color;
attribute float a_objid;

varying float v_objid;

void main(void) {
    gl_Position = u_projection * u_view * u_model * vec4(a_position, 1.0);
    v_objid = a_objid;
}
"""


picking_fragment = """

varying float v_objid;

int int_mod(int x, int y) { 
    int z = x / y;
    return x - y*z;
}

void main(void) {
    // ints are only required to be 7bit...
    int int_objid = int(v_objid + 0.5);
    int red = int_mod(int_objid, 256);
    int_objid /= 256;
    int green = int_mod(int_objid, 256);
    int_objid /= 256;
    int blue = int_mod(int_objid, 256);
    gl_FragColor = vec4(float(red), float(green), float(blue), 255.0)/255.0;
}


"""


def get_polar(x, y):
    r = math.sqrt(x * x + y * y)
    if x != 0.0:
        theta = math.atan(y / float(x))
    else:
        if y > 0:
            theta = math.pi / 2
        else:
            theta = -math.pi / 2
    if x < 0:
        if y > 0:
            theta += math.pi
        else:
            theta -= math.pi
    return r, theta


class Console:
    """Debug text overlay for displaying atom information on hover."""

    def __init__(self, size, init_str=""):
        self.text = Text(
            init_str,
            bold=True,
            color=(0.7, 1.0, 0.3, 1.0),
            font_size=10,
            pos=(0, 0),
            anchor_y="bottom",
            anchor_x="center",
        )
        self.size = size
        self.x = 0
        self.y = 0

    def draw(self):
        viewport = gl.glGetIntegerv(gl.GL_VIEWPORT)
        size = viewport[2:4]
        x_view_offset = (size[0] - self.size[0]) // 2
        y_view_offset = (size[1] - self.size[1]) // 2
        x = self.x + x_view_offset
        y = self.y + 15 + y_view_offset
        gl.glViewport(x, y, self.size[0], self.size[1])
        self.text.pos = (0, 0)
        self.text.draw()


class MolecularViewerCanvas(app.Canvas):
    """
    Main application window with interactive 3D protein visualization.

    Provides multiple rendering modes (cartoon, cylinders, arrows, ball-and-stick),
    interactive camera control (rotate/zoom), and atom picking for inspection.
    Press 's' to toggle sidechains, 'q' to quit.
    """

    def __init__(self, fname):
        app.Canvas.__init__(self, title="Molecular viewer", keys="interactive")
        size = (800, 600)
        self.size = size
        self.program = Program(semilight_vertex, semilight_fragment)
        self.picking_program = Program(picking_vertex, picking_fragment)

        soup = parse.load_soup(fname)
        rendered_soup = RenderedSoup(soup)
        self.rendered_soup = rendered_soup

        print("Building arrows...")
        self.arrow_buffer = make_calpha_arrow_mesh(rendered_soup, rendered_soup.trace)

        print("Building cylindrical trace...")
        self.cylinder_index_buffer, self.cylinder_vertex_buffer = make_cylinder_trace_mesh(
            rendered_soup, rendered_soup.pieces
        )

        print("Building cartoon...")
        self.cartoon_index_buffer, self.cartoon_vertex_buffer = make_carton_mesh(
            rendered_soup, rendered_soup.pieces
        )

        print("Building ball&sticks...")
        self.ballstick_index_buffer, self.ballstick_vertex_buffer = make_ball_and_stick_mesh(
            rendered_soup
        )

        self.draw_style = "sidechains"
        self.camera = Camera()
        self.camera.resize(*size)
        self.camera.set_center(rendered_soup.center)
        self.camera.rezoom(2.0 / rendered_soup.scale)
        self.new_camera = Camera()
        self.n_step_animate = 0
        self.console = Console(size)
        self.text = self.console.text
        self.timer = app.Timer(1.0 / 30)
        self.timer.connect(self.on_timer)
        self.timer.start()

        print(
            f"Camera center: {rendered_soup.center}, scale: {rendered_soup.scale}, zoom: {self.camera.zoom}"
        )

        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glDepthFunc(gl.GL_LESS)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glClearColor(0.0, 0.0, 0.0, 1.0)

        self.draw_style = "sidechains"
        self.camera = Camera()
        self.camera.resize(*size)
        self.camera.set_center(rendered_soup.center)
        self.camera.rezoom(2.0 / rendered_soup.scale)
        self.new_camera = Camera()
        self.n_step_animate = 0
        self.console = Console(size)
        self.text = self.console.text
        self.timer = app.Timer(1.0 / 30)
        self.timer.connect(self.on_timer)
        self.timer.start()

        print(
            f"Camera center: {rendered_soup.center}, scale: {rendered_soup.scale}, zoom: {self.camera.zoom}"
        )

        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glDepthFunc(gl.GL_LESS)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glClearColor(0.0, 0.0, 0.0, 1.0)

    def on_key_press(self, event):
        if event.text == " ":
            if self.timer.running:
                self.timer.stop()
            else:
                self.timer.start()
        if event.text == "s":
            if self.draw_style == "sidechains":
                self.draw_style = "no-sidechains"
            else:
                self.draw_style = "sidechains"
            self.update()
        if event.text == "c":
            self.draw_style = "no-sidechains"
            self.update()
        if event.text == "b":
            self.draw_style = "sidechains"
            self.update()

    def on_mouse_press(self, event):
        self.mouse_press_pos = event.pos

    def on_mouse_move(self, event):
        if event.is_dragging and hasattr(self, "mouse_press_pos"):
            dx = event.pos[0] - self.mouse_press_pos[0]
            dy = event.pos[1] - self.mouse_press_pos[1]
            self.camera.rotate(dx * 0.5, dy * 0.5, 0)
            self.mouse_press_pos = event.pos
            self.update()

    def on_mouse_wheel(self, event):
        self.camera.rezoom(-event.delta[1] * 2)
        self.update()

    def on_timer(self, event):
        if self.n_step_animate > 0:
            diff = self.new_camera.center - self.camera.center
            fraction = 1.0 / float(self.n_step_animate)
            new_center = self.camera.center + fraction * diff
            self.camera.set_center(new_center)
            self.n_step_animate -= 1
            self.update()

    def on_draw(self, event):
        width, height = self.physical_size

        if width != self.camera.width or height != self.camera.height:
            self.camera.resize(width, height)

        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glViewport(0, 0, width, height)

        self.program["u_light_position"] = [100.0, 100.0, 500.0]
        self.program["u_is_lighting"] = True
        self.program["u_model"] = self.camera.model
        self.program["u_normal"] = self.camera.rotation
        self.program["u_view"] = self.camera.view
        self.program["u_projection"] = self.camera.projection
        self.program["u_is_fog"] = self.camera.is_fog
        self.program["u_fog_far"] = self.camera.fog_far
        self.program["u_fog_near"] = self.camera.fog_near
        self.program["u_fog_color"] = self.camera.fog_color

        self.draw_buffers(self.program)

        self.console.draw()
        self.last_draw = "screen"

    def draw_buffers(self, program):
        if self.draw_style == "sidechains":
            program.bind(self.ballstick_vertex_buffer)
            program.draw("triangles", self.ballstick_index_buffer)
        elif self.draw_style == "no-sidechains":
            program.bind(self.cylinder_vertex_buffer)
            program.draw("triangles", self.cylinder_index_buffer)
        else:
            raise Exception("Unknown draw style")


def main(fname):
    print(f"Creating canvas for {fname}...")
    mvc = MolecularViewerCanvas(fname)
    print("Canvas created, showing window...")
    mvc.show()
    print("Window shown, starting event loop...")
    app.run()
    print("Event loop finished")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: pyball.py pdb")
    else:
        main(sys.argv[1])

"""
Rendering components for 3D protein visualization.

Contains camera control, mesh generation, shader programs, and rendering
functions for different protein representations (arrows, cylinders, cartoons,
ball-and-stick).
"""

import itertools
import math

import numpy as np
from pdbstruct import vector3d as v3
from vispy import gloo
from vispy.util.transforms import perspective, rotate, translate

from . import geometry
from .structures import SplineTrace, SubTrace


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
    arrow = geometry.Arrow(length, width, thickness)

    n_point = len(trace.points)
    triangle_store = TriangleStore(n_point * len(arrow.indices))

    for i_point in range(n_point):
        orientate = arrow.get_orientate(trace.tangents[i_point], trace.ups[i_point], 1.0)

        for indices in group(arrow.indices, 3):
            points = [arrow.vertices[i] for i in indices]

            normal = v3.cross_product_vec(points[1] - points[0], points[2] - points[0])
            normal = np.dot(orientate[:3, :3], normal)  # Matrix-vector multiply with numpy

            res_idx = trace.residue_indices[i_point]
            color = (
                rendered_soup.color_by_residue_idx.get(res_idx, [0.4, 1.0, 0.4])
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
    cylinder = geometry.Cylinder(coil_detail)

    n_point = sum(len(piece.points) for piece in pieces)
    triangle_store = TriangleStore(2 * n_point * cylinder.n_vertex)

    for piece in pieces:
        points = piece.points

        for i_point in range(len(points) - 1):
            tangent = 0.5 * (points[i_point + 1] - points[i_point])

            res_idx1 = piece.residue_indices[i_point]
            res_idx2 = piece.residue_indices[i_point + 1]
            color1 = (
                rendered_soup.color_by_residue_idx.get(res_idx1, [0.4, 1.0, 0.4])
                if res_idx1 is not None
                else [0.4, 1.0, 0.4]
            )
            color2 = (
                rendered_soup.color_by_residue_idx.get(res_idx2, [0.4, 1.0, 0.4])
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
    rect = geometry.RectProfile(width, 0.15)
    circle = geometry.CircleProfile(coil_detail, 0.3)

    builders = []
    for piece in pieces:
        spline = SplineTrace(piece, 2 * spline_detail)

        n_point = len(piece.points)

        i_point = 0
        j_point = 1
        while i_point < n_point:
            res_idx = piece.residue_indices[i_point]
            ss = rendered_soup.ss_by_residue_idx.get(res_idx, "C") if res_idx is not None else "C"
            color = (
                rendered_soup.color_by_residue_idx.get(res_idx, [0.4, 1.0, 0.4])
                if res_idx is not None
                else [0.4, 1.0, 0.4]
            )
            color = [min(1.0, 1.2 * c) for c in color]
            profile = circle if ss == "C" else rect

            while j_point < n_point:
                res_idx_j = piece.residue_indices[j_point]
                ss_j = (
                    rendered_soup.ss_by_residue_idx.get(res_idx_j, "C")
                    if res_idx_j is not None
                    else "C"
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

            builders.append(geometry.TubeBuilder(sub_spline, profile, color))

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
    sphere = geometry.Sphere(sphere_stack, sphere_arc)
    cylinder = geometry.Cylinder(4)

    n_vertex = len(rendered_soup.draw_to_screen_atoms) * sphere.n_vertex
    n_vertex += 2 * len(rendered_soup.bonds) * cylinder.n_vertex
    triangle_store = TriangleStore(n_vertex)

    for atom_idx in rendered_soup.draw_to_screen_atoms:
        atom_proxy = rendered_soup.soup.get_atom_proxy(atom_idx)
        res_idx = rendered_soup.residue_idx_by_atom_idx.get(atom_idx)
        color = (
            rendered_soup.color_by_residue_idx.get(res_idx, [0.4, 1.0, 0.4])
            if res_idx is not None
            else [0.4, 1.0, 0.4]
        )
        objid = rendered_soup.objid_by_atom_idx.get(atom_idx, atom_idx)

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

        res1_idx = rendered_soup.residue_idx_by_atom_idx.get(bond.atom1_idx)
        res2_idx = rendered_soup.residue_idx_by_atom_idx.get(bond.atom2_idx)
        color1 = (
            rendered_soup.color_by_residue_idx.get(res1_idx, [0.4, 1.0, 0.4])
            if res1_idx is not None
            else [0.4, 1.0, 0.4]
        )
        color2 = (
            rendered_soup.color_by_residue_idx.get(res2_idx, [0.4, 1.0, 0.4])
            if res2_idx is not None
            else [0.4, 1.0, 0.4]
        )
        objid1 = rendered_soup.objid_by_atom_idx.get(bond.atom1_idx, bond.atom1_idx)
        objid2 = rendered_soup.objid_by_atom_idx.get(bond.atom2_idx, bond.atom2_idx)

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

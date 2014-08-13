# -*- coding: utf-8 -*-

"""
Classes and objects that render geometric objects into
polygon triangle strips to form meshes for OpenGL display.
"""

import math

from ctypes import sizeof, c_float, c_void_p, c_uint

from pdbremix import v3numpy as v3

from OpenGL import GL as gl



class TriangleVertexBuffer:
  """
  Creates an indexed interleaved vertex buffer consisting of vertex, 
  normal and colors. Used to display triangle strips or triangles.

  This is a draw_object.

  Works closely with the shaders in shader.py as the
  variable name for the vertex attributes must be consistent.
  """
  
  def __init__(self, n_vertex):
    self.vertex_attribs = ['vertex', 'normal', 'color', 'objectid']
    self.n_floats = { 'vertex': 3, 'normal': 3, 'color': 3, 'objectid': 1 }
    self.offsets = { 'vertex': 0, 'normal': 3, 'color': 6 , 'objectid': 9 }
    self.n_float_in_vertex = sum(self.n_floats.values())
    self.size_vertex = self.n_float_in_vertex * sizeof(c_float)

    self.n_vertex = n_vertex
    self.size_vertex_buffer = self.n_float_in_vertex*n_vertex
    self.i_vertex_in_buffer = 0
    self.vertex_buffer = (c_float*self.size_vertex_buffer)()
    self.size_strip_list = []
    self.indices = []
    self.index_buffer = None

  def setup_next_strip(self, indices):
    """
    Add triangular indices relative to self.i_vertex_in_buffer
    """
    indices = [i + self.i_vertex_in_buffer for i in indices]
    self.indices.extend(indices)


  def add_vertex(self, vertex, normal, color, objectid):
    if self.i_vertex_in_buffer >= self.n_vertex:
      raise "Buffer overflow"
    i = self.i_vertex_in_buffer*self.n_float_in_vertex
    self.vertex_buffer[i  ] = vertex[0]
    self.vertex_buffer[i+1] = vertex[1]
    self.vertex_buffer[i+2] = vertex[2]
    self.vertex_buffer[i+3] = normal[0]
    self.vertex_buffer[i+4] = normal[1]
    self.vertex_buffer[i+5] = normal[2]
    self.vertex_buffer[i+6] = color[0]
    self.vertex_buffer[i+7] = color[1]
    self.vertex_buffer[i+8] = color[2]
    self.vertex_buffer[i+9] = objectid
    self.i_vertex_in_buffer += 1

  def bind_to_draw_context(self):
    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, gl.glGenBuffers(1))
    gl.glBufferData(gl.GL_ARRAY_BUFFER, self.vertex_buffer, gl.GL_STATIC_DRAW)

    for attrib in self.vertex_attribs:
      i = self.shader.locations[attrib]
      gl.glEnableVertexAttribArray(i)
      offset = self.offsets[attrib] * sizeof(c_float)
      gl.glVertexAttribPointer(
          i, self.n_floats[attrib], gl.GL_FLOAT, False, 
          self.size_vertex, c_void_p(offset))

    if self.index_buffer is None:
      self.index_buffer = (c_uint*len(self.indices))(*self.indices)

    gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, gl.glGenBuffers(1))
    gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER, self.index_buffer, gl.GL_STATIC_DRAW)

  def release_from_draw_context(self):
    for attrib in self.vertex_attribs:
      gl.glDisableVertexAttribArray(self.shader.locations[attrib])

  def draw_triangle_strips(self):
    n_index = len(self.indices)*3
    gl.glDrawElements(
        gl.GL_TRIANGLES,
        n_index, 
        gl.GL_UNSIGNED_INT,
        c_void_p(0))

  
  def draw(self, shader):
    self.shader = shader
    self.bind_to_draw_context()
    self.draw_triangle_strips()
    self.release_from_draw_context()



##################################################
# Profiles and extrusions for cartoon view

def calc_cyclic_normals(arcs, tangent):
  n = len(arcs)
  normals = []
  for i in range(len(arcs)):
    j = i-1 if i>0 else n-1
    k = i+1 if i<n-1 else 0
    arc_diff = arcs[k] - arcs[j]
    normals.append(v3.norm(v3.cross(arc_diff, tangent)))
  return normals


class CircleProfile():
  def __init__(self, n_arc=10, radius=1.0):
    self.n_arc = n_arc
    tangent = v3.vector(0, 0, 1)
    rotator = v3.vector(0, radius, 0)
    self.arcs = []
    angle = v3.radians(360/n_arc)
    rotation = v3.rotation(tangent, angle)
    self.arcs.append(v3.vector(rotator))
    for i_segment in range(n_arc-1):
      rotator = v3.transform(rotation, rotator)
      self.arcs.append(rotator)
    self.normals = calc_cyclic_normals(self.arcs, tangent)


class RectProfile():
  def __init__(self, width=1.5, thickness=0.3):
    tangent = v3.vector(0, 0, 1)
    up = v3.vector(0, width, 0)
    right = v3.cross(up, tangent)
    right = thickness/width*right
    self.arcs = [
        right + up,
              + up,
      - right + up, 
      - right - up,
              - up,
        right - up]
    self.normals = calc_cyclic_normals(self.arcs, tangent)


def get_xy_face_transform(tangent, up, scale):
  m = v3.identity()
  left = v3.cross(up, tangent)
  m[:3,0] = scale*v3.norm(left)
  m[:3,1] = scale*v3.norm(up)
  m[:3,2] = scale*v3.norm(tangent)
  return m


class TubeRender():
  def __init__(self, trace, profile, color):
    self.trace = trace
    self.profile = profile
    self.color = color
    n_arc = len(self.profile.arcs)
    n_slice = len(self.trace.points)
    self.n_vertex = n_arc*(n_slice + 2)

  def render_to_buffer(self, vertex_buffer):
    n_point = len(self.trace.points)
    n_arc = len(self.profile.arcs)

    # draw front face
    indices = []
    for i_arc in range((n_arc-1)/2):
      indices.extend([i_arc, i_arc+1, n_arc-1-i_arc])
      indices.extend([n_arc-1-i_arc, i_arc+1, n_arc-2-i_arc])
    vertex_buffer.setup_next_strip(indices)

    m = get_xy_face_transform(
        self.trace.tangents[0], self.trace.ups[0], 1.0)
    arcs = [
        v3.transform(m, a) + self.trace.points[0] 
        for a in self.profile.arcs]
    normal = -self.trace.tangents[0]
    for i_arc in range(n_arc):
      vertex_buffer.add_vertex(
          arcs[i_arc], normal, self.color, self.trace.objids[0])

    # draw extrusion in tube
    indices = []
    for i_point in range(n_point-1):
      i_slice_offset = i_point*n_arc
      j_slice_offset = (i_point+1)*n_arc
      for i_arc in range(n_arc):
        j_arc = (i_arc+1) % n_arc
        indices.append(i_slice_offset + i_arc)
        indices.append(j_slice_offset + i_arc)
        indices.append(i_slice_offset + j_arc)
        indices.append(i_slice_offset + j_arc)
        indices.append(j_slice_offset + i_arc)
        indices.append(j_slice_offset + j_arc)
    vertex_buffer.setup_next_strip(indices)

    for i in range(n_point):
      m = get_xy_face_transform(
          self.trace.tangents[i], self.trace.ups[i], 1.0)
      arcs = [
          v3.transform(m, a) + self.trace.points[i] 
          for a in self.profile.arcs]
      normals = [v3.transform(m, n) for n in self.profile.normals]
      for i_arc in range(n_arc):
        vertex_buffer.add_vertex(
            arcs[i_arc], normals[i_arc], self.color, self.trace.objids[i])

    # draw back face
    indices = []
    for i_arc in range((n_arc-1)/2):
      indices.extend([i_arc, i_arc+1, n_arc-1-i_arc])
      indices.extend([n_arc-1-i_arc, i_arc+1, n_arc-2-i_arc])
    vertex_buffer.setup_next_strip(indices)

    i_point = n_point - 1
    m = get_xy_face_transform(
        self.trace.tangents[i_point], self.trace.ups[i_point], 1.0)
    arcs = [
        v3.transform(m, a) + self.trace.points[i_point] 
        for a in self.profile.arcs]
    normal = self.trace.tangents[i_point]
    for i_arc in reversed(range(n_arc)):
      vertex_buffer.add_vertex(
          arcs[i_arc], normal, self.color, self.trace.objids[i_point])


##################################################
# Faces and polygons from standard 3D shapes

class ArrowShape():
  def __init__(self, length=2.0, width=0.8, thickness=1.0):
    arrow_face_in_zx = [
      v3.vector(0, 0, length),
      v3.vector(0, -width, -length),
      v3.vector(0, width, -length),
      ]

    self.points = []
    self.points.extend(
        [p + v3.vector(thickness, 0, 0) for p in arrow_face_in_zx])
    self.points.extend(
        [p + v3.vector(-thickness, 0, 0) for p in arrow_face_in_zx])

    self.point_indices = [0, 1, 2, 4, 3, 5]
    n_arc = 3
    for i in range(n_arc):
      j = (i+1) % n_arc
      self.point_indices.extend([
         i, i + n_arc, j, i + n_arc, j + n_arc, j ])
 
    self.n_triangle = 2 + 2*n_arc
    self.n_vertex = 3*self.n_triangle

  def render_to_center(
      self, vertex_buffer, center, tangent, up, scale, color, objid):
    m = get_xy_face_transform(tangent, up, scale)
    vertex_buffer.setup_next_strip(range(self.n_vertex))
    for i_triangle in range(self.n_triangle):
      triangle_indices = self.point_indices[i_triangle*3:i_triangle*3+3]
      triangle_points = [self.points[i] for i in reversed(triangle_indices)]
      p0, p1, p2 = triangle_points
      normal = v3.transform(m, v3.cross(p1 - p0, p0 - p2))
      for point in triangle_points:
        vertex = v3.transform(m, point) + center
        vertex_buffer.add_vertex(vertex, normal, color, objid)


class SphereShape:
  def __init__(self, n_stack=5, n_arc=5, scale=1.0):
    self.arcs = []
    self.stacks = []
    self.indices = []
    self.points = []
    vert_angle = math.pi / (n_stack - 1);
    arc_angle = math.pi * 2.0 / n_arc
    for i_stack in range(n_stack):
      radius = math.sin(i_stack * vert_angle);
      z = -scale * math.cos(i_stack * vert_angle);
      for i_arc in range(n_arc):
        x = scale * radius * math.cos(i_arc * arc_angle);
        y = scale * radius * math.sin(i_arc * arc_angle);
        self.points.append(v3.vector(x,y,z))
    for i_stack in range(n_stack-1):
      for i_arc in range(n_arc):
        j_arc = (i_arc+1) % n_arc
        self.indices.extend([
          (i_stack) * n_arc     + i_arc,
          (i_stack+1) * n_arc   + i_arc,
          (i_stack+1) * n_arc   + j_arc,
          (i_stack+1) * n_arc   + j_arc,
          (i_stack) * n_arc     + j_arc,
          (i_stack) * n_arc     + i_arc,
          ])
    self.n_vertex = n_stack*n_arc

  def render_to_center(
      self, vertex_buffer, center, tangent, up, scale, color, objid):
    m = v3.scaling_matrix(scale, scale, scale)
    vertex_buffer.setup_next_strip(self.indices)
    for i in range(self.n_vertex):
      normal = self.points[i]
      vertex = v3.transform(m, self.points[i]) + center
      vertex_buffer.add_vertex(vertex, normal, color, objid)


class CylinderShape:
  def __init__(self, n_arc, radius=1.0):
    self.points = []
    self.normals = []
    self.indices = []
    arc_angle = math.pi * 2.0 / n_arc
    for z in [0.0, 1.0]:
      for j in range(n_arc):
        x = radius * math.cos(j * arc_angle)
        y = radius * math.sin(j * arc_angle)
        self.normals.append(v3.vector(x, y, 0.0))
        self.points.append(v3.vector(x, y, z))
    for i_arc in range(n_arc):
      j_arc = (i_arc+1) % n_arc
      self.indices.extend(
          [j_arc, i_arc, n_arc + i_arc, 
           j_arc, n_arc + i_arc, n_arc + j_arc])
    self.n_vertex = 2*n_arc

  def render_to_center(
      self, vertex_buffer, center, tangent, up, scale, color, objid):
    v = v3.mag(tangent)/scale
    l = v3.scaling_matrix(1, 1, v)
    m = get_xy_face_transform(tangent, up, scale)
    n = v3.combine(m, l)
    vertex_buffer.setup_next_strip(self.indices)
    for i in range(self.n_vertex):
      vertex = v3.transform(n, self.points[i]) + center
      normal = v3.transform(n, self.normals[i])
      vertex_buffer.add_vertex(vertex, normal, color, objid)


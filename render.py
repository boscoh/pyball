# -*- coding: utf-8 -*-

"""
Classes and objects that render geometric objects into
polygon triangle strips to form meshes for OpenGL display.
"""

import math

import numpy as np

from pdbremix import v3numpy as v3

from vispy import gloo


class TriangleStore:
  def __init__(self, n_vertex):
    self.data = np.zeros(
      n_vertex, 
      [('a_position', np.float32, 3),
       ('a_normal', np.float32, 3),
       ('a_color', np.float32, 3),
       ('a_objid', np.float32, 1)])
    self.i_vertex = 0
    self.n_vertex = n_vertex
    self.indices = []

  def add_vertex(self, vertex, normal, color, objid):
    self.data['a_position'][self.i_vertex,:] = vertex
    self.data['a_normal'][self.i_vertex,:] = normal
    self.data['a_color'][self.i_vertex,:] = color
    self.data['a_objid'][self.i_vertex] = objid
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




def get_xy_face_transform(tangent, up, scale):
  m = v3.identity()
  left = v3.cross(up, tangent)
  m[:3,0] = scale*v3.norm(left)
  m[:3,1] = scale*v3.norm(up)
  m[:3,2] = scale*v3.norm(tangent)
  return m



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



class TubeBuilder():
  def __init__(self, trace, profile, color):
    self.trace = trace
    self.profile = profile
    self.color = color
    n_arc = len(self.profile.arcs)
    n_slice = len(self.trace.points)
    self.n_vertex = n_arc*(n_slice + 2)

  def build_triangles(self, vertex_buffer):
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




class Arrow:
  def __init__(self, length, width, thickness):
    self.vertices = [
      v3.vector( thickness,      0,  length),
      v3.vector( thickness, -width, -length),
      v3.vector( thickness,  width, -length),
      v3.vector(-thickness,      0,  length),
      v3.vector(-thickness, -width, -length),
      v3.vector(-thickness,  width, -length),
    ]
    self.indices = [1, 0, 2, 3, 4, 5]
    n_arc = 3
    for i in range(n_arc):
      j = (i+1) % n_arc
      self.indices.extend([i+n_arc, i, j, j+n_arc, i+n_arc, j])

  def get_orientate(self, tanget, up, scale):
    return get_xy_face_transform(tanget, up, scale)



class Cylinder:
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

  def get_orientate(self, tangent, up, scale):
    """
    Transform stretches the length to the length of tangent and
    the radius to scale.
    """
    return v3.combine(
        get_xy_face_transform(tangent, up, scale),
        v3.scaling_matrix(1.0, 1.0, v3.mag(tangent)/scale))



class Sphere:
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

  def get_orientate(self, scale):
    return v3.scaling_matrix(scale, scale, scale)



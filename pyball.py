#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
pyball - a pure-Python protein viewer in OpenGL ES 2.0 subset, 
must provide own shaders, matrices and lighting.

This implements a new ribbon representation with arrows at
every Calpha position.

Copyright (c) 2014, Bosco Ho

With ideas from Renaud Blanch's PyOpenGl tutorials
& Marco Biasini's javascript PV viewer.
"""


# imports ####################################################################

import os
import sys
from math import exp
from ctypes import sizeof, c_float, c_void_p, c_uint, string_at
from pprint import pprint

from OpenGL.GLUT import *
from OpenGL.GL import *

from pdbremix import pdbatoms
from pdbremix import v3



############################################
# Define camera and projection matrices


def matrix_to_c_floats(m, n_dim):
  n_float = n_dim*n_dim
  data = (c_float*n_float)()
  for i in range(n_dim):
    for j in range(n_dim):
      data[j*n_dim + i] = m[i,j]
  return data


def ortho(l, r, b, t, n, f):
  m = v3.identity()
  w, h, d = r-l, t-b, f-n
  m[0,:] = [2./w, 0.,    0.,   -(r+l)/w]
  m[1,:] = [0.,   2./h,  0.,   -(t+b)/h]
  m[2,:] = [0.,   0.,   -2./d, -(f+n)/d]
  m[3,:] = [0.,   0.,    0.,    1.     ]
  return m


def frustum(l, r, b, t, n, f):
  m = v3.identity()
  w, h, d = r-l, t-b, f-n
  m[0,:] = [2.*n/w, 0.,      (r+l)/w,  0.      ]
  m[1,:] = [0.,     2.*n/h,  (t+b)/h,  0.      ]
  m[2,:] = [0.,     0.,     -(f+n)/d, -2.*f*n/d]
  m[3,:] = [0.,     0.,     -1.,       0.      ]
  return m


class Camera:
  def __init__(self):
    self.rotation = v3.identity()
    self.center = v3.vector()
    self.scale = 1.0
    self.is_perspective = True
    self.is_lighting = True
    self.light_position = [100., 100., 500., 0.]
    self.half_width = 1.0
    self.half_height = 1.0

  def rotate_xy(self, rad_around_x, rad_around_y):
    rot_along_x = v3.rotation(v3.vector(0, 1, 0), rad_around_x)
    rot_along_y = v3.rotation(v3.vector(-1, 0, 0), rad_around_y)
    new_rotation = v3.combine(rot_along_x, rot_along_y)
    self.rotation = v3.combine(new_rotation, self.rotation)

  def rescale(self, new_scale):
    self.scale *= new_scale 
    
  def set_screen(self, width, height):
    self.half_width = width
    self.half_height = height

  def set_center(self, center):
    self.center = center

  def modelview_matrix(self):
    modelview = v3.identity()

    modelview = v3.translation(-self.center)
    modelview = v3.combine(self.rotation, modelview)

    s = self.scale
    scaling_matrix = v3.scaling_matrix(s, s, s)
    modelview = v3.combine(scaling_matrix, modelview)

    # w, h and values around 1.0 to indicate skew of rectangle
    w, h = self.half_width, self.half_height
    if self.is_perspective:
      projection = v3.combine(
          frustum(-w, w, -h, h, 8, 16), 
          v3.translation(v3.vector(0, 0, -12)))
      projection = v3.combine(
          projection, 
          v3.scaling_matrix(1.5, 1.5, 1.5))
    else:
      projection = ortho(-w, w, -h, h, -2, 2)
    modelview = v3.combine(projection, modelview)

    return matrix_to_c_floats(modelview, 4)

  def normal_matrix(self):
    # normal_matrix = m.transpose(m.inverse(modelview[:3,:3]))
    # since we only do uniform scaling, this reduces to
    return matrix_to_c_floats(self.rotation[:3,:3], 3)

  def bind_to_shader(self, shader):
    shader.set_matrix44(
        "modelview_matrix", self.modelview_matrix())
    shader.set_matrix33(
        "normal_matrix", self.normal_matrix())
    shader.set_boolean('lighting', self.is_lighting)
    shader.set_vec4('light_position', self.light_position)


#####################################################
# Vertex Buffer handler for triangle strips

class IndexedVertexBuffer:
  """
  Creates an indexed interleaved vertex buffer consisting of vertex, 
  normal and colors. Used to display triangle strips or triangles.
  """
  
  def __init__(self):
    self.vertex_attribs = ['vertex', 'normal', 'color']
    self.n_floats = { 'vertex': 3, 'normal': 3, 'color': 3 }
    self.offsets = { 'vertex': 0, 'normal': 3, 'color': 6 }
    self.size_vertex = 9 * sizeof(c_float)

    self.n_vertex = 0
    self.size_vertex_buffer = 0
    self.i_vertex_in_buffer = 0
    self.vertex_buffer = None
    self.size_strip_list = []
    self.indices = []
    self.index_buffer = None

  def set_size(self, n_vertex):
    self.n_vertex = n_vertex
    self.size_vertex_buffer = 9*n_vertex
    self.i_vertex_in_buffer = 0
    self.vertex_buffer = (c_float*self.size_vertex_buffer)()
    self.size_strip_list = []
    self.indices = []
    self.index_buffer = None

  def setup_next_strip(self, indices):
    self.size_strip_list.append(len(indices))
    if indices is None:
      indices = range(n_vertex)
    i_vertex = self.i_vertex_in_buffer
    indices = [i + i_vertex for i in indices]
    self.indices.extend(indices)

  def add_vertex(self, vertex, normal, color):
    if self.i_vertex_in_buffer >= self.n_vertex:
      raise "Buffer overflow"
    i = self.i_vertex_in_buffer*9
    self.vertex_buffer[i  ] = vertex[0]
    self.vertex_buffer[i+1] = vertex[1]
    self.vertex_buffer[i+2] = vertex[2]
    self.vertex_buffer[i+3] = normal[0]
    self.vertex_buffer[i+4] = normal[1]
    self.vertex_buffer[i+5] = normal[2]
    self.vertex_buffer[i+6] = color[0]
    self.vertex_buffer[i+7] = color[1]
    self.vertex_buffer[i+8] = color[2]
    self.i_vertex_in_buffer += 1

  def bind_to_draw_context(self):
    glBindBuffer(GL_ARRAY_BUFFER, glGenBuffers(1))
    glBufferData(GL_ARRAY_BUFFER, self.vertex_buffer, GL_STATIC_DRAW)

    for attrib in self.vertex_attribs:
      i = self.shader.locations[attrib]
      glEnableVertexAttribArray(i)
      offset = self.offsets[attrib] * sizeof(c_float)
      glVertexAttribPointer(
          i, self.n_floats[attrib], GL_FLOAT, False, 
          self.size_vertex, c_void_p(offset))

    if self.index_buffer is None:
      self.index_buffer = (c_uint*len(self.indices))(*self.indices)

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glGenBuffers(1))
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, self.index_buffer, GL_STATIC_DRAW)

  def release_from_draw_context(self):
    for attrib in self.vertex_attribs:
      glDisableVertexAttribArray(self.shader.locations[attrib])

  def draw_triangle_strips(self):
    offset_strip = 0
    for size_strip in self.size_strip_list:
      glDrawElements(
          GL_TRIANGLE_STRIP,
          size_strip, 
          GL_UNSIGNED_INT,
          c_void_p(offset_strip))
      offset_strip += size_strip*sizeof(c_uint)
  
  def draw(self, shader):
    self.shader = shader
    self.bind_to_draw_context()
    self.draw_triangle_strips()
    self.release_from_draw_context()



############################################
# Shader objects

vertex_shader_source = b"""

uniform mat4 modelview_matrix;
uniform mat4 texture_matrix;
uniform mat3 normal_matrix;

uniform bool lighting;
uniform vec4 light_position;

attribute vec3 vertex;
attribute vec3 normal;
attribute vec3 color;

varying vec3 N, L, S;

void main() {
  gl_Position = modelview_matrix * vec4(vertex, 1.);
  
  if (lighting) {
    N = normalize(normal_matrix*normal);
    L = normalize(light_position.xyz);
    S = normalize(L+vec3(0, 0, 1));
  }
  gl_FrontColor = vec4(color, 1.);
}

"""


fragment_shader_source = b"""

const vec4 acs = vec4(.2, .2, .2, 1.); // ambient color of scene
const vec4 di0 = vec4(1., 1., 1., 1.); // diffuse intensity of light 0
const vec4 white = vec4(1., 1., 1., 1.);
const float shininess = 100.;
const float alpha_threshold = .9;

uniform bool lighting;

varying vec3 N, L, S;

void main() {
  vec4 color = gl_Color;
  if(lighting) {
    vec4 ambient = color * acs;
    vec4 diffuse = color * di0;
    vec4 specular = white;
    float d = max(0., dot(N, L));
    float s = pow(max(0., dot(N, S)), shininess);
    color = clamp(ambient + diffuse * d + specular * s, 0., 1.);
  }
  gl_FragColor = color;
}

"""


class Shader():
  """
  Objects to hold compiled openGL shaders, with the
  locations to variables linked properly and accessed
  through the locations attribute.
  """
  def __init__(
      self, vertex_shader, fragment_shader,
      vertex_attribs, uniforms):

    self.program = glCreateProgram()

    glAttachShader(
       self.program,
       self.compile(GL_VERTEX_SHADER, vertex_shader))
    glAttachShader(
      self.program,
      self.compile(GL_FRAGMENT_SHADER, fragment_shader))

    self.locations = {}

    # before linking, link names of vertex attribs
    for i, vertex_attrib in enumerate(vertex_attribs):
      self.locations[vertex_attrib] = i
      glBindAttribLocation(self.program, i, vertex_attrib)

    glLinkProgram(self.program)
    if glGetProgramiv(self.program, GL_LINK_STATUS) != GL_TRUE:
      raise RuntimeError(glGetProgramInfoLog(self.program))

    # after linking, link names of uniforms
    for uniform in uniforms:
       i = glGetUniformLocation(self.program, uniform)
       self.locations[uniform] = i

  def compile(self, shader_type, source):
    shader = glCreateShader(shader_type)
    glShaderSource(shader, source)
    glCompileShader(shader)
    if glGetShaderiv(shader, GL_COMPILE_STATUS) != GL_TRUE:
      raise RuntimeError(glGetShaderInfoLog(shader))
    return shader

  def set_boolean(self, var_name, is_state):
    if is_state:
      glUniform1i(self.locations[var_name], 1)
    else:
      glUniform1i(self.locations[var_name], 0)

  def set_matrix33(self, var_name, c_float_9):
    glUniformMatrix3fv(
        self.locations[var_name], 
        1, # number of matrices
        GL_FALSE, # transpose
        c_float_9)

  def set_matrix44(self, var_name, c_float_16):
    glUniformMatrix4fv(
        self.locations[var_name], 
        1, # number of matrices
        GL_FALSE, # transpose
        c_float_16)

  def set_vec4(self, var_name, c_float_4):
    glUniform4f(self.locations[var_name], *c_float_4)


#########################################
# Main OpenGl Handler

class OpenGlHandler:
  """
  Class to bind arrays to OpenGL and tell it to display it.
  """

  def __init__(self):
    # name of vertex attributes used in shaders 
    self.vertex_attribs = [
      "vertex", 
      "normal", 
      "color"
    ]

    # name of uniform(global) variables used in shaders
    self.uniforms = [
      "lighting", 
      "light_position",
      "modelview_matrix", 
      "normal_matrix"
    ]

    self.shaders = {
      'default': 
        Shader(
          vertex_shader_source, 
          fragment_shader_source,
          self.vertex_attribs, 
          self.uniforms)
    }

    self.shader = self.shaders['default']
    self.camera = Camera()
    self.vertex_buffer = IndexedVertexBuffer()

  def setup_draw(self):
    glEnable(GL_DEPTH_TEST)
    glDepthFunc(GL_LEQUAL)

  def draw(self):
    self.setup_draw()
    glUseProgram(self.shader.program)
    self.camera.bind_to_shader(self.shader)
    self.vertex_buffer.draw(self.shader)


def save_to_png(name, width, height):
  size = width*height*3
    
  pixel_buffer = glGenBuffers(1)
  glBindBuffer(GL_PIXEL_PACK_BUFFER, pixel_buffer)
  glBufferData(GL_PIXEL_PACK_BUFFER, size, None, GL_STREAM_READ)
  
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, c_void_p(0))
  data = string_at(glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY), size)

  glUnmapBuffer(GL_PIXEL_PACK_BUFFER)
  glBindBuffer(GL_PIXEL_PACK_BUFFER, 0)
  glDeleteBuffers(1, [pixel_buffer])

  import png

  png.write(open(name, "wb"), width, height, 3, data)


#########################################
# GLUT Window Handler

class GlutWindowApp:
  def __init__(self, argv):
    self.width = 640
    self.height = 480

    self.is_mouse_left_down = False
    self.is_mouse_right_down  = False

    self.save_mouse_x = 0.0
    self.save_mouse_y = 0.0

    # glut initialization
    glutInit(argv)
    glutInitWindowSize(self.width, self.height)
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH)
    
    window_title = argv[0].encode()
    glutCreateWindow(window_title)
    
    glutReshapeFunc(self.reshape)
    glutDisplayFunc(self.display)
    glutKeyboardFunc(self.keyboard)
    glutMouseFunc(self.mouse)
    glutMotionFunc(self.motion)

    self.opengl = OpenGlHandler()

  def reshape(self, width, height):
    """window reshape callback."""
    glViewport(0, 0, width, height)
    radius = .5 * min(width, height)
    self.opengl.camera.set_screen(width/radius, height/radius)

  def display(self):
    """window redisplay callback."""
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    self.opengl.draw()
    glutSwapBuffers()

  def get_window_dims(self):
    self.width = glutGet(GLUT_WINDOW_WIDTH)
    self.height = glutGet(GLUT_WINDOW_HEIGHT)

  def screen_shot(self, name="screen_shot.png"):
    """window screenshot."""
    self.get_window_dims()
    save_to_png(name, self.width, self.height)

  def keyboard(self, c, x=0, y=0):
    """keyboard callback."""
    if c == b'p':
      is_perspective = not self.opengl.camera.is_perspective
      self.opengl.camera.is_perspective = is_perspective
      self.reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT))
    elif c == b'l':
      is_lighting = not self.opengl.camera.is_lighting
      self.opengl.camera.is_lighting = is_lighting
    elif c == b's':
      self.screen_shot()
    elif c == b'q':
      sys.exit(0)
    glutPostRedisplay()

  def screen2space(self, x, y):
    self.get_window_dims()
    radius = min(self.width, self.height)*self.opengl.camera.scale
    return (2.*x-self.width)/radius, -(2.*y-self.height)/radius

  def mouse(self, button, state, x, y):
    if button == GLUT_LEFT_BUTTON:
      self.is_mouse_left_down = (state == GLUT_DOWN)
      self.save_mouse_x, self.save_mouse_y = x, y
    elif button == GLUT_RIGHT_BUTTON:
      self.is_mouse_right_down = (state == GLUT_DOWN)
      self.save_mouse_x, self.save_mouse_y = x, y

  def motion(self, x1, y1):
    if self.is_mouse_left_down:
      old_x, old_y = self.screen2space(self.save_mouse_x, self.save_mouse_y)
      x, y = self.screen2space(x1, y1)
      self.opengl.camera.rotate_xy(0.1*(x-old_x), 0.1*(y-old_y))
    if self.is_mouse_right_down:
      diff = (x1-self.save_mouse_x)-(y1-self.save_mouse_y)
      new_scale = exp((diff)*.01)
      self.opengl.camera.rescale(new_scale)
    self.save_mouse_x, self.save_mouse_y = x1, y1
    glutPostRedisplay()

  def run(self):
    return glutMainLoop()


#########################################################
# Convert PDB structure into smooth pieces of secondary 
# structure

class PieceCAlphaTrace:
  def __init__(self, points, ups, tangents, ss):
    self.points = points
    self.ups = ups
    self.tangents = tangents
    self.ss = ss


class CAlphaTrace():
  def __init__(self, pdb):
    self.protein = pdbatoms.Soup(pdb)
    self.find_points()
    self.find_ss()
    self.find_ss_pieces()

  def find_points(self):
    self.points = []
    tops = []
    for residue in self.protein.residues():
      if residue.has_atom('CA'):
        self.points.append(residue.atom('CA').pos)
        tops.append(residue.atom('C').pos - residue.atom('O').pos)

    n_point = len(self.points)

    self.center = v3.get_center(self.points)
    centered_points = [p - self.center for p in self.points]
    self.scale = 1.0/max(map(max, centered_points))

    for i in range(1, n_point):
      if v3.dot(tops[i-1], tops[i]) < 0:
         tops[i] = -tops[i]

    self.tangents = []
    for i in range(0, n_point):
      if i == 0:
        tangent = self.points[1] - self.points[0]
      elif i == n_point-1:
        tangent = self.points[i] - self.points[i-1]
      else:
        tangent = self.points[i+1] - self.points[i-1]
      self.tangents.append(v3.norm(tangent))

    self.ups = []
    for i in range(n_point):
      up = tops[i]
      if i > 0:
        up = up + tops[i-1]
      elif i < n_point-1:
        up = up + tops[i+1]
      up = v3.perpendicular(up, self.tangents[i])
      self.ups.append(v3.norm(up))

  def find_ss(self):
    
    def zhang_skolnick_test(i, template_dists, delta):
      for j in range(max(0, i-2), i+1):
        for diff in range(2, 5):
          k = j + diff
          if k >= len(self.points):
            continue
          dist = v3.distance(self.points[j], self.points[k])
          if abs(dist - template_dists[diff]) > delta:
            return False
      return True

    helix_distances = { 2:5.45, 3:5.18, 4:6.37 }
    helix_delta = 2.1
    sheet_distances = { 2:6.1, 3:10.4, 4:13.0 }
    sheet_delta = 1.42
    self.ss = []
    for i in range(len(self.points)):
      if zhang_skolnick_test(i, helix_distances, helix_delta):
        self.ss.append('H')
      elif zhang_skolnick_test(i, sheet_distances, sheet_delta):
        self.ss.append('E')
      else:
        self.ss.append('C')

  def find_ss_pieces(self):
    self.pieces = []
    ss = self.ss[0]
    i = 0
    n_point = len(self.points)
    for j in range(1, n_point+1):
      is_new_piece = False
      is_break = False
      if j == n_point:
        is_new_piece = True
      elif j < n_point and self.ss[j] != ss:
        is_new_piece = True
      else:
        dist = v3.distance(self.points[j-1], self.points[j]) 
        cutoff = 5
        if dist > cutoff:
          is_new_piece = True
          is_break = True
      if is_new_piece:
        piece = PieceCAlphaTrace(
            self.points[i:j], 
            self.ups[i:j], 
            self.tangents[i:j],
            ss)
        self.pieces.append(piece)
        if j < n_point:
          ss = self.ss[j]
        if is_break:
          i = j
        else:
          i = j-1


def spline(n, p1, p2, p3, p4):
  """
  Returns n points between p2 and p3 that has a 
  tangent to p1 and p4 using Catmull-Rom spline.
  """
  vecs = []
  for i in range(n-1):
    t = i/float(n)
    vecs.append(
        (   t*((2-t)*t    - 1)  * p1
          + (t*t*(3*t - 5) + 2) * p2
          + t*((4 - 3*t)*t + 1) * p3
          + (t-1)*t*t           * p4 ) / 2
    )
  vecs.append(p3)
  return vecs


class ExpandedTrace():
  def __init__(self, trace, n_division):
    n_guide_point = len(trace.points)

    self.points = []
    for i in range(n_guide_point-1):
      if i == 0:
        prev_point = trace.points[0] - trace.tangents[0]
      else:
        prev_point = trace.points[i-1]
      if i == n_guide_point - 2:
        next_point = trace.points[i+1] + trace.tangents[-1]
      else:
        next_point = trace.points[i+2]
      spline_points = spline(
        n_division, prev_point, trace.points[i], trace.points[i+1], next_point)
      if i < n_guide_point - 2:
        self.points.extend(spline_points[:-1])
      else:
        self.points.extend(spline_points)

    self.tangents = []
    n_point = len(self.points)
    for i in range(n_point):
      if i == 0:
        tangent = trace.tangents[0]
      elif i == n_point-1:
        tangent = trace.tangents[-1]
      else:
        tangent = self.points[i+1] - self.points[i-1]
      self.tangents.append(tangent)

    self.ups = []
    for i in range(n_guide_point-1):
      if i == 0:
        prev_up = trace.ups[0]
      else:
        prev_up = trace.ups[i-1]
      if i == n_guide_point - 2:
        next_up = trace.ups[-1]
      else:
        next_up = trace.ups[i+2]
      spline_ups = spline(
        n_division, prev_up, trace.ups[i], trace.ups[i+1], next_up)
      if i < n_guide_point - 2:
        self.ups.extend(spline_ups[:-1])
      else:
        self.ups.extend(spline_ups)



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
  def __init__(self, n_arc=10):
    self.n_arc = n_arc
    tangent = v3.vector(0, 0, 1)
    rotator = v3.vector(0, 0.2, 0)
    self.arcs = []
    angle = v3.radians(360/n_arc)
    rotation = v3.rotation(tangent, angle)
    self.arcs.append(v3.vector(rotator))
    for i_segment in range(n_arc-1):
      rotator = v3.transform(rotation, rotator)
      self.arcs.append(rotator)
    self.normals = calc_cyclic_normals(self.arcs, tangent)


class RectProfile():
  def __init__(self):
    tangent = v3.vector(0, 0, 1)
    up = v3.vector(0, 1.2, 0)
    right = v3.cross(up, tangent)
    right = 0.1*right
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
  
  def get_n_vertex(self):
    n_arc = len(self.profile.arcs)
    n_slice = len(self.trace.points)
    return n_arc*n_slice

  def render(self, vertex_buffer):
    n_point = len(self.trace.points)
    n_arc = len(self.profile.arcs)

    indices = []
    for i_point in range(n_point-1):
      i_slice_offset = i_point*n_arc
      j_slice_offset = (i_point+1)*n_arc
      for i_arc in range(n_arc):
        if i_arc < n_arc-1:
          j_arc = i_arc + 1
        else:
          j_arc = 0
        indices.append(i_slice_offset + i_arc)
        indices.append(j_slice_offset + i_arc)
        indices.append(i_slice_offset + j_arc)
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
            arcs[i_arc], normals[i_arc], self.color)


class CartoonTraceRenderer():
  def __init__(self, trace, coil_detail=4, spline_detail=6):
    self.trace = trace
    rect = RectProfile()
    circle = CircleProfile(coil_detail)
    color_by_ss = {
      'C': (0.8, 0.8, 0.8),
      'H': (1.0, 0.6, 0.6),
      'E': (0.6, 0.6, 1.0)
    }
    self.piece_renderers = []
    for piece in trace.pieces:
      profile = circle if piece.ss == "C" else rect
      color = color_by_ss[piece.ss]
      trace_piece = ExpandedTrace(piece, spline_detail)
      self.piece_renderers.append(
          TubeRender(trace_piece, profile, color))

  def get_n_vertex(self):
    return sum(p.get_n_vertex() for p in self.piece_renderers)

  def render(self, vertex_buffer):
    for r in self.piece_renderers:
      r.render(vertex_buffer)


##################################################
# Basic Geometry function

class Arrow():
  def __init__(self, w=0.6):
    arrow_face_in_zx = [
      v3.vector(0, 0, 1.2),
      v3.vector(0, 0.4, -1.2),
      v3.vector(0, -0.4, -1.2)]

    i = 0
    self.points = []
    self.faces = []
    n_point_in_arrow_face = len(arrow_face_in_zx)
    points = [p + v3.vector(0.3, 0, 0) for p in arrow_face_in_zx]
    points = [p*w for p in points]
    self.points.extend(points)
    face = [i + j for j in range(len(points))] 
    self.faces.append(list(reversed(face)))
    i += n_point_in_arrow_face

    points = [p + v3.vector(-0.4, 0, 0) for p in arrow_face_in_zx]
    points = [p*w for p in points]
    self.points.extend(points)
    face = [i + j for j in range(len(points))] 
    self.faces.append(face)
    i += n_point_in_arrow_face

    for j in range(n_point_in_arrow_face):
      if j == n_point_in_arrow_face-1:
        k = 0
      else:
        k = j+1
      face = [k, k+n_point_in_arrow_face, j, j+n_point_in_arrow_face]
      self.faces.append(face)

  def render_to_center(
      self, vertex_buffer, center, tangent, up, scale, color):
    m = get_xy_face_transform(tangent, up, scale)
    for indexes in self.faces:
      n_vertex = len(indexes)
      vertex_buffer.setup_next_strip(range(n_vertex))
      p0, p1, p2 = [self.points[i] for i in indexes[:3]]
      normal = v3.transform(m, v3.cross(p0 - p1, p0 - p2))
      for index in indexes:
        vertex = v3.transform(m, self.points[index]) + center
        vertex_buffer.add_vertex(vertex, normal, color)


class ShapeTraceRenderer():
  def __init__(self, trace, shape):
    self.trace = trace
    self.shape = shape

  def get_n_vertex(self):
    n_vertex_in_shape = sum(len(f) for f in self.shape.faces)
    n_point = len(self.trace.points)
    return n_vertex_in_shape*n_point

  def render(self, vertex_buffer):
    n_point = len(self.trace.points)
    color_by_ss = {
      'C': (0.5, 0.5, 0.5),
      'H': (1.0, 0.4, 0.4),
      'E': (0.4, 0.4, 1.0)
    }
    scale = 1.0
    for i in range(n_point):
      self.shape.render_to_center(
          vertex_buffer,
          self.trace.points[i], 
          self.trace.tangents[i],
          self.trace.ups[i],
          scale,
          color_by_ss[self.trace.ss[i]])


def main():
  pdb = sys.argv[1] if len(sys.argv) > 1 else '1cph.pdb'
  c_alpha_trace = CAlphaTrace(pdb)

  app = GlutWindowApp(sys.argv)

  camera = app.opengl.camera
  camera.rescale(c_alpha_trace.scale)
  camera.set_center(c_alpha_trace.center)

  renderers = [
    CartoonTraceRenderer(c_alpha_trace),
    ShapeTraceRenderer(c_alpha_trace, Arrow()),
  ]

  n_all_vertex = sum(r.get_n_vertex() for r in renderers)
  vertex_buffer = app.opengl.vertex_buffer
  vertex_buffer.set_size(n_all_vertex)
  for r in renderers: 
    r.render(vertex_buffer)

  app.run()


if __name__ == "__main__":
  main()




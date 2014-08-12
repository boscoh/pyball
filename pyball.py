#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
pyball - a pure-Python protein viewer in OpenGL ES 2.0 subset, 
must provide own shaders, matrices and lighting.

This implements a new ribbon representation with arrows at
every Calpha position.

Copyright (c) 2014, Bosco Ho

With bits adapted from Renaud Blanch's PyOpenGl tutorials
& Marco Biasini's javascript PV viewer.
"""


import os
import sys
from math import exp
import math
from pprint import pprint
import time

from ctypes import c_float

import OpenGL.GLUT as glut
import OpenGL.GL as gl

from pdbremix import pdbatoms
from pdbremix import v3numpy as v3
from pdbremix.data import backbone_atoms

import camera
import shader
import render
from spacehash import SpaceHash



#########################################################
# Convert PDB structure into smooth pieces of secondary 
# structure and geometrial objects that use
# render functions to turn into polygon


next_objid = 1
def get_next_objid():
  global next_objid
  objid = next_objid
  next_objid += 1
  return objid


def atom_name(atom):
  return atom.res_tag() + '-' + atom.res_type  + '-' + atom.type


class PieceCalphaTrace:
  def __init__(self):
    self.points = []
    self.ups = []
    self.tangents = []
    self.objids = []

  def get_prev_point(self, i):
    if i > 0:
      return self.points[i-1]
    else:
      return self.points[i] - self.tangents[i]

  def get_next_point(self, i):
    if i < len(self.points)-1:
      return self.points[i+1]
    else:
      return self.points[i] + self.tangents[i]

  def get_prev_up(self, i):
    if i > 0:
      return self.ups[i-1]
    else:
      return self.ups[i]

  def get_next_up(self, i):
    if i < len(self.points)-1:
      return self.ups[i+1]
    else:
      return self.ups[i]


def spline(t, p1, p2, p3, p4):
  """
  Returns a point at fraction t between p2 and p3 
  using Catmull-Rom spline.
  """
  return \
      0.5 * (   t*((2-t)*t    - 1)  * p1
              + (t*t*(3*t - 5) + 2) * p2
              + t*((4 - 3*t)*t + 1) * p3
              + (t-1)*t*t           * p4 )


class SplineTrace(PieceCalphaTrace):
  """
  A Spline Trace is used to draw ribbons and tubes.
  It's essentially a collection of points, objids, tangents and ups.
  """
  def __init__(self, trace, n_division):
    PieceCalphaTrace.__init__(self)

    n_guide_point = len(trace.points)
    delta = 1/float(n_division)

    for i in range(n_guide_point-1):
      n = n_division
      j = i+1
      if j == n_guide_point - 1:
        n += 1
      for k in range(n):
        self.points.append(
          spline(
            k*delta, 
            trace.get_prev_point(i), 
            trace.points[i],
            trace.points[j], 
            trace.get_next_point(j)))
        if k/float(n) < 0.5:
          i_objid = i
        else:
          i_objid = i+1
        self.objids.append(trace.objids[i_objid])

    n_point = len(self.points)
    for i in range(n_point):
      if i == 0:
        tangent = trace.tangents[0]
      elif i == n_point-1:
        tangent = trace.tangents[-1]
      else:
        tangent = self.points[i+1] - self.points[i-1]
      self.tangents.append(tangent)

    for i in range(n_guide_point-1):
      if i == 0:
        prev_up = trace.ups[0]
      else:
        prev_up = trace.ups[i-1]
      if i == n_guide_point - 2:
        next_up = trace.ups[-1]
      else:
        next_up = trace.ups[i+2]
      n = n_division
      if i == n_guide_point - 2:
        n += 1
      for k in range(n):
        self.ups.append(
          spline(
             k*delta, 
             trace.get_prev_up(i), 
             trace.ups[i], 
             trace.ups[i+1], 
             trace.get_next_up(i)))


class Bond():
  def __init__(self, atom1, atom2):
    self.atom1 = atom1
    self.atom2 = atom2


class RenderedSoup():
  def __init__(self, soup):
    self.soup = soup
    self.bonds = []
    self.points = []
    self.objids = []
    self.tops = []
    self.pieces = []
    self.ss_pieces = []
    self.objid_ref = {}
    self.residues = []

    self.build_objids()
    self.find_points()
    self.find_ss()
    self.find_pieces()
    self.find_bonds()

  def build_objids(self):
    for atom in self.soup.atoms():
      objid = get_next_objid()
      self.objid_ref[objid] = atom
      atom.objid = objid
      atom.res_objid = 0
    for residue in self.soup.residues():
      if residue.has_atom('CA'):
        res_objid = residue.atom('CA').objid
      else:
        res_objid = residue.atoms()[0].objid
      residue.ss = '-'
      residue.color = [0.4, 1.0, 0.4]
      residue.objid = res_objid
      for atom in residue.atoms():
        atom.res_objid = atom.objid
        atom.residue = residue

  def find_points(self):
    for residue in self.soup.residues():
      if residue.has_atom('CA'):
        self.residues.append(residue)
        atom = residue.atom('CA')
        self.points.append(atom.pos)
        self.tops.append(residue.atom('C').pos - residue.atom('O').pos)
        self.objids.append(atom.objid)
    # FIXIT: remove the alternate conformation not in residue
    atoms = self.soup.atoms()
    n = len(atoms)
    for i in reversed(range(n)):
      atom = atoms[i]
      if not hasattr(atom, 'residue'):
        del atoms[i]
    n_point = len(self.points)
    for i in range(1, n_point):
      if v3.dot(self.tops[i-1], self.tops[i]) < 0:
         self.tops[i] = -self.tops[i]
    self.center = v3.get_center(self.points)
    centered_points = [p - self.center for p in self.points]
    self.scale = 1.0/max(map(max, centered_points))

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
    for i in range(len(self.points)):
      if zhang_skolnick_test(i, helix_distances, helix_delta):
        ss = 'H'
      elif zhang_skolnick_test(i, sheet_distances, sheet_delta):
        ss = 'E'
      else:
        ss = 'C'
      self.residues[i].ss = ss

    color_by_ss = {
      'C': (0.5, 0.5, 0.5),
      'H': (0.8, 0.4, 0.4),
      'E': (0.4, 0.4, 0.8)
    }
    for residue in self.residues:
      residue.color = color_by_ss[residue.ss]

  def find_pieces(self):
    self.pieces = []
    i = 0
    n_point = len(self.points)
    for j in range(1, n_point+1):
      is_new_piece = False
      is_break = False
      if j == n_point:
        is_new_piece = True
      else:
        dist = v3.distance(self.points[j-1], self.points[j]) 
        cutoff = 5.5
        if dist > cutoff:
          is_new_piece = True
      if is_new_piece:
        piece = PieceCalphaTrace()
        piece.points = self.points[i:j]
        piece.objids = self.objids[i:j]
        piece.residues = self.residues[i:j]
        n_point_piece = len(piece.points)
        for k in range(n_point_piece):
          if k == 0:
            tangent = piece.points[1] - piece.points[0]
          elif k == n_point_piece-1:
            tangent = piece.points[k] - piece.points[k-1]
          else:
            tangent = piece.points[k+1] - piece.points[k-1]
          piece.tangents.append(v3.norm(tangent))

        piece.ups = []
        for k in range(n_point_piece):
          k_full = k + i
          up = self.tops[k_full]
          if k > 0:
            up = up + self.tops[k_full-1]
          elif k < n_point_piece-1:
            up = up + self.tops[k_full+1]
          up = v3.perpendicular(up, piece.tangents[k])
          piece.ups.append(v3.norm(up))

        self.pieces.append(piece)
        i = j

  def find_bonds(self):
    self.draw_to_screen_atoms = self.soup.atoms()
    backbone_atoms.remove('CA')
    self.draw_to_screen_atoms = [a for a in self.draw_to_screen_atoms if a.type not in backbone_atoms and a.element!="H"]
    vertices = [a.pos for a in self.draw_to_screen_atoms]
    self.bonds = []
    print "Finding bonds..."
    for i, j in SpaceHash(vertices).close_pairs():
      atom1 = self.draw_to_screen_atoms[i]
      atom2 = self.draw_to_screen_atoms[j]
      d = 2
      if atom1.element == 'H' or atom2.element == 'H':
        continue
      if v3.distance(atom1.pos, atom2.pos) < d:
        if atom1.alt_conform != " " and atom2.alt_conform != " ":
          if atom1.alt_conform != atom2.alt_conform:
            continue
        bond = Bond(atom1, atom2)
        bond.tangent = atom2.pos - atom1.pos
        bond.up = v3.cross(atom1.pos, bond.tangent)
        self.bonds.append(bond)


def make_carton_mesh(
    rendered_soup, coil_detail=5, spline_detail=3, 
    width=1.6, thickness=0.2):

  # extrusion cross sections
  rect = render.RectProfile(width, thickness)
  circle = render.CircleProfile(coil_detail, thickness)

  renderers = []
  for piece in rendered_soup.pieces:
    spline = SplineTrace(piece, 2*spline_detail)

    n_residue = len(piece.points)
    i_residue = 0
    j_residue = 1
    while i_residue < n_residue:

      ss = piece.residues[i_residue].ss
      color = piece.residues[i_residue].color
      color = [min(1.0, 1.6*c) for c in color]
      profile = circle if ss == "C" else rect  

      while j_residue < n_residue and piece.residues[j_residue].ss == ss:
        j_residue += 1

      sub_spline = PieceCalphaTrace()
      i_spline = i_residue * 2*spline_detail - spline_detail
      if i_spline < 0:
        i_spline = 0
      j_spline = (j_residue-1) * 2*spline_detail + spline_detail + 1
      if j_spline > len(spline.points) - 1:
        j_spline = len(spline.points) - 1
      sub_spline.points = spline.points[i_spline:j_spline] 
      sub_spline.objids = spline.objids[i_spline:j_spline] 
      sub_spline.ups = spline.ups[i_spline:j_spline] 
      sub_spline.tangents = spline.tangents[i_spline:j_spline]

      renderers.append(
          render.TubeRender(sub_spline, profile, color))

      i_residue = j_residue
      j_residue = i_residue + 1

  n_vertex = sum(r.n_vertex for r in renderers)
  vertex_buffer = render.IndexedVertexBuffer(n_vertex)

  for r in renderers:
      r.render_to_buffer(vertex_buffer)

  return vertex_buffer


def make_arrow_mesh(rendered_soup):
  shape = render.ArrowShape(0.6, 0.3, 0.24)
  n_point = 0
  for piece in rendered_soup.pieces:
    n_point += len(piece.points)
  n_vertex = shape.n_vertex*n_point
  vertex_buffer = render.IndexedVertexBuffer(n_vertex)
  for piece in rendered_soup.pieces:
    n_point = len(piece.points)
    for i in range(n_point):
      shape.render_to_center(
          vertex_buffer,
          piece.points[i], 
          piece.tangents[i],
          piece.ups[i],
          1.0,
          piece.residues[i].color,
          piece.objids[i])
  return vertex_buffer


def make_cylinder_trace_mesh(
    rendered_soup, coil_detail=4, spline_detail=6, sphere_detail=4):
  cylinder = render.CylinderShape(coil_detail)
  n_point = sum(len(piece.points) for piece in rendered_soup.pieces)
  n_vertex = cylinder.n_vertex * 2* n_point
  vertex_buffer = render.IndexedVertexBuffer(n_vertex)
  for piece in rendered_soup.pieces:
    points = piece.points
    for i_segment in range(len(points) - 1):
      tangent = 0.5*(points[i_segment+1] - points[i_segment])
      up = piece.ups[i_segment] + piece.ups[i_segment+1]
      cylinder.render_to_center(
          vertex_buffer,
          points[i_segment],
          tangent,
          up,
          0.3,
          piece.residues[i_segment].color,
          piece.objids[i_segment])
      cylinder.render_to_center(
          vertex_buffer,
          points[i_segment+1],
          -tangent,
          up,
          0.3,
          piece.residues[i_segment+1].color,
          piece.objids[i_segment+1])
  return vertex_buffer


def make_ball_and_stick_mesh(rendered_soup):
  sphere = render.SphereShape(5, 5)
  cylinder = render.CylinderShape(5)
  n_vertex = len(rendered_soup.draw_to_screen_atoms)*sphere.n_vertex
  n_vertex += len(rendered_soup.bonds)*cylinder.n_vertex
  vertex_buffer = render.IndexedVertexBuffer(n_vertex)
  for atom in rendered_soup.draw_to_screen_atoms:
    point = atom.pos
    objid = atom.res_objid
    color = atom.residue.color
    sphere.render_to_center(
        vertex_buffer,
        point,
        point,
        point,
        0.2,
        color,
        objid)
  for bond in rendered_soup.bonds:
    color = bond.atom1.residue.color
    objid = bond.atom1.res_objid
    cylinder.render_to_center(
        vertex_buffer,
        bond.atom1.pos,
        bond.tangent,
        bond.up,
        0.2,
        color,
        objid)
  return vertex_buffer


class OpenGlHandler():
  def __init__(self, width, height):
    self.camera = camera.Camera()
    self.new_camera = camera.Camera()
    self.reshape(width, height)
    self.shader_catalog = shader.ShaderCatalog()
    self.shader = self.shader_catalog.shader
    self.draw_objects = []
    self.bg_rgb = (0.0, 0.0, 0.0, 0.0)
    self.is_outline = False

  def draw_objects_with_shader(self, shader_type):
    self.shader = self.shader_catalog.catalog[shader_type]
    gl.glUseProgram(self.shader.program)
    self.shader.bind_camera(self.camera)
    for draw_object in self.draw_objects:
      draw_object.draw(self.shader)

  def draw_to_screen(self):
    bg_color = self.camera.fog_color + [1.]
    gl.glClearColor(*bg_color)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT|gl.GL_DEPTH_BUFFER_BIT)
    gl.glEnable(gl.GL_BLEND)
    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glDepthFunc(gl.GL_LEQUAL)

    if self.is_outline:
      gl.glCullFace(gl.GL_BACK);
      gl.glEnable(gl.GL_CULL_FACE);
      self.draw_objects_with_shader('outline')

    gl.glCullFace(gl.GL_FRONT)
    gl.glEnable(gl.GL_CULL_FACE)
    self.draw_objects_with_shader('default')

  def pick(self, x, y):
    gl.glDisable(gl.GL_BLEND)
    gl.glClearColor(0.0, 0.0, 0.0, 0.0)
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
    gl.glCullFace(gl.GL_FRONT)
    gl.glEnable(gl.GL_CULL_FACE)
    self.draw_objects_with_shader('select')
    pixels = (c_float*4)()
    y_screen = self.height - y # screen and OpenGL y coord flipped
    gl.glReadPixels(x, y_screen, 1, 1, gl.GL_RGBA, gl.GL_FLOAT, pixels)
    objid = int(pixels[2]*255*256*256)
    objid += int(pixels[1]*255*256)
    objid += int(pixels[0]*255)
    return objid

  def reshape(self, width, height):
    gl.glViewport(0, 0, width, height)
    radius = .5 * min(width, height)
    self.camera.set_screen(width/radius, height/radius)
    self.width = width
    self.height = height



class PyBall:
  """
  Ties everything together.
  """

  def __init__(self, title, pdb):
    self.width = 640
    self.height = 480

    self.is_mouse_left_down = False
    self.is_mouse_right_down  = False

    self.x_mouse_save = 0.0
    self.y_mouse_save = 0.0

    self.init_glut(title)

    self.opengl = OpenGlHandler(self.width, self.height)

    self.soup = pdbatoms.Soup(pdb)
    self.rendered_soup = RenderedSoup(self.soup)

    self.opengl.camera.rescale(self.rendered_soup.scale)
    self.opengl.camera.set_center(self.rendered_soup.center)

    self.opengl.new_camera.rescale(self.rendered_soup.scale)

    self.n_step_animate = 0

    print "Building cylindrical trace..."
    self.ribbon_draw_object = make_cylinder_trace_mesh(self.rendered_soup, 6, 5, 10)
    print "Building ribbons..."
    self.ribbon_draw_object = make_carton_mesh(self.rendered_soup)
    print "Building arrows..."
    self.arrow_draw_object = make_arrow_mesh(self.rendered_soup)
    print "Building ball and sticks..."
    self.ball_stick_draw_object = make_ball_and_stick_mesh(self.rendered_soup)
    self.is_stick = False
    self.opengl.draw_objects.append(self.ribbon_draw_object)
    self.opengl.draw_objects.append(self.arrow_draw_object)

    self.set_callbacks()

    self.last = time.time()

  def init_glut(self, title):
    glut.glutInit()
    glut.glutInitWindowSize(self.width, self.height)
    glut.glutInitDisplayMode(glut.GLUT_RGBA|glut.GLUT_DOUBLE|glut.GLUT_DEPTH)
    glut.glutCreateWindow(title)
    
  def draw_to_screen(self):
    self.opengl.draw_to_screen()
    glut.glutSwapBuffers()

  def reshape(self, width, height):
    """window reshape callback."""
    self.opengl.reshape(width, height)
    self.opengl.draw_to_screen() 

  def get_window_dims(self):
    return glut.glutGet(glut.GLUT_WINDOW_WIDTH), glut.glutGet(glut.GLUT_WINDOW_HEIGHT)

  def keyboard(self, c, x=0, y=0):
    """keyboard callback."""
    if c == b'p':
      self.opengl.camera.is_perspective = not self.opengl.camera.is_perspective
      self.reshape(*self.get_window_dims())
    elif c == b'l':
      self.opengl.camera.is_lighting = not self.opengl.camera.is_lighting
    elif c == b'o':
      self.opengl.is_outline = not self.opengl.is_outline
    elif c == b'f':
      self.opengl.camera.is_fog = not self.opengl.camera.is_fog
    elif c == b's':
      self.is_stick = not self.is_stick
      self.opengl.draw_objects = []
      self.opengl.draw_objects.append(self.ribbon_draw_object)
      self.opengl.draw_objects.append(self.arrow_draw_object)
      if self.is_stick:
        self.opengl.draw_objects.append(self.ball_stick_draw_object)
    elif c.lower() == b'c':
      if glut.glutGetModifiers() == glut.GLUT_ACTIVE_CTRL:
        print "Exit"
        glut.glutLeaveMainLoop()
        sys.exit(1)
    elif c == b'q':
      sys.exit(0)
    glut.glutPostRedisplay()

  def mouse(self, button, state, x, y):
    if button == glut.GLUT_LEFT_BUTTON:
      self.objid = self.opengl.pick(x, y)
      if state == glut.GLUT_DOWN:
        self.save_objid = self.objid
      else:
        if self.save_objid == self.objid and self.objid > 0:
          atom = self.rendered_soup.objid_ref[self.objid]
          self.opengl.new_camera.center = atom.pos
          self.n_step_animate = 10
          print x, y, atom_name(atom)
      self.is_mouse_left_down = (state == glut.GLUT_DOWN)
      self.x_mouse_save, self.y_mouse_save = x, y
    elif button == glut.GLUT_RIGHT_BUTTON:
      self.is_mouse_right_down = (state == glut.GLUT_DOWN)
      self.x_mouse_save, self.y_mouse_save = x, y

  def get_polar(self, x_mouse, y_mouse):
    x = x_mouse - self.opengl.width/2
    y = y_mouse - self.opengl.height/2
    r = math.sqrt(x*x + y*y)
    if x != 0.0:
      theta = math.atan(y/float(x))
    else:
      if y > 0:
        theta = math.pi/2
      else:
        theta = -math.pi/2
    if x<0:
      if y>0:
        theta += math.pi
      else:
        theta -= math.pi
    return r, theta

  def motion(self, x_mouse, y_mouse):
    if self.is_mouse_left_down:
      scale = 0.5/self.opengl.camera.scale/math.sqrt(self.width**2 + self.height**2)
      x_diff = scale*(x_mouse - self.x_mouse_save)
      y_diff = scale*-(y_mouse - self.y_mouse_save)
      self.opengl.camera.rotate_xy(x_diff, y_diff)
    if self.is_mouse_right_down:
      r_save, theta_save = self.get_polar(self.x_mouse_save, self.y_mouse_save)
      r, theta = self.get_polar(x_mouse, y_mouse) 
      self.opengl.camera.rescale(exp(0.01*(r-r_save)))
      # self.opengl.camera.change_zoom(-0.1*(r-r_save))
      self.opengl.camera.rotate_z(theta - theta_save)
    self.x_mouse_save, self.y_mouse_save = x_mouse, y_mouse
    glut.glutPostRedisplay()

  def idle(self):
    now = time.time()
    elapsed = now - self.last
    time_step = 0.02
    if self.n_step_animate > 0:
      n_step = min(int(elapsed/time_step), self.n_step_animate)
      if n_step > 0:
        diff_center = self.opengl.new_camera.center - self.opengl.camera.center
        fraction = n_step/float(self.n_step_animate)
        self.n_step_animate -= n_step
        self.opengl.camera.center += fraction*diff_center
        self.draw_to_screen()
        self.last = now
    else:
      self.last = now

  def set_callbacks(self):
    glut.glutIdleFunc(self.idle)
    glut.glutReshapeFunc(self.reshape)
    glut.glutDisplayFunc(self.draw_to_screen)
    glut.glutKeyboardFunc(self.keyboard)
    glut.glutMouseFunc(self.mouse)
    glut.glutMotionFunc(self.motion)

  def run(self):
    return glut.glutMainLoop()


if __name__ == "__main__":
  pdb = '1ssx.pdb'
  if len(sys.argv) > 1:
    pdb = sys.argv[1]
  title = sys.argv[0].encode()
  pyball = PyBall(title, pdb)
  pyball.run()





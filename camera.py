# -*- coding: utf-8 -*-

"""
Camera object to hold all relevant variables to define
the camera of a 3D scene, and methods to convert this
to an OpenGL matrix for rendering.
"""


import pdbremix.v3numpy as v3
from ctypes import c_float



def matrix_to_c_floats(m, n_dim):
  n_float = n_dim*n_dim
  data = (c_float*n_float)()
  for i in range(n_dim):
    for j in range(n_dim):
      data[j*n_dim + i] = m[i,j]
  return data


class Camera:
  """
  This class holds all the variables to describe the
  camera for the viewer.
  """
  def __init__(self):
    self.rotation = v3.identity()
    self.center = v3.vector()
    self.scale = 1.0
    self.is_perspective = True
    self.is_lighting = True
    self.fog_near = -3
    self.fog_far = 20
    self.is_fog = True
    self.light_position = [100., 100., 500., 0.]
    self.half_width = 1.0
    self.half_height = 1.0

  def rotate_xy(self, rad_around_x, rad_around_y):
    rot_along_x = v3.rotation(v3.vector(0, 1, 0), rad_around_x)
    rot_along_y = v3.rotation(v3.vector(-1, 0, 0), rad_around_y)
    new_rotation = v3.combine(rot_along_x, rot_along_y)
    self.rotation = v3.combine(new_rotation, self.rotation)

  def rotate_z(self, rad_around_z):
    rot_along_z = v3.rotation(v3.vector(0, 0, -1), rad_around_z)
    self.rotation = v3.combine(rot_along_z, self.rotation)

  def rescale(self, new_scale):
    self.scale *= new_scale 
    
  def set_screen(self, width, height):
    self.half_width = width
    self.half_height = height

  def set_center(self, center):
    self.center = center

  def frustum(self):
    w, h = self.half_width, self.half_height
    l, r, b, t, n, f = -w, w, -h, h, 8, 16
    m = v3.identity()
    w, h, d = r-l, t-b, f-n
    m[0,:] = [2.*n/w, 0.,      (r+l)/w,  0.      ]
    m[1,:] = [0.,     2.*n/h,  (t+b)/h,  0.      ]
    m[2,:] = [0.,     0.,     -(f+n)/d, -2.*f*n/d]
    m[3,:] = [0.,     0.,     -1.,       0.      ]
    return m

  def ortho(self):
    w, h = self.half_width, self.half_height
    l, r, b, t, n, f = -w, w, -h, h, -2, 2 
    m = v3.identity()
    w, h, d = r-l, t-b, f-n
    m[0,:] = [2./w, 0.,    0.,   -(r+l)/w]
    m[1,:] = [0.,   2./h,  0.,   -(t+b)/h]
    m[2,:] = [0.,   0.,   -2./d, -(f+n)/d]
    m[3,:] = [0.,   0.,    0.,    1.     ]
    return m

  def modelview_cfloat16(self):
    modelview = v3.identity()

    modelview = v3.translation(-self.center)
    modelview = v3.combine(self.rotation, modelview)

    s = self.scale
    scaling_matrix = v3.scaling_matrix(s, s, s)
    modelview = v3.combine(scaling_matrix, modelview)

    if self.is_perspective:
      projection = v3.combine(
          self.frustum(), v3.translation(v3.vector(0, 0, -12)))
      projection = v3.combine(
          projection, v3.scaling_matrix(1.5, 1.5, 1.5))
    else:
      projection = ortho(-w, w, -h, h, -2, 2)
    modelview = v3.combine(projection, modelview)

    return matrix_to_c_floats(modelview, 4)

  def normal_cfloat9(self):
    # normal_matrix = m.transpose(m.inverse(modelview[:3,:3]))
    # since we only do uniform scaling, this reduces to
    return matrix_to_c_floats(self.rotation[:3,:3], 3)



    

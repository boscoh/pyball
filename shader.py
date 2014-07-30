# -*- coding: utf-8 -*-


"""
Shader objects
"""

from ctypes import sizeof, c_float, c_void_p, c_uint, string_at
from OpenGL.GLUT import *
from OpenGL.GL import *



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
  if (lighting) {
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



select_vertex_shader = """

uniform mat4 modelview_matrix;

attribute vec3 vertex;
attribute float objectid;

varying float objId;

void main(void) {
  gl_Position = modelview_matrix * vec4(vertex, 1.0);
  objId = objectid;
}
"""


select_fragment_shader = """

varying float objId;

int intMod(int x, int y) { 
  int z = x/y;
  return x - y*z;
}
void main(void) {
  // ints are only required to be 7bit...
  int integralObjId = int(objId + 0.5);
  int red = intMod(integralObjId, 256);
  integralObjId /= 256;
  int green = intMod(integralObjId, 256);
  integralObjId /= 256;
  int blue = intMod(integralObjId, 256);
  gl_FragColor = vec4(float(red), float(green), float(blue), 255.0)/255.0;
}
"""


class Shader():
  """
  Objects to hold compiled OpenGL shader_catalog, with the
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

  def set_matrix33(self, var_name, cfloat9):
    glUniformMatrix3fv(
        self.locations[var_name], 
        1, # number of matrices
        GL_FALSE, # transpose
        cfloat9)

  def set_matrix44(self, var_name, cfloat16):
    glUniformMatrix4fv(
        self.locations[var_name], 
        1, # number of matrices
        GL_FALSE, # transpose
        cfloat16)

  def set_vec4(self, var_name, c_float_4):
    glUniform4f(self.locations[var_name], *c_float_4)

  def bind_camera(self, camera):
    self.set_matrix44("modelview_matrix", camera.modelview_cfloat16())
    self.set_matrix33("normal_matrix", camera.normal_cfloat9())
    self.set_boolean('lighting', camera.is_lighting)
    self.set_vec4('light_position', camera.light_position)



class ShaderCatalog:
  def __init__(self):
    vertex_attribs = [
      "vertex", 
      "normal", 
      "color",
      "objectid",
    ]
    uniforms = [
      "modelview_matrix", 
      "normal_matrix", 
      "lighting", 
      "light_position", 
    ]
    self.catalog = {
      'default': 
        Shader(
          vertex_shader_source, 
          fragment_shader_source,
          vertex_attribs, 
          uniforms),
      'select': 
        Shader(
          select_vertex_shader, 
          select_fragment_shader,
          vertex_attribs, 
          uniforms)
    }
    self.shader = self.catalog['default']

# -*- coding: utf-8 -*-


"""
Shader objects
"""

from OpenGL import GL as gl


vertex_shader_source = b"""
uniform mat4 modelview_matrix;
uniform mat4 texture_matrix;
uniform mat3 normal_matrix;

uniform bool is_lighting;
uniform vec4 light_position;

attribute vec3 vertex;
attribute vec3 normal;
attribute vec3 color;

varying vec3 N, L, S;

void main() {
  gl_Position = modelview_matrix * vec4(vertex, 1.);
  
  if (is_lighting) {
    N = normalize(normal_matrix*normal);
    L = normalize(light_position.xyz);
    S = normalize(L+vec3(0, 0, 1));
  }
  gl_FrontColor = vec4(color, 1.);
}
"""

fragment_shader_source = b"""
const vec4 ambient_color = vec4(.2, .2, .2, 1.);
const vec4 diffuse_intensity = vec4(1., 1., 1., 1.); 

uniform bool is_lighting;
uniform bool is_fog;
uniform float fog_near;
uniform float fog_far;
uniform vec3 fog_color;

varying vec3 N, L, S;

void main() {
  vec4 color = gl_Color;
  if (is_lighting) {
    vec4 ambient = color * ambient_color;
    vec4 diffuse = color * diffuse_intensity;
    float d = max(0., dot(N, L));
    color = clamp(ambient + diffuse * d, 0., 1.);
  }

  gl_FragColor = color;

  if (is_fog) {
    float depth = gl_FragCoord.z / gl_FragCoord.w;
    float fog_factor = smoothstep(fog_near, fog_far, depth);
    gl_FragColor = mix(
        gl_FragColor, vec4(fog_color, gl_FragColor.w), fog_factor);
  }
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


# outline vertex shader. Expands vertices along the (in-screen) xy
# components of the normals.
outline_vertex_shader = """

attribute vec3 vertex;
attribute vec3 normal;
attribute vec3 color;
                                                                       
uniform mat4 modelview_matrix;
uniform mat3 normal_matrix;
uniform float outline_width;

varying vec3 N;

void main(void) {
  gl_Position = modelview_matrix * vec4(vertex, 1.0);
  N = normal_matrix * normalize(normal);
  gl_Position.xy += N.xy*outline_width;
}
"""

# outline shader. mixes outlineColor with fogColor
outline_fragment_shader = """

uniform vec3 outline_color;
uniform float fog_near;
uniform float fog_far;
uniform vec3 fog_color;
uniform bool is_fog;

void main() {
  gl_FragColor = vec4(outline_color, 1.0);
  float depth = gl_FragCoord.z / gl_FragCoord.w;
  if (is_fog) { 
    float fog_factor = smoothstep(fog_near, fog_far, depth);
    gl_FragColor = mix(
        gl_FragColor, vec4(fog_color, gl_FragColor.w), fog_factor);
  }
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

    self.program = gl.glCreateProgram()

    gl.glAttachShader(
       self.program,
       self.compile(gl.GL_VERTEX_SHADER, vertex_shader))
    gl.glAttachShader(
      self.program,
      self.compile(gl.GL_FRAGMENT_SHADER, fragment_shader))

    self.locations = {}

    # before linking, link names of vertex attribs
    for i, vertex_attrib in enumerate(vertex_attribs):
      self.locations[vertex_attrib] = i
      gl.glBindAttribLocation(self.program, i, vertex_attrib)

    gl.glLinkProgram(self.program)
    if gl.glGetProgramiv(self.program, gl.GL_LINK_STATUS) != gl.GL_TRUE:
      raise RuntimeError(gl.glGetProgramInfoLog(self.program))

    # after linking, link names of uniforms
    for uniform in uniforms:
       i = gl.glGetUniformLocation(self.program, uniform)
       self.locations[uniform] = i

  def compile(self, shader_type, source):
    shader = gl.glCreateShader(shader_type)
    gl.glShaderSource(shader, source)
    gl.glCompileShader(shader)
    if gl.glGetShaderiv(shader, gl.GL_COMPILE_STATUS) != gl.GL_TRUE:
      raise RuntimeError(gl.glGetShaderInfoLog(shader))
    return shader

  def set_boolean(self, var_name, is_state):
    if is_state:
      gl.glUniform1i(self.locations[var_name], 1)
    else:
      gl.glUniform1i(self.locations[var_name], 0)

  def set_matrix33(self, var_name, cfloat9):
    gl.glUniformMatrix3fv(
        self.locations[var_name], 
        1, # number of matrices
        gl.GL_FALSE, # transpose
        cfloat9)

  def set_matrix44(self, var_name, cfloat16):
    gl.glUniformMatrix4fv(
        self.locations[var_name], 
        1, # number of matrices
        gl.GL_FALSE, # transpose
        cfloat16)

  def set_vec4(self, var_name, c_float_4):
    gl.glUniform4f(self.locations[var_name], *c_float_4)

  def set_vec3(self, var_name, c_float_3):
    gl.glUniform3f(self.locations[var_name], *c_float_3)

  def set_float(self, var_name, c_float):
    gl.glUniform1f(
        self.locations[var_name],
        c_float)

  def bind_camera(self, camera):
    self.set_matrix44("modelview_matrix", camera.modelview_cfloat16())
    self.set_matrix33("normal_matrix", camera.normal_cfloat9())
    self.set_boolean('is_lighting', camera.is_lighting)
    self.set_boolean('is_fog', camera.is_fog)
    self.set_vec4('light_position', camera.light_position)
    self.set_float('fog_near', camera.fog_near)
    self.set_float('fog_far', camera.fog_far)
    self.set_vec3('fog_color', camera.fog_color)
    self.set_vec3('outline_color', camera.outline_color)
    self.set_float('outline_width', camera.outline_width)



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
      "is_lighting", 
      "light_position", 
      "fog_near",
      "fog_far",
      "is_fog",
      "fog_color",
      "outline_color",
      "outline_width",
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
          uniforms),
      'outline': 
        Shader(
          outline_vertex_shader, 
          outline_fragment_shader,
          vertex_attribs, 
          uniforms)
    }
    self.shader = self.catalog['default']

"""
Interactive 3D protein structure viewer.

Contains the main application window (MolecularViewerCanvas) with camera controls,
rendering modes, and debug console overlay.
"""

import OpenGL.GL as gl
from pdbstruct import parse
from vispy import app
from vispy.scene.visuals import Text

from vispy.gloo import Program

from .rendering import (
    Camera,
    make_ball_and_stick_mesh,
    make_calpha_arrow_mesh,
    make_carton_mesh,
    make_cylinder_trace_mesh,
    picking_fragment,
    picking_vertex,
    semilight_fragment,
    semilight_vertex,
)
from .structures import RenderedSoup


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

# -*- coding: utf-8 -*-

"""
PNG writer.

Copyright (c) 2010, Renaud Blanch <rndblnch at gmail dot com>
Licence: GPLv3 or higher <http://www.gnu.org/licenses/gpl.html>
"""


# imports ####################################################################

import ctypes
import struct
import zlib


# minimal png encoder ########################################################

def cat(gen):
	def joined(*args):
		return b"".join(iter(gen(*args)))
	return joined

@cat
def lines(width, height, depth, data):
	width *= depth
	for i in reversed(range(height)):
		yield struct.pack(">B%ss" % width, 0, data[i*width:(i+1)*width])

@cat
def chunk(chunk_type, data=b""):
	yield struct.pack(">I", len(data))
	yield chunk_type
	yield data
	yield struct.pack(">I", zlib.crc32(chunk_type + data) & 0xffffffff)

@cat
def image(width, height, depth, data):
	data = ctypes.string_at(data, width*height*depth)
	color = {3: 2, 4: 6}[depth]
	yield b"\x89PNG\r\n\x1a\n"
	yield chunk(b"IHDR", struct.pack(">2I5B", width, height, 8, color, 0, 0, 0))
	yield chunk(b"IDAT", zlib.compress(lines(width, height, depth, data)))
	yield chunk(b"IEND")

def write(output, width, height, depth, data):
	return output.write(image(width, height, depth, data))

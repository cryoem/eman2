#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine


# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from past.utils import old_div
from builtins import range
from builtins import object
from EMAN2 import *
import copy
from .emapplication import EMGLWidget
from .emdataitem3d import *
from .emglobjects import get_default_gl_colors
from .emitem3d import EMItem3D, EMItem3DInspector
from .empdbitem3d import *
from .emshapeitem3d import *
from libpyGLUtils2 import GLUtil
import math
import os
import pickle
from .valslider import ValSlider, EMLightControls, CameraControls, EMSpinWidget, EMQTColorWidget, EMANToolButton
import weakref

import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL import GLU
from OpenGL.GL import *
from PyQt5 import QtCore, QtGui, QtWidgets, QtOpenGL
from PyQt5.QtCore import Qt


#from emdataitem3d import EMDataItem3D, EMIsosurface, EMSliceItem3D, EMVolumeItem3D
# XPM format Cursors
visibleicon = [
    '16 12 3 1',
    'a c #0000ff',
    'b c #000000',
    'c c None',
    'cccccccccccccccc',
    'ccccccbbbbcccccc',
    'ccccbbbbbbbbcccc',
    'ccbbbccccccbbbcc',
    'cbbccccaaccccbbc',
    'cbccccaaaaccccbc',
    'cbccccaaaaccccbc',
    'cbbccccaaccccbbc',
    'ccbbbccccccbbbcc',
    'ccccbbbbbbbbcccc',
    'ccccccbbbbcccccc',
    'cccccccccccccccc'
]
    
invisibleicon = [
    '16 12 2 1',
    'b c #000000',
    'c c None',
    'cbbcccccccccbbcc',
    'ccbbcccccccbbccc',
    'cccbbcccccbbcccc',
    'ccccbbcccbbccccc',
    'cccccbbcbbcccccc',
    'ccccccbbbccccccc',
    'ccccccbbbccccccc',
    'cccccbbcbbcccccc',
    'ccccbbcccbbccccc',
    'cccbbcccccbbcccc',
    'ccbbcccccccbbccc',
    'cbbcccccccccbbcc'
]

zrotatecursor = [
    '15 14 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccc',
    'ccccbbbbbbccbcc',
    'ccbbbbbbbbbbbbc',
    'cbbbcccccbbbbbc',
    'bbbcccccbbbbbbb',
    'bbcccccbbbbbbbc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'cbbbbbbbcccccbb',
    'bbbbbbbcccccbbb',
    'cbbbbbcccccbbbc',
    'cbbbbbbbbbbbbcc',
    'ccbccbbbbbbcccc',
    'ccccccccccccccc'
]

xyrotatecursor = [
    '14 13 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccccccccc',
    'ccccbbbbbbcccc',
    'ccbbbbbbbbbbcc',
    'cbbbccccccbbbc',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'cbbbccccccbbbc',
    'ccbbbbbbbbbbcc',
    'ccccbbbbbbcccc',
    'cccccccccccccc'
]

crosshairscursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'ccbccccbbccccbccc',
    'cbbccccbbccccbbcc',
    'bbbbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbbbb',
    'cbbccccbbccccbbcc',
    'ccbccccbbccccbccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
]   

zhaircursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'cccbccccccccccccc',
    'ccbbbcccccccccccc',
    'cbcbcbccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccbbccccc',
    'cccbccccbbccccccc',
    'cccbccbbccccccccc',
    'cccbbbccccccccccc',
    'cccbccccccccccccc',
    'ccccbbccccccccccc',
    'ccccccbbccccccccc',
    'ccccccccbbccccccc',
    'ccccccccccbbccccc',
    'ccccccccccccccccc'
]

scalecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'bbbbbbbbcccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccbbbbbbccccc',
    'bccccbccccbccccc',
    'bccccbccccbccccc',
    'cccccbccccbccccb',
    'cccccbccccbccccb',
    'cccccbbbbbbccccb',
    'cccccccccccbcccb',
    'ccccccccccccbccb',
    'cccccccccccccbcb',
    'ccccccccccccccbb',
    'ccccccccbbbbbbbb'
]   

selectorcursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cbbbbbbbbbcccccc',
    'bcccccccccbccccc',
    'cbbbbbbbcccbcccc',
    'cccbccccccccbccc',
    'ccccbbbbccccbccc',
    'cccbccccccccbccc',
    'ccccbbbbcccbcbcc',
    'cccbccccccbcccbc',
    'ccccbbbbbbcccccb',
    'cccccccbccbcccbc',
    'ccccccccbccccbcc',
    'cccccccccbccbccc',
    'ccccccccccbbcccc',
    'cccccccccccccccc',
    'cccccccccccccccc',
    'cccccccccccccccc'
]

cubecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccccccccccccc',
    'cccccccccccccccccc',
    'ccccccbbbbbbbbbccc',
    'cccccbbbbbbbbbbccc',
    'ccccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbccccc',
    'ccbbbbbbbbbbcccccc',
    'ccbbbbbbbbbccccccc',
    'cccccccccccccccccc',
    'cccccccccccccccccc'
]

spherecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'cccccccbbbccccccc',
    'cccccbbbbbbbccccc',
    'ccccbbbbbbbbbcccc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccccbbbbbbbbbcccc',
    'cccccbbbbbbbccccc',
    'cccccccbbbccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

cylindercursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccbbbbbbbccccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

conecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccbbbbbbcccccc',
    'ccccbbbbbbbbccccc',
    'ccccbbbbbbbbccccc',
    'cccbbbbbbbbbbcccc',
    'cccbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbccc',
    'cbbbbbbbbbbbbbbcc',
    'cbbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbcccc',
    'cccccbbbbbbcccccc',
    'ccccccccccccccccc'
]

linecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccccccccccbbcc',
    'ccccccccccccbbccc',
    'cccccccccccbbcccc',
    'ccccccccccbbccccc',
    'cccccccccbbcccccc',
    'ccccccccbbccccccc',
    'cccccccbbcccccccc',
    'ccccccbbccccccccc',
    'cccccbbcccccccccc',
    'ccccbbccccccccccc',
    'cccbbcccccccccccc',
    'ccbbccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

datacursor = [
    '17 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbbcccbbb',
    'bbccccccbbbbcbbbb',
    'bbccccccbbbbbbbbb',
    'bbccccccbbcbbbcbb',
    'bbbbbcccbbcccccbb',
    'bbbbbcccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

textcursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbcccbbbcccbbcc',
    'ccbccccbbbccccbcc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'ccccccbbbbbcccccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

appcursor = [
    '17 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'caaacccaaacccaaac',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'aaaaacaaaaccaaaac',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

cubeicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccccccccccccc',
    'cccccccccccccccccc',
    'ccccccbbbbbbbbbccc',
    'cccccbbbbbbbbbbccc',
    'ccccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbccccc',
    'ccbbbbbbbbbbcccccc',
    'ccbbbbbbbbbccccccc',
    'cccccccccccccccccc',
    'cccccccccccccccccc'
] 

sphereicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'cccccccbbbccccccc',
    'cccccbbbbbbbccccc',
    'ccccbbbbbbbbbcccc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccccbbbbbbbbbcccc',
    'cccccbbbbbbbccccc',
    'cccccccbbbccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

cylindericon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccbbbbbbbccccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

coneicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccbbbbbbcccccc',
    'ccccbbbbbbbbccccc',
    'ccccbbbbbbbbccccc',
    'cccbbbbbbbbbbcccc',
    'cccbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbccc',
    'cbbbbbbbbbbbbbbcc',
    'cbbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbcccc',
    'cccccbbbbbbcccccc',
    'ccccccccccccccccc'
]

lineicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccccccccccbbcc',
    'ccccccccccccbbccc',
    'cccccccccccbbcccc',
    'ccccccccccbbccccc',
    'cccccccccbbcccccc',
    'ccccccccbbccccccc',
    'cccccccbbcccccccc',
    'ccccccbbccccccccc',
    'cccccbbcccccccccc',
    'ccccbbccccccccccc',
    'cccbbcccccccccccc',
    'ccbbccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

texticon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbcccbbbcccbbcc',
    'ccbccccbbbccccbcc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'ccccccbbbbbcccccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

dataicon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbbcccbbb',
    'bbccccccbbbbcbbbb',
    'bbccccccbbbbbbbbb',
    'bbccccccbbcbbbcbb',
    'bbbbbcccbbcccccbb',
    'bbbbbcccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

appicon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'caaacccaaacccaaac',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'aaaaacaaaaccaaaac',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

rotateicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccbbbbbbccbcc',
    'ccbbbbbbbbbbbbc',
    'cbbbcccccbbbbbc',
    'bbbcccccbbbbbbb',
    'bbcccccbbbbbbbc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'cbbbbbbbcccccbb',
    'bbbbbbbcccccbbb',
    'cbbbbbcccccbbbc',
    'cbbbbbbbbbbbbcc',
    'ccbccbbbbbbcccc',
    'ccccccccccccccc'
]

crosshairsicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'ccbccccbbccccbccc',
    'cbbccccbbccccbbcc',
    'bbbbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbbbb',
    'cbbccccbbccccbbcc',
    'ccbccccbbccccbccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
] 

multiselectoricon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cbcbcbcbcbcbcbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'cccccccccccccbbbc',
    'cbcbcbcbcbcbcbbbc',
    'cccccccccccccbbbb',
    'ccccccccccccccccb',
    'ccccccccccccccccc'
] 

scaleicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'bbbbbbbbcccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccbbbbbbccccc',
    'bccccbccccbccccc',
    'bccccbccccbccccc',
    'cccccbccccbccccb',
    'cccccbccccbccccb',
    'cccccbbbbbbccccb',
    'cccccccccccbcccb',
    'ccccccccccccbccb',
    'cccccccccccccbcb',
    'ccccccccccccccbb',
    'ccccccccbbbbbbbb'
]

ztransicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'cccbccccccccccccc',
    'ccbbbcccccccccccc',
    'cbcbcbccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccbbccccc',
    'cccbccccbbccccccc',
    'cccbccbbccccccccc',
    'cccbbbccccccccccc',
    'cccbccccccccccccc',
    'ccccbbccccccccccc',
    'ccccccbbccccccccc',
    'ccccccccbbccccccc',
    'ccccccccccbbccccc',
    'ccccccccccccccccc'
] 

selectionicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccbcccccccccccccc',
    'ccbbccccccccccccc',
    'ccbcbbccccccccccc',
    'ccbcccbcccccccccc',
    'ccbccccbbcccccccc',
    'ccbccccccbccccccc',
    'ccbcccccccbbccccc',
    'ccbcccccbbbbbcccc',
    'ccbccbccbcccccccc',
    'ccbcbcbccbccccccc',
    'ccbbccbccbccccccc',
    'ccbccccbccbcccccc',
    'ccccccccbbbcccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

rulericon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccbbbbbbbbbbbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbbbbbbbccc',
    'ccccccccccccccccc'
] 

class EMScene3D(EMItem3D, EMGLWidget):
	"""
	Widget for rendering 3D objects. Uses a scene graph for rendering
	"""
	sgmousepress = QtCore.pyqtSignal(float, float)
	sgmousemove = QtCore.pyqtSignal(float, float)
	sgmouserelease = QtCore.pyqtSignal(float, float)
	name = "SG"
	def __init__(self, parentwidget=None, SGactivenodeset=set(), scalestep=0.5):
		"""
		@param parent: The parent of the widget
		@param SGnodelist: a list enumerating all the SGnodes
		@param SGactivenodeset: a set enumerating the list of active nodes
		@param scalestep: The step to increment the object scaling
		"""
		EMItem3D.__init__(self, parent=None, transform=Transform())
		EMGLWidget.__init__(self,parentwidget)
		QtOpenGL.QGLFormat().setDoubleBuffer(True)
		QtOpenGL.QGLFormat().setDepth(True)
		self.setSelectedItem(True)			# The root is selected by default
		self.currentselecteditem = self
		self.camera = EMCamera(1.0, 500.0)		# Default near,far
		self.clearcolor = [0.0, 0.0, 0.0, 0.0]		# Back ground color	
		self.main_3d_inspector = None			# No inspector by default
		self.apix = None				# No angstrom per pixel to begin with
		self.item_inspector = None			# Get the inspector GUI
		self.reset_camera = False			# Toggle flag to determine if the clipping plane has changed and needs redrawing
		self.zslicemode = False				# Set to true when making a Z slice
		#self.SGactivenodeset = SGactivenodeset		# A set of all active nodes (currently not used)
		self.scalestep = scalestep			# The scale factor stepsize
		self.toggle_render_selectedarea = False		# Don't render the selection box by default
		self.mousemode = "rotate"			# The mouse mode
		self.initialviewportdims = None			# Default is not to set viewport dims, the objects will appear as actual sizes
		self.zrotatecursor = QtGui.QCursor(QtGui.QPixmap(zrotatecursor),-1,-1)
		self.xyrotatecursor = QtGui.QCursor(QtGui.QPixmap(xyrotatecursor),-1,-1)
		self.crosshaircursor = QtGui.QCursor(QtGui.QPixmap(crosshairscursor),-1,-1)
		self.scalecursor = QtGui.QCursor(QtGui.QPixmap(scalecursor),-1,-1)
		self.zhaircursor = QtGui.QCursor(QtGui.QPixmap(zhaircursor),-1,-1)
		self.selectorcursor = QtGui.QCursor(QtGui.QPixmap(selectorcursor),-1,-1)
		self.linecursor = QtGui.QCursor(QtGui.QPixmap(linecursor),-1,-1)
		self.cubecursor = QtGui.QCursor(QtGui.QPixmap(cubecursor),-1,-1)
		self.spherecursor = QtGui.QCursor(QtGui.QPixmap(spherecursor),-1,-1)
		self.cylindercursor = QtGui.QCursor(QtGui.QPixmap(cylindercursor),-1,-1)
		self.conecursor = QtGui.QCursor(QtGui.QPixmap(conecursor),-1,-1)
		self.datacursor = QtGui.QCursor(QtGui.QPixmap(datacursor),-1,-1)
		self.textcursor = QtGui.QCursor(QtGui.QPixmap(textcursor),-1,-1)
		self.appcursor = QtGui.QCursor(QtGui.QPixmap(appcursor),-1,-1)
		
	def getEvalString(self):
		"""
		Return a string that after eval can reinstantiate the object
		"""
		return "SG"
		
	def initializeGL(self):
		glClearColor(self.clearcolor[0], self.clearcolor[1], self.clearcolor[2], self.clearcolor[3])		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_NORMALIZE)
		self.firstlight = EMLight(GL_LIGHT0)
		self.firstlight.enableLighting()
        
	def paintGL(self):
		#if self.reset_camera: self.camera.update()
		self.camera.update()
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		#Call rendering
		self.renderSelectedArea() 	# Draw the selection box if needed
		self.render()			# SG nodes must have a render method
		glFlush()			# Finish rendering
		self.reset_camera = False
		
	def resizeGL(self, width, height):
		self.camera.update(width, height)
		self.updateInspector()
	
	def renderNode(self):
		pass
	
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""	
		if not self.item_inspector: self.item_inspector = EMSGNodeInspector("All Objects", self)
		return self.item_inspector
		
	def setInspector(self, inspector):
		"""
		Set the main 3d inspector
		"""
		self.main_3d_inspector = inspector
	
	def showInspector(self):
		"""
		Show the SG inspector
		"""
		if not self.main_3d_inspector:
			self.main_3d_inspector = EMInspector3D(self)
			self.main_3d_inspector.updateTree()
			self.main_3d_inspector.updateInspector()
			self.main_3d_inspector.show()
			
			
	def pickItem(self):
		"""
		Pick an item on the screen using openGL's selection mechanism
		"""
		self.makeCurrent()	# This is needed so we use openGL for the correct widget!!! An not say the lights or camera 3D widget
		viewport = glGetIntegerv(GL_VIEWPORT)
		glSelectBuffer(1024)	# The buffer size for selection
		glRenderMode(GL_SELECT)
		glInitNames()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		
		# Find the selection box. Go from Volume view coords to viewport coords. sa = selection area
		dx = self.sa_xf - self.sa_xi
		dy = self.sa_yf - self.sa_yi
		x = (self.sa_xi + old_div(viewport[2],2)) 
		y = (old_div(viewport[3],2) - self.sa_yi)
		if dx < 2 and dx > -2: dx = 2
		if dy < 2 and dy > -2: dy = 2

		# Apply selection box, and center it in the green box
		GLU.gluPickMatrix(x + old_div(dx,2), (viewport[3]- 2*self.camera.getPseudoFovyHeight() - y) + old_div(dy,2), int(math.fabs(dx)), int(math.fabs(dy)), viewport)
		self.camera.setProjectionMatrix()
		
		#drawstuff, but first we need to remove the influence of any previous xforms which ^$#*$ the selection
		glMatrixMode(GL_MODELVIEW)
		glPushMatrix()
		glLoadIdentity()
		self.camera.setCameraPosition(sfactor=1) # Factor of two to compensate for the camera already being set
		self.render()
		glPopMatrix()
		
		# Return to default state
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		records = glRenderMode(GL_RENDER)
		
		# process records
		self.processSelection(records)
	
	def selectArea(self, xi, xf, yi, yf, togglearea=True):
		"""
		Set an area for selection. Need to switch between viewport coords, where (0,0 is bottom left) to
		volume view coords where 
0,0) is center of the screen.
		"""
		self.sa_xi = xi - old_div(self.camera.getWidth(),2) - self.camera.getPseudoFovyWidth()
		self.sa_xf = xf - old_div(self.camera.getWidth(),2) - self.camera.getPseudoFovyWidth()
		self.sa_yi = -yi + old_div(self.camera.getHeight(),2) + self.camera.getPseudoFovyHeight()
		self.sa_yf = -yf + old_div(self.camera.getHeight(),2) + self.camera.getPseudoFovyHeight()
		self.toggle_render_selectedarea = togglearea
		
	def deselectArea(self):
		"""
		Turn off selection box
		"""
		self.sa_xi = 0.0
		self.sa_xf = 0.0
		self.sa_yi = 0.0
		self.sa_yf = 0.0
		self.toggle_render_selectedarea = False
		
	def renderSelectedArea(self):
		"""
		Draw the selection box, box is always drawn orthographically
		"""
		if self.toggle_render_selectedarea: 
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			glMatrixMode(GL_PROJECTION)
			glPushMatrix()
			glLoadIdentity()
			self.camera.setOrthoProjectionMatrix()
			glColor3f(0.0,1.0,0.0)
			glMaterialfv(GL_FRONT, GL_AMBIENT, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_DIFFUSE, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_SPECULAR, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,1.0,0.0,1.0])
			glBegin(GL_LINE_LOOP)
			# set the box just in front of the clipping plane
			z = -self.camera.getClipNear() - self.camera.getZclip()
			sfw = self.camera.getViewPortWidthScaling()
			sfh = self.camera.getViewPortHeightScaling()
			glVertex3f(sfw*(self.sa_xi+self.camera.getPseudoFovyWidth()), sfh*(self.sa_yi-self.camera.getPseudoFovyHeight()), z)
			glVertex3f(sfw*(self.sa_xi+self.camera.getPseudoFovyWidth()), sfh*(self.sa_yf-self.camera.getPseudoFovyHeight()), z)
			glVertex3f(sfw*(self.sa_xf+self.camera.getPseudoFovyWidth()), sfh*(self.sa_yf-self.camera.getPseudoFovyHeight()), z)
			glVertex3f(sfw*(self.sa_xf+self.camera.getPseudoFovyWidth()), sfh*(self.sa_yi-self.camera.getPseudoFovyHeight()), z)
			glEnd()
			glPopMatrix()
			glMatrixMode(GL_MODELVIEW)
			glPopAttrib()
	
	def _IsAncestorSelected(self, item):
		"""
		If an Ancestor is selected don't deselect by selecting a descendant
		"""
		if item.getSelectedAncestorNodes():
			return True
		else:
			return False
		
	def processSelection(self, records):
		"""
		Process the selection records
		"""
		# Remove old selection if not in append mode
		if (not records and not self.toggleselection) or (self.multiselect and not self.appendselection):
			self.clearSelection()
			
		# Select the desired items	
		closestitem = None
		bestdistance = 1.0
		for record in records:
			if self.multiselect:
				selecteditem = EMItem3D.selection_idx_dict[record.names[len(record.names)-1]]()
				# If an ancestor is selected consider THIS item to be selected, WHY????
				#if self._IsAncestorSelected(selecteditem): return
				# Select the data itself and not the isosurface, slice, etc 
				if selecteditem.nodetype == "DataChild": selecteditem = selecteditem.parent
				selecteditem.setSelectedItem(True)
				# Inspector tree management
				self.updateTreeSelVis(selecteditem)
				self.setCurrentSelection(selecteditem)
			else:
				if record.near < bestdistance:
					bestdistance = record.near
					closestitem = record
		if closestitem:
			selecteditem = EMItem3D.selection_idx_dict[closestitem.names[len(closestitem.names)-1]]()
			# If an ancestor is selected consider THIS item to be selected, WHY?????
			#if self._IsAncestorSelected(selecteditem): return
			# Select the data itself and no the isosurface, slice, etc
			if selecteditem.nodetype == "DataChild": selecteditem = selecteditem.parent
			if not selecteditem.isSelectedItem() and not self.appendselection and not self.toggleselection: self.clearSelection()
			if self.toggleselection and selecteditem.isSelectedItem():
				selecteditem.setSelectedItem(False)
			else:
				selecteditem.setSelectedItem(True)
			# Inspector tree management
			self.updateTreeSelVis(selecteditem)
			self.setCurrentSelection(selecteditem)
	
	def clearSelection(self):
		"""
		Clear all selected items to deselect
		"""
		for selected in self.getAllSelectedNodes():
			selected.setSelectedItem(False)
			# Inspector tree management
			self.updateTreeSelVis()

	# Event subclassing
	def mousePressEvent(self, event):
		"""
		QT event handler. Records the coords when a mouse button is pressed and sets the cursor depending on what button(s) are pressed
		"""
		# The previous x,y records where the mouse was prior to mouse move, rest upon mouse move
		self.previous_x = event.x()
		self.previous_y = event.y()
		# The first x,y records where the mouse was first pressed
		self.first_x = self.previous_x
		self.first_y = self.previous_y
		# Process mouse events
		if (event.buttons()&Qt.LeftButton and self.mousemode == "app"):
			QtWidgets.qApp.setOverrideCursor(self.appcursor)
			self.sgmousepress.emit(event.x(), event.y())
		if (event.buttons()&Qt.LeftButton and self.mousemode == "data"):
			QtWidgets.qApp.setOverrideCursor(self.datacursor)
			filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Get file', os.getcwd())[0]
			if not filename: return
			name = os.path.basename(str(filename))
			self.newnode = EMDataItem3D(filename, transform=self._gettransformbasedonscreen(event, rescale=False))
			self.insertNewNode(name, self.newnode)
			self.newnode.setTransform(self.newnode.getParentMatrixProduct().inverse()*self.newnode.getTransform())
			self.isonode = EMIsosurface(self.newnode, transform=Transform())
			self.clearSelection()
			self.isonode.setSelectedItem(True)
			self.insertNewNode("Isosurface", self.isonode, parentnode=self.newnode)
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "text"):
			QtWidgets.qApp.setOverrideCursor(self.textcursor)
			text, ok = QtWidgets.QInputDialog.getText(self, 'Enter Text', '')
			if ok:
				self.newnode = EM3DText(str(text), 32.0, transform=self._gettransformbasedonscreen(event))
				self._insert_shape(text, self.newnode)
				self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "line"):
			QtWidgets.qApp.setOverrideCursor(self.linecursor)
			self.newnode = EMLine(0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 20.0, transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Line", self.newnode)
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cube"):
			QtWidgets.qApp.setOverrideCursor(self.cubecursor)
			self.newnode = EMCube(2.0, transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Cube", self.newnode)
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "sphere"):
			QtWidgets.qApp.setOverrideCursor(self.spherecursor)
			self.newnode = EMSphere(2.0, transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Sphere", self.newnode)
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cylinder"):
			QtWidgets.qApp.setOverrideCursor(self.cylindercursor)
			self.newnode = EMCylinder(2.0,2.0, transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Cylinder", self.newnode)
			self.newnode.updateMatrices([90,1,0,0], "rotate")
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cone"):
			QtWidgets.qApp.setOverrideCursor(self.conecursor)
			self.newnode = EMCone(2.0,2.0, transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Cone", self.newnode)
			self.newnode.updateMatrices([90,1,0,0], "rotate")
			self.updateSG()	
		if (event.buttons()&Qt.LeftButton and self.mousemode == "rotate"):
			if  event.y() > 0.95*self.size().height(): # The lowest 5% of the screen is reserved from the Z spin virtual slider
				QtWidgets.qApp.setOverrideCursor(self.zrotatecursor)
				self.zrotate = True
			else:
				QtWidgets.qApp.setOverrideCursor(self.xyrotatecursor)
				self.zrotate = False
		if (event.buttons()&Qt.LeftButton and self.mousemode == "ruler"):
			self.newnode = EMRuler(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.getAPix(), self.camera.getViewPortWidthScaling(), transform=self._gettransformbasedonscreen(event))
			self._insert_shape("Ruler", self.newnode, clearsel=False)
			#self.newnode.updateMatrices([90,1,0,0], "rotate")
			self.updateSG()
		if event.buttons()&Qt.LeftButton and self.mousemode == "scale":
			QtWidgets.qApp.setOverrideCursor(self.scalecursor)
		if (event.buttons()&Qt.LeftButton and self.mousemode == "selection"): 
			#QtWidgets.qApp.setOverrideCursor(self.selectorcursor)
			self.multiselect = False
			self.appendselection = False
			self.toggleselection = False
			if event.modifiers()&Qt.ShiftModifier:
				self.appendselection = True
			if event.modifiers()&Qt.ControlModifier:
				self.toggleselection = True
			# 5 seems a big enough selection box
			self.selectArea(event.x(), event.x() + 5, event.y(), event.y() + 5, togglearea=False)
			self.pickItem()
			self.updateSG()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "multiselection"):
			#QtWidgets.qApp.setOverrideCursor(self.selectorcursor)
			self.multiselect = True
			self.appendselection = False
			self.toggleselection = False
			self.selectArea((event.x()-2.0), event.x(), (event.y()-2.0), event.y()) # To enable selection just by clicking
			if event.modifiers()&Qt.ShiftModifier:
				self.appendselection = True
		if (event.buttons()&Qt.LeftButton and self.mousemode == "ztranslate"):
			QtWidgets.qApp.setOverrideCursor(self.zhaircursor)
		if event.buttons()&Qt.RightButton or (event.buttons()&Qt.LeftButton and self.mousemode == "xytranslate"):
			QtWidgets.qApp.setOverrideCursor(self.crosshaircursor)
		if event.buttons()&Qt.MidButton or (event.buttons()&Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.showInspector()
			
		if (event.buttons()&Qt.LeftButton and event.modifiers()&Qt.ControlModifier and event.modifiers()&Qt.ShiftModifier):
			### so one can open it from a linux remote controlled by a mac...
			self.showInspector()
	
	def _insert_shape(self, name, node, clearsel=True):
		""" Helper function for mousePressEvent to reduce code duplication"""
		if clearsel: self.clearSelection()
		node.setSelectedItem(True)
		node.setHiddenSelected(self.getHiddenSelected()) # inherit hidden state from SG
		self.insertNewNode(name, node)
		node.setTransform(self.newnode.getParentMatrixProduct().inverse()*self.newnode.getTransform()) # so the object is not modified by parent transform upon insertion
		
	def _gettransformbasedonscreen(self, event, rescale=True):
		""" Helper function to for mousePressEvent"""
		x = self.camera.getViewPortWidthScaling()*(event.x() - old_div(self.camera.getWidth(),2))
		y = self.camera.getViewPortHeightScaling()*(-event.y() + old_div(self.camera.getHeight(),2))
		
		if rescale:
			return Transform({"type":"eman","tx":x,"ty":y,"scale":self.camera.getViewPortWidthScaling()})
		else:
			return Transform({"type":"eman","tx":x,"ty":y})
		
	def mouseMoveEvent(self, event):
		"""
		Qt event handler. Scales the SG depending on what mouse button(s) are pressed when dragged
		"""
		dx = (event.x() - self.previous_x)*self.camera.getViewPortWidthScaling()
		dy = (event.y() - self.previous_y)*self.camera.getViewPortHeightScaling()
		x = event.x()
		y = event.y()
		if (event.buttons()&Qt.LeftButton and self.mousemode == "app"):
			self.sgmousemove.emit(event.x(), event.y())
		if (event.buttons()&Qt.LeftButton and self.mousemode == "line"):
			self.newnode.setEndAndWidth(0.0, 0.0, 0.0, x - self.first_x, self.first_y - y, 0.0, 20.0)
		if (event.buttons()&Qt.LeftButton and self.mousemode == "ruler"):
			self.newnode.setRuler(0.0, 0.0, 0.0, x - self.first_x, self.first_y - y, 0.0)
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cube"):
			self.newnode.setSize(math.sqrt((x - self.first_x)**2 + (y - self.first_y)**2))
		if (event.buttons()&Qt.LeftButton and self.mousemode == "sphere"):
			self.newnode.setRadius(math.sqrt((x - self.first_x)**2 + (y - self.first_y)**2))
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cylinder"):
			self.newnode.setRadiusAndHeight(math.fabs(x - self.first_x), math.fabs(y - self.first_y))
		if (event.buttons()&Qt.LeftButton and self.mousemode == "cone"):
			self.newnode.setRadiusAndHeight(math.fabs(x - self.first_x), math.fabs(y - self.first_y))
		if (event.buttons()&Qt.LeftButton and self.mousemode == "rotate"):
			magnitude = math.sqrt(dx*dx + dy*dy)
			# We want to remove the effect of self.camera.getViewPortWidthScaling() for rotation. For everything else the effect is desired
			#Check to see if the cursor is in the 'virtual slider panel'
			try:
				if  self.zrotate: # The lowest 5% of the screen is reserved from the Z spin virtual slider
					self.updateMatrices([old_div(dx,self.camera.getViewPortWidthScaling()),0,0,-1], "rotate")
				else:
					self.updateMatrices([magnitude/self.camera.getViewPortWidthScaling(),-dy/magnitude,-dx/magnitude,0], "rotate")
			except ValueError: pass
			except ZeroDivisionError: pass
		if (event.buttons()&Qt.LeftButton and self.mousemode == "selection") and not event.modifiers()&Qt.ControlModifier:
			self.updateMatrices([dx,-dy,0], "translate")
		if (event.buttons()&Qt.LeftButton and self.mousemode == "multiselection"):
			self.selectArea(self.first_x, event.x(), self.first_y, event.y())
		if (event.buttons()&Qt.LeftButton and self.mousemode == "ztranslate"):
			self.updateMatrices([0,0,(-dy)], "translate")
		if event.buttons()&Qt.RightButton or (event.buttons()&Qt.LeftButton and self.mousemode == "xytranslate"):
			self.updateMatrices([dx,-dy,0], "translate")
		if event.buttons()&Qt.LeftButton and self.mousemode == "scale":
			self.updateMatrices([self.scalestep*0.1*(dx+dy)], "scale")
		self.previous_x =  x
		self.previous_y =  y
		self.updateSG()	
			
	def mouseReleaseEvent(self, event):
		"""
		Qt event handler. Returns the cursor to arrow upon mouse button release
		"""
		if (event.buttons()&Qt.LeftButton and self.mousemode == "app"):
			self.sgmouserelease.emit([event.x(), event.y()])
			
		QtWidgets.qApp.setOverrideCursor(Qt.ArrowCursor)
		# Select using the selection box
		if self.toggle_render_selectedarea:
			self.pickItem()
			self.deselectArea()
			self.updateSG()
			
	def wheelEvent(self, event):
		"""
		QT event handler. Scales the SG upon wheel movement, does so by changing fovy or orthographic equivalent
		"""
		# Originally the wheel scaled by zoom the viewport, but that caused all sorts of issues, so I now just scale the SG
		# The 25 is a fudge factor that controls the speed of scaling, lower if slower scaling
		if event.angleDelta().y() > 0:
			if self.camera.getUseOrtho():
				self.camera.setPseudoFovy(self.camera.getPseudoFovyWidth()+old_div((self.camera.getPseudoFovyWidth()+self.camera.getWidth()),25))
			else:
				self.camera.setFovy(self.camera.getFovy()+1.0)
		else:
			if self.camera.getUseOrtho():
				self.camera.setPseudoFovy(self.camera.getPseudoFovyWidth()-old_div((self.camera.getPseudoFovyWidth()+self.camera.getWidth()),25))
			else:
				self.camera.setFovy(self.camera.getFovy()-1.0)
		self.updateSG()
		self.cameraNeedsanUpdate()
			
			
	def mouseDoubleClickEvent(self,event):
		print("Mouse Double Click Event")
	
	# Grab and release the keyboard
	def enterEvent(self, event):
		self.grabKeyboard()	
		
	def leaveEvent(self, event):
		self.releaseKeyboard()
		
	def keyPressEvent(self,event):
		"""Process keyboard stuff"""
		# Delete nodes
		if event.key() == QtCore.Qt.Key_Delete:
			for node in self.getAllSelectedNodes():
				self.deleteNode(node)
			self.setCurrentSelection(self)
			self.updateSG()
		# Hide selections
		if event.key() == QtCore.Qt.Key_H:
			for node in self.getAllNodes():
				node.setHiddenSelected(not node.isSelectionHidded())
			
			self.updateSG()
		# Now toolbar hot keys
		if event.key() == QtCore.Qt.Key_Escape:
			self.setMouseMode("selection")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_M:
			self.setMouseMode("multiselection")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_T:
			self.setMouseMode("xytranslate")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_Z:
			self.setMouseMode("ztranslate")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_R:
			self.setMouseMode("rotate")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_S:
			self.setMouseMode("scale")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_L:
			self.setMouseMode("line")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_U and not event.modifiers()&Qt.ControlModifier:
			self.setMouseMode("ruler")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_C and not event.modifiers()&Qt.ControlModifier:
			self.setMouseMode("cube")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_P:
			self.setMouseMode("sphere")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_Y:
			self.setMouseMode("cylinder")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_O:
			self.setMouseMode("cone")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_X:
			self.setMouseMode("text")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_E:
			self.setMouseMode("data")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_A:
			self.setMouseMode("app")
			self.updateInspector()
		if event.key() == QtCore.Qt.Key_C and event.modifiers()&Qt.ControlModifier:
			nodes = self.getAllSelectedNodes()
			if not nodes: return
			self.copyNodes(nodes)
		if event.key() == QtCore.Qt.Key_G and event.modifiers()&Qt.ControlModifier:
			nodes = self.getAllSelectedNodes()
			if not nodes: return
			self.groupNodes(nodes)
		if event.key() == QtCore.Qt.Key_U and event.modifiers()&Qt.ControlModifier:
			nodes = self.getAllSelectedNodes()
			if not nodes: return
			self.unGroupNodes(nodes)
	
	def deleteNode(self, node):
		"""
		Remove a node from the SG
		@param node, the EMItem3D to remove
		"""
		if node.parent:
			node.parent.removeChild(node)
			self.updateTree()	# Update the inspector tree
			
	def copyNodes(self, nodes, parentnode=None, parentidx=None):
		"""
		Copy a set nof node(s), parentnode and parentidx are only used for recursion
		"""
		for node in nodes:
			if node.name == "Data" or node.name == "SG" or node.nodetype == "DataChild": continue	# Can't copy data nodes
			copiednode = eval(node.getEvalString())
			copiednode.setTransform(copy.deepcopy(node.getTransform()))
			name = "Copy"
			if node.getLabel(): name = node.getLabel()+name
			try:
				copiednode.setAmbientColor(node.ambient[0], node.ambient[1], node.ambient[2])
				copiednode.setSpecularColor(node.specular[0], node.specular[1], node.specular[2])
				copiednode.setDiffuseColor(node.diffuse[0], node.diffuse[1], node.diffuse[2])
			except:
				pass
			
			self.insertNewNode(name, copiednode, parentnode=parentnode, parentidx=parentidx)
			if node.getChildren(): 
				self.copyNodes(node.getChildren(), parentnode=copiednode, parentidx=0)
	
	def groupNodes(self, nodes):
		"""
		Group nodes into a new base node
		"""
		skipped = 0
		newbasenode = EMItem3D(transform=Transform())
		self.insertNewNode("Node", newbasenode)
		for i,node in enumerate(nodes):
			if node.name == "SG": 
				skipped +=1
				continue	# Don't group these 
			node.parent.getChildren().remove(node)
			newbasenode.addChild(node)
			self.updateTree(self.getCurrentSelection())	# Update the inspector tree
	
	def unGroupNodes(self, nodes):
		"""
		UnGroup nodes
		"""
		skipped = 0
		parentnode = None
		for i,node in enumerate(nodes):
			parentnode = node.parent
			if parentnode.name == "SG": 
				skipped +=1
				continue
			parentnode.getChildren().remove(node)
			parentnode.parent.addChild(node)
		
		if not parentnode.getChildren():
			parentnode.parent.getChildren().remove(parentnode)
		self.updateTree(self.getCurrentSelection())	# Update the inspector tree
			
	def setMouseMode(self, mousemode):
		"""
		Sets the mouse mode, used by the inspector
		"""
		self.mousemode = mousemode
		
	def getMouseMode(self):
		"""
		Return the mouse mode
		"""
		return self.mousemode
		
	def setZclip(self, zclip):
		""" Set the Z clipping plane """
		self.camera.setZclip(zclip)
		self.reset_camera = True
		
	def setZslice(self):
		"""
		Get a Z slice to display in the camera widget, only works for orthographic mode
		"""
		# Getting the Z slice will have problems when using perspective viewing
		self.setAutoBufferSwap(False)
		oldnear = self.camera.getClipNear() 
		oldfar = self.camera.getClipFar()
		self.zslicemode = True
		# We want to see the full volume data rather than a clip, so move clipping planes to BIG
		# BIG is a bit different for perspective and orthographic volumes
		if self.camera.usingortho:
			self.camera.setClipNear(-1000)
			self.camera.setClipFar(1000)
		else:
			self.camera.setClipNear(1)
			self.camera.setClipFar(1000)
		self.reset_camera = True
		current_xform = self.getTransform()
		current_xform_side = Transform({"type":"spin","omega":90,"n1":0,"n2":1,"n3":0})*current_xform
		self.setTransform(current_xform_side)
		QtOpenGL.QGLWidget.updateGL(self)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		pixeldata = glReadPixels(1,1,self.camera.width,self.camera.height,GL_RGBA,GL_UNSIGNED_BYTE)
		# Then move back
		self.setTransform(current_xform)
		self.camera.setClipNear(oldnear)
		self.camera.setClipFar(oldfar)
		self.reset_camera = True	# Reset the camera when rendering the real scene
		self.setAutoBufferSwap(True)
		self.zslicemode = False
		self.pixels = []
		self.pixels.append(1)
		self.pixels.append(1)
		self.pixels.append(self.camera.width)
		self.pixels.append(self.camera.height)
		self.pixels.append(pixeldata)
		
	def saveSnapShot(self, filename, format="tiff"):
		"""
		Save the frame buffer to an image file
		@param filename The Filename you want to save to
		@param format The image file format
		"""
		image = self.grabFrameBuffer()
		fregex = re.compile('\.\w{3,4}$')
		if re.findall(fregex, filename):
			image.save(filename, re.findall(fregex, filename)[0][1:])
		else:
			filename = "%s.%s"%(filename,format)
			image.save(filename, format)
		print("Saved %s to disk"%os.path.basename(str(filename)))
	
	def saveMovie_model(self, filename):

		tmpnames=[]
		fmt= str(filename).split(':')
		fmtstr=""
		oscillate=False
		if len(fmt)>1:
			filename=fmt[0]
			for f in fmt[1:]:
				if f=="oscillate":
					oscillate=True
				else:
					f=f.replace('=',' ')
					fmtstr+=" -{}".format(f)
		
		nimgs=len(self.getChildren())
		for mdl in range(nimgs):
		
			for i,child in enumerate(self.getChildren()):
				if i==mdl: child.setVisibleItem(True)
				else: child.setVisibleItem(False)
			
			#self.main_3d_inspector.updateTreeSelVis()
			self.updateSG()
			self.paintGL()
			image = self.grabFrameBuffer()
			tmpname = "tmp.{:03d}.png".format(mdl)
			image.save(tmpname, "png")
			tmpnames.append(tmpname)
			if oscillate:
				tmpname = "tmp.{:03d}.png".format(2*nimgs-mdl-1)
				image.save(tmpname, "png")
				tmpnames.append(tmpname)
		
		cmd="convert tmp.???.png {} {}".format(fmtstr, filename)

		ret= os.system(cmd)
		for t in tmpnames:
			os.unlink(t)
		#print "Saved %s to disk"%os.path.basename(str(filename))
		
		
	def saveMovie_rotate(self, filename):
		
		
		dt=3.
		tmpnames=[]
		fmt= str(filename).split(':')
		fmtstr=""
		oscillate=False
		rg=360
		if len(fmt)>1:
			filename=fmt[0]
			for f in fmt[1:]:
				if f=="oscillate":
					oscillate=True
				elif f.startswith("dt="):
					dt=float(f[3:])
				elif f.startswith("range="):
					rg=float(f[6:])
				else:
					f=f.replace('=',' ')
					fmtstr+=" -{}".format(f)
		
		nimgs=rg//dt
		print("will make {} images with {:.1f} angular step covering {:d} degrees".format(nimgs, dt, int(rg)))
		for i in range(nimgs):
		
			self.updateMatrices([dt,0,1,0], "rotate")
			
			self.updateSG()
			self.paintGL()
			image = self.grabFrameBuffer()
			tmpname = "tmp.{:03d}.png".format(i)
			image.save(tmpname, "png")
			tmpnames.append(tmpname)
			if oscillate:
				tmpname = "tmp.{:03d}.png".format(2*nimgs-i-1)
				image.save(tmpname, "png")
				tmpnames.append(tmpname)
		
		cmd="convert tmp.???.png {} {}".format(fmtstr, filename)

		ret= os.system(cmd)
		for t in tmpnames:
			os.unlink(t)
		#print "Saved %s to disk"%os.path.basename(str(filename))
		
	def insertNewNode(self, name, node, parentnode=None, parentidx=None):
		"""
		Insert a new node in the SG, also takes care of inspector
		if parent node is specified the node is inserted as a child of the parent
		This function should be used to add nodes to the tree b/c it determines where in the tree the node should be inserted
		@param name The node name
		@param node the node itself
		@param parentnode, the parent node if there is one
		@param the index to insert the child in. NOne, means append child to end of children 
		"""
		insertionpoint = None
		node.setLabel(name)
		if parentnode:
			if parentidx:
				parentnode.insertChild(node, parentidx)
			else:
				parentnode.addChild(node)
		else:
			insertionpointitem = self.getCurrentSelection()
			if insertionpointitem.getParent():
				insertionpointitem.getParent().insertChild(node, insertionpointitem.getParent().getChildIndex(insertionpointitem)+1)
			else:
				insertionpointitem.addChild(node)
		self.updateTree(node)
				
	def loadSession(self, filename):
		"""
		Loads a session
		@param filename Name of the session file
		"""
		rfile = open(filename, 'rb')
		try:
			tree = pickle.load(rfile)
		except:
			print("ERROR!!! Couldn't load the session file")
			rfile.close()
			return	
		rfile.close()
		
		self.makeCurrent()
		# Clear old tree
		children = tuple(self.getChildren()) # Need to create immutable so that list chaos does not ensue with pruning down the list
		for child in children:
			self.removeChild(child)
		
		self.parentnodestack = [self]
		#print tree
		self._process_session_load(tree)
		self.updateTree()
		self.updateSG()
		
	def _process_session_load(self, line):
		"""
		Helper function to load a session. 
		This is a recursive function
		"""
		if type([]) == type(line):
			if len(line) > 1:
				for cline in line:
					self._process_session_load(cline)
				self.parentnodestack.pop()
				
			else:
				# These nodes are leaves
				item3dobject = eval(line[0]["CONSTRUCTOR"])
				if line[0]["NODETYPE"] == "DataChild": # Data child need to have a parent set
					item3dobject.setParent(self.parentnodestack[-1:][0])
					item3dobject.dataChanged()
				item3dobject.setUsingDictionary(line[0])
				self.parentnodestack[-1:][0].addChild(item3dobject)

		else:
			# These nodes have children to process
			if line["CONSTRUCTOR"] != "SG": 
				item3dobject = eval(line["CONSTRUCTOR"])
				item3dobject.setUsingDictionary(line)
				self.parentnodestack[-1:][0].addChild(item3dobject)
				self.parentnodestack.append(item3dobject)
			else:
				# This is the SG, so treat special
				self.setUsingDictionary(line)
	
	def saveSession(self, filename):
		"""
		Save the SG as a session
		@param filename The Filename you want to save to
		"""
		if ".eman" not in filename: filename = "%s.%s"%(filename,"eman")
		wfile = open(filename, 'wb')
		treelist = self.getTreeAsList(self)
		pickle.dump(treelist, wfile)
		wfile.close()
		
	def getTreeAsList(self, item):
		"""
		Returns the SG as a list of lists along with all metadata. Used for saving a session
		item is the item you want to start at.
		This is a recursive function.
		@param item, the emItem that you want to transverse (grabs its subtree)
		"""
		childlist = []
		# Now do the recursion to build tree
		childlist.append(item.getItemDictionary()) # This gets a data dictionary of attributes pertaining to each item
		for ichild in item.getChildren():
			ichildlist = self.getTreeAsList(ichild)
			childlist.append(ichildlist)
		
		return childlist
	
	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMScene3D, self).getItemDictionary()
		dictionary.update({"AMBIENTLIGHT":self.firstlight.getAmbient(),"ANGULARLIGHTPOSITION":self.firstlight.getAngularPosition(),"CAPCOLOR":self.camera.getCapColor(),"LINKINGMODE":self.camera.getLinkingMode(),"CAPPINGMODE":self.camera.getCappingMode(),"CAMERACLIP":[self.camera.getClipNear(), self.camera.getClipFar()],"CAMERAPM":self.camera.getUseOrtho(),"CLEARCOLOR":self.getClearColor(),"PSEUDOFOVY":self.camera.getPseudoFovyWidth(),"FOVY":self.camera.getFovy()})
		return dictionary
	
	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMScene3D, self).setUsingDictionary(dictionary)
		# Set up the light
		self.firstlight.setAmbient(dictionary["AMBIENTLIGHT"][0],dictionary["AMBIENTLIGHT"][1],dictionary["AMBIENTLIGHT"][2],dictionary["AMBIENTLIGHT"][3])
		self.firstlight.setAngularPosition(dictionary["ANGULARLIGHTPOSITION"][0], dictionary["ANGULARLIGHTPOSITION"][1])
		# Set up the Camera
		self.camera.setClipNear(dictionary["CAMERACLIP"][0])
		self.camera.setClipFar(dictionary["CAMERACLIP"][1])
		if dictionary["CAMERAPM"]: 
			self.camera.useOrtho(self.camera.getZclip(disregardvv=True))
		else:
			objectsize = 50
			self.camera.usePerspective(objectsize, 0.25, 60.0)
		try:
			self.camera.setCapColor(*dictionary["CAPCOLOR"])
			self.camera.setCappingMode(dictionary["CAPPINGMODE"])
			self.camera.setLinkingMode(dictionary["LINKINGMODE"])
			self.camera.setFovy(dictionary["FOVY"])
			self.camera.setPseudoFovy(dictionary["PSEUDOFOVY"])
		except:
			pass
		self.cameraNeedsanUpdate()
		# Set up the utils
		self.setClearColor(dictionary["CLEARCOLOR"][0],dictionary["CLEARCOLOR"][1],dictionary["CLEARCOLOR"][2])
		
	def cameraNeedsanUpdate(self):
		"""
		Tell the SG to reset the camera
		"""
		self.reset_camera = True
		
	def setClearColor(self, r, g, b, a=0.0):
		"""
		Set the background color ambient
		"""
		self.clearcolor = [r, g, b, a]
		glClearColor(r, g, b, a)
	
	def getClearColor(self):
		"""
		Return the clear color
		"""
		return self.clearcolor
	
	def setCurrentSelection(self, item):
		"""
		This is called by the inspector or the selector to record the current selected item
		"""
		self.currentselecteditem = item
	
	def getCurrentSelection(self):
		"""
		Return the currently selected item
		"""
		return self.currentselecteditem
		
	def setAPix(self, apix):
		"""
		Set the sg angstrom per pixeldata
		"""
		self.apix = apix
		
	def getAPix(self):
		"""
		Return the sg angstrom per pixeldata. Only for ortho mode as apix doesn't have any meaning in perspective mode
		"""
		if self.camera.getUseOrtho():
			return self.apix
	
	def initialViewportDims(self, dims):
		"""
		Scale viewport to dims when the EMScene is first displayed. This only works for ortho mode, so support is needed for perspective mode
		The initial size is set when the showEvent is called, and the scaling is to the viewport width
		"""
		self.initialviewportdims = dims
		
	def updateSG(self):
		"""
		Update the SG
		"""
		self.update()
		self.updateInspector()
	
	def updateInspector(self):
		"""
		Update the inspector, if there is one
		"""
		if self.main_3d_inspector: self.main_3d_inspector.updateInspector()
		
	def updateTreeSelVis(self, selecteditem=None):
		"""
		Update the inspector Tree for vis sel attributes. Much faster than updateTree, but only checks vis sel attrs
		"""
		if self.main_3d_inspector: self.main_3d_inspector.updateTreeSelVis(selecteditem)
	
	def updateTree(self, currentnode=None):
		"""
		Update the entire inspector tree if there is one
		"""
		if self.main_3d_inspector: self.main_3d_inspector.updateTree(currentnode)
	
	def closeEvent(self, event):
		"""
		Close the main inspector
		"""
		if self.main_3d_inspector: self.main_3d_inspector.close()
		EMGLWidget.closeEvent(self,event)
	
	def showEvent(self, event):
		QtWidgets.QWidget.showEvent(self, event)
		# if desired scale viewport to a chosen scale so that initialviewportdims fills the screen
		if self.initialviewportdims != None:
			self.camera.scaleToDims(self.initialviewportdims)
			self.camera.update()
		
	# Maybe add methods to control the lights

class EMLight(object):
	def __init__(self, light):
		"""
		@type light: GL_LIGHTX, where 0 =< X <= 8
		@param light: an OpenGL light
		The light properties are set to reasonable defaults.
		"""
		self.light = light
		self.setAmbient(0.3, 0.3, 0.3, 1.0)		# Default ambient color is light grey
		self.setGlobalAmbient(0.0, 0.0, 0.0, 1.0)		# Default global ambient
		self.setDiffuse(1.0, 1.0, 1.0, 1.0)		# Default diffuse color is white
		self.setSpecualar(1.0, 1.0, 1.0, 1.0)		# Default specular color is white
		self.setAngularPosition(-30.0, 45.0)		# Default angular position
		if not glIsEnabled(GL_LIGHTING):
			glEnable(GL_LIGHTING)

	def setAmbient(self, r, g, b, a):
		"""
		@param r: the red component of the ambient light
		@param g: the green component of the ambient light
		@param b: the blue component of the ambient light
		@param a: the alpha component of the ambient light
		Set the ambient light color
		"""
		self.colorambient = [r, g, b, a]
		glLightfv(self.light, GL_AMBIENT, self.colorambient)
		

	def setGlobalAmbient(self, r, g, b, a):
		"""
		@param r: the red component of the ambient light
		@param g: the green component of the ambient light
		@param b: the blue component of the ambient light
		@param a: the alpha component of the ambient light
		Set the global ambient light color
		"""
		self.colorglobalambient = [r, g, b, a]
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, self.colorglobalambient)	# We want no global ambient light
		
	def setDiffuse(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the diffuse light color
		"""
		self.colordiffuse = [r, g, b, a]
		glLightfv(self.light, GL_DIFFUSE, self.colordiffuse)
		
	def setSpecualar(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the specular light color
		"""
		self.colorspecular = [r, g, b, a]
		glLightfv(self.light, GL_SPECULAR, self.colorspecular)

	def setPosition(self, x, y, z, w):
		"""
		@param x: The x component of the light position
		@param y: The y component of the light position
		@param z: The z component of the light position
		@param w: The w component of the light position
		Set the light position, in homogeneous coords
		"""
		self.position = [x, y, z, w]
		glLightfv(self.light, GL_POSITION, self.position)
		
	def setAngularPosition(self, theta, phi):
		"""
		@param theta: The theta component of the light position in spherical coords
		@param phi: The theta component of the light position in spherical coords
		Set the light position in spherical coords. This is only used by the lightwidget in the inspector
		"""
		z = math.sin(math.radians(theta + 90))*math.cos(math.radians(phi))
		y = math.sin(math.radians(theta + 90))*math.sin(math.radians(phi))
		x = math.cos(math.radians(theta + 90))
		self.setPosition(x, y, z, 0.0)
		self.angularposition = [theta, phi]
		
	def getAngularPosition(self):
		"""
		Return the light position as spherical coords
		"""
		return self.angularposition
		
	def getAmbient(self):
		"""
		Return the ambient lighting
		"""
		return self.colorambient
	
	def getGlobalAmbient(self):
		"""
		Return the global ambient color
		"""
		return self.colorglobalambient
		
	def getDiffuse(self):
		"""
		Return the diffuse lighting
		"""
		return self.colordiffuse
		
	def getSpecular(self):
		"""
		Return the specular lighting
		"""
		returnself.colorspecular
	
	def getPosition(self):
		"""
		Return the light position
		"""
		return self.position
		
	def enableLighting(self):
		"""
		Enables this light
		"""
		if not glIsEnabled(self.light):
			glEnable(self.light)

	def disableLighting(self):
		"""
		Disables this light
		"""
		if glIsEnabled(self.light):
			glDisable(self.light)
			
	def updateGLlighting(self):
		"""
		With open GL, Stupid is a stupid does :)
		Sometimes we need to re-issue glLighting commands to get GLto execute them, I don't know why this happens......
		"""
		self.setAmbient(self.colorambient[0], self.colorambient[1], self.colorambient[2], self.colorambient[3])
		self.setGlobalAmbient(self.colorglobalambient[0], self.colorglobalambient[1], self.colorglobalambient[2], self.colorglobalambient[3])
		
class EMCamera(object):
	"""Implementation of the camera"""
	def __init__(self, near, far, usingortho=True, fovy=60.0, boundingbox=50.0, screenfraction=0.5):
		"""
		@param fovy: The field of view angle
		@param near: The volume view near position
		@param far: The volume view far position
		@param usingortho: Use orthographic projection
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		"""
		self.pseudofovy = 0
		self.far = far
		self.near = near
		self.fovy = fovy
		self.maxviewport = glGetIntegerv(GL_MAX_VIEWPORT_DIMS)
		self.setCappingMode(False)
		self.setCapColor(*(get_default_gl_colors()["bluewhite"]['ambient']))
		self.setLinkingMode(False)
		zclip = old_div((self.near-self.far),2.0)	# Puts things in the center of the viewing volume
		if usingortho:
			self.useOrtho(zclip)
		else:
			self.usePerspective(boundingbox, screenfraction, fovy)

	def update(self, width=None, height=None):
		"""
		@param width: The width of the window in pixels
		@param height: The height of the window in pixels
		updates the camera and viewport after window resize
		"""
		if width: self.width = width
		if height: self.height = height
		self.aspectratio = old_div(float(self.height),float(self.width))
		if self.usingortho:
			# this deals with maxviewport better than the previous solution
			# the whole pseudofov concept is really messed up. Needs a complete rethink
			if self.width+2*self.pseudofovy>=self.maxviewport[0] :
				self.pseudofovy=(self.maxviewport[0]-self.width)//2-1
			self.setViewPort(-self.pseudofovy, int(-self.pseudofovy*self.aspectratio), int(self.width+2*self.getPseudoFovyWidth()), int(self.height+2*self.getPseudoFovyHeight()))
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.setOrthoProjectionMatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			self.setCameraPosition()
		else:
			# This may need some work to get it to behave
			self.setViewPort(0,0,self.width,self.height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.setPerspectiveProjectionMatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			self.setCameraPosition()
		
		
	def scaleToDims(self, dims):
		"""
		Scale in, ortho mode to the given dimensions. Scale is to width
		"""
		self.setPseudoFovy(old_div((self.width*self.width),(2*dims)) - (old_div(self.width,2)))
		
	def setViewPort(self, x, y, vpwidth, vpheight):
		"""Set the viewport subject to openGL constraints """
		if (vpwidth < self.maxviewport[0] and vpheight < self.maxviewport[1]):
			glViewport(x, y, vpwidth, vpheight)
	#		print(x,y,vpwidth, vpheight)
		
	def setCameraPosition(self, sfactor=1):
		"""
		Set the default camera position
		"""
		glTranslate(0,0,sfactor*self.getZclip())
		
	def setProjectionMatrix(self):
		"""
		Set the projection matrix
		"""
		if self.usingortho:
			self.setOrthoProjectionMatrix()
		else:
			self.setPerspectiveProjectionMatrix()
			
	def setOrthoProjectionMatrix(self):
		"""
		Set the orthographic projection matrix. Volume view origin (0,0) is center of screen
		"""
		glOrtho(-self.width//2, self.width//2, -self.height//2, self.height//2, self.near, self.far)
		
	def setPerspectiveProjectionMatrix(self):
		"""
		Set the perspective projection matrix. Volume view origin (0,0) is center of screen
		"""
		GLU.gluPerspective(self.fovy, (float(self.width)/float(self.height)), self.near, self.far)
			
	def usePerspective(self, boundingbox, screenfraction, fovy=60.0):
		""" 
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		Changes projection matrix to perspective
		"""
		self.fovy = fovy
		self.perspective_z = old_div(-(old_div(boundingbox*2,screenfraction)),(2*math.tan(math.radians(old_div(self.fovy,2)))))  + boundingbox
		self.usingortho = False
		
	def useOrtho(self, zclip):
		"""
		Changes projection matrix to orthographic
		"""
		self.usingortho = True
		self.zclip = zclip
		
	def getUseOrtho(self):
		"""
		Returns the projection state
		"""
		return self.usingortho
		
	def setClipFar(self, far):
		"""
		Set the far aspect of the viewing volume
		"""
		self.far = far

	def setClipNear(self, near):
		"""
		Set the near aspect of the viewing volume
		"""
		self.near = near

	def getClipNear(self):
		"""
		Get the near clipping plane
		"""
		return self.near
		
	def getClipFar(self):
		"""
		Get the far clipping plane
		"""
		return self.far
		
	def setFovy(self, fovy):
		"""
		Set the field of view angle aspect of the viewing volume
		"""
		self.fovy = fovy
	
	def getFovy(self):
		""" 
		Return FOVY
		"""
		return self.fovy
		
	def setPseudoFovy(self, pseudofovy):
		"""
		Set PseudoFovy, a sort of fovy for orthographic projections, do bounds checking
		"""
		#if ((self.width+2*pseudofovy) > 0 and (self.height+2*pseudofovy) > 0):
		if ((int(self.width+pseudofovy) < self.maxviewport[0] and int(self.height+pseudofovy*self.aspectratio) < self.maxviewport[1]) or pseudofovy < self.pseudofovy): 
			self.pseudofovy = pseudofovy
		else:
			# Set to max zoom
			self.pseudofovy = (self.maxviewport[0] - self.width)
	
		#print(pseudofovy,self.pseudofovy,self.maxviewport)

	def getPseudoFovyWidth(self):
		"""
		Return PseudoFovy, a sort of fovy for orthographic projections
		"""
		return self.pseudofovy
	
	def getPseudoFovyHeight(self):
		"""
		Return PseudoFovy, a sort of fovy for orthographic projections
		"""
		return self.pseudofovy*self.aspectratio
		
	def getViewPortWidthScaling(self):
		"""
		Return the scaling necessary to insert move from the viewport to the actual viewport. Need this b/c of the crazy scaling scheme
		"""
		return old_div(float(self.getWidth()),float(self.getWidth() + 2*self.getPseudoFovyWidth()))
		
	def getViewPortHeightScaling(self):
		"""
		Return the scaling necessary to insert move from the viewport to the actual viewport. Need this b/c of the crazy scaling scheme
		"""
		return old_div(float(self.getHeight()),float(self.getHeight() + 2*self.getPseudoFovyHeight()))
		
	def getHeight(self):
		""" Get the viewport height """
		return self.height
	
	def getWidth(self):
		""" Get the viewport width """
		return self.width
		
	def getZclip(self, disregardvv=False):
		""" Get the zclip """
		if self.usingortho or disregardvv:
			return self.zclip
		else:
			return self.perspective_z
			
	def setZclip(self, clip):
		""" Set the Z clip """
		if self.usingortho:
			self.zclip = clip
		else:
			self.perspective_z = clip
			
	def getCappingMode(self):
		""" Return the capping mode """
		return self.cappingmode
		
	def setCappingMode(self, mode):
		self.cappingmode = mode
		
	def setCapColor(self, r, g, b, a=1.0):
		""" Set the cap color """
		self.capcolor = [r, g, b, a]
		
	def getCapColor(self):
		""" Return the cap color """
		return self.capcolor
		
	def setLinkingMode(self, mode):
		""" Set the clipping plane linking mode"""
		self.linkingmode = mode
		
	def getLinkingMode(self):
		""" Return the linking mode """
		return self.linkingmode
		
###################################### Inspector Code #########################################################################################

class EMInspector3D(QtWidgets.QWidget):
	def __init__(self, scenegraph):
		"""
		The inspector for the 3D widget. The inspector is a strict observer of the SceneGraph, and is updated by calling update inspector
		"""
		QtWidgets.QWidget.__init__(self)
		self.scenegraph = weakref.ref(scenegraph)
		self.mintreewidth = 250		# minimum width of the tree
		self.mincontrolwidth = 0
		
		vbox = QtWidgets.QVBoxLayout(self)
		self.inspectortab = QtWidgets.QTabWidget()
		self.inspectortab.addTab(self.getTreeWidget(), "Tree View")
		self.inspectortab.addTab(self.getLightsWidget(), "Lights")
		self.inspectortab.addTab(self.getCameraWidget(), "Camera")
		self.inspectortab.addTab(self.getUtilsWidget(), "Utils")
		toolframe = QtWidgets.QFrame()
		toolframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
		toolframe.setLayout(self._get_toolbox_layout())
		vbox.addWidget(self.inspectortab)
		vbox.addWidget(toolframe)
		
		self.inspectortab.currentChanged[int].connect(self._on_load_camera)
		
		self.setLayout(vbox)
		self.updateGeometry()
	
	def closeEvent(self, event):
		""" Clean up and close """
		self.scenegraph().main_3d_inspector = None
		# There is a BUG in QStackedWidget @^#^&#, so it thinks that widgets have been deleted when they haven't!!! (It thinks that when you delete the stacked widget all widgets in the stack have been removed when in fact that is not always the case)
		for node in self.scenegraph().getAllNodes():
			node.item_inspector = None
			
	def getTreeWidget(self):
		"""
		This returns the treeview-control panel widget
		"""
		widget = QtWidgets.QWidget()
		hbox = QtWidgets.QHBoxLayout(widget)
		treeframe = QtWidgets.QFrame()
		treeframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
		treeframe.setLayout(self._get_tree_layout(widget))
		treeframe.setMinimumWidth(self.mintreewidth)
		hbox.addWidget(treeframe)
		self.stacked_widget = QtWidgets.QStackedWidget()
		self.stacked_widget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		widget.setLayout(hbox)
		
		return widget
		
	def _get_tree_layout(self, parent):
		"""
		Returns the tree layout
		"""
		tvbox = QtWidgets.QVBoxLayout()
		self.tree_widget = EMQTreeWidget(parent)
		self.tree_widget.setHeaderLabel("Choose a item")
		tvbox.addWidget(self.tree_widget)
		self.tree_node_button_add = QtWidgets.QPushButton("Add Object")
		self.tree_node_button_remove = QtWidgets.QPushButton("Remove Object")
		self.tree_node_slider = ValSlider(label="Seq:")
		self.tree_node_slider.setIntonly(True)
		self.tree_node_slider.setRange(0,1)
		self.tree_node_slider.setValue(0)
		tvbox.addWidget(self.tree_node_button_add)
		tvbox.addWidget(self.tree_node_button_remove)
		tvbox.addWidget(self.tree_node_slider)
		
		self.tree_widget.itemClicked[QtWidgets.QTreeWidgetItem, int].connect(self._tree_widget_click)
		self.tree_widget.visibleItem[QtWidgets.QTreeWidgetItem].connect(self._tree_widget_visible)
		self.tree_widget.editItem[QtWidgets.QTreeWidgetItem].connect(self._tree_widget_edit)
		self.tree_node_button_remove.clicked.connect(self._tree_widget_remove)
		self.tree_node_button_add.clicked.connect(self._on_add_button)
		self.tree_node_slider.valueChanged.connect(self._slider_change)
		
		return tvbox
	
	def _recursiveupdatetreeselvis(self, item):
		item.setSelectionStateBox()
		item.getVisibleState()
		for childidx in range(item.childCount()):
			self._recursiveupdatetreeselvis(item.child(childidx))
			
	def updateTreeSelVis(self, selecteditem=None):
		"""
		Update the selection and visibility states. Makes the Sel Vis states an observer of the SG
		"""
		# Update the tree
		self._recursiveupdatetreeselvis(self.tree_widget.topLevelItem(0))
		# Set the desired item if desired
		if selecteditem:
			try:
				self.stacked_widget.setCurrentWidget(selecteditem.getItemInspector())
				self.tree_widget.setCurrentItem(selecteditem.EMQTreeWidgetItem)
				if selecteditem: self.scenegraph().setCurrentSelection(selecteditem)
			except:
				pass
			# Unsure unique selection
			self.ensureUniqueTreeLevelSelection(selecteditem)
			
	def ensureUniqueTreeLevelSelection(self, item):
		"""
		Make sure that we don't select both an ancestor and child at the same time
		"""
		for ancestor in item.getSelectedAncestorNodes():
			if ancestor.EMQTreeWidgetItem:			# Not al ancestors are listed on the inspector tree (such as a data node)
				ancestor.EMQTreeWidgetItem.setSelectionState(False)
		for child in item.getAllSelectedNodes()[1:]: 	# Lop the node itself off
			child.EMQTreeWidgetItem.setSelectionState(False)
	
	def _slider_change(self):
		mdl=int(self.tree_node_slider.getValue())
		for i,child in enumerate(self.scenegraph().getChildren()):
			if i==mdl: child.setVisibleItem(True)
			else: child.setVisibleItem(False)
		
		self.updateTreeSelVis()
		self.scenegraph().updateSG()

	def _recursiveAdd(self, parentitem, parentnode,depth=0):
		"""
		Helper function to load the SG
		"""
		for child in parentnode.getChildren():
			if not child.getLabel(): child.setLabel(child.name)
			addeditem = self.addTreeNode(child.getLabel(), child, parentitem)
			self._recursiveAdd(addeditem, child,depth+1)
		# Expand the data items
		if parentitem.childCount() > 0: parentitem.setExpanded(True)
		self.tree_node_slider.setRange(0,len(parentnode.getChildren())-1)
		
	def loadSG(self):
		"""
		Load the SG
		"""
		rootitem = self.addTreeNode("All Objects", self.scenegraph())
		self._recursiveAdd(rootitem, self.scenegraph())
		
	def addTreeNode(self, name, item3d, parentitem=None, insertionindex=-1):
		"""
		Add a node (item3d) to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the node, so it can talk to the node.
		You can think of this as a three way conversation (the alternative it to use a mediator, but that is not worth it w/ only three players)
		"""
		tree_item = EMQTreeWidgetItem(name, item3d, parentitem)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI
		item3d.setEMQTreeWidgetItem(tree_item)				# Reference to the EMQTreeWidgetItem
		item_inspector = item3d.getItemInspector()				# Get the node GUI controls 
		#return tree_item
		item_inspector.setInspector(self)					# Associate the item GUI with the inspector
		self.stacked_widget.addWidget(item_inspector)			# Add a widget to the stack
		item3d.setLabel(name)						# Set the label
		# Set icon status
		tree_item.setSelectionStateBox()
		# Set parent if one exists	
		if not parentitem:
			self.tree_widget.insertTopLevelItem(0, tree_item)
		else:
			if insertionindex >= 0:
				parentitem.insertChild(insertionindex, tree_item)
			else:
				parentitem.addChild(tree_item)
		return tree_item
	
	def removeTreeNode(self, parentitem, childindex):
		# I am using the parent item rather than the item itself b/c the stupid widget has no , remove self function...
		# Remove both the QTreeWidgetItem and the widget from the WidgetStack, otherwise we'll get memory leaks 
		if parentitem.child(childindex).item3d():
			self.stacked_widget.removeWidget(parentitem.child(childindex).item3d().getItemInspector())
		parentitem.takeChild(childindex)
	
	def clearTree(self):
		"""
		Clear the entire tree
		"""
		if self.tree_widget.topLevelItem(0):
			self.tree_widget.topLevelItem(0).removeAllChildren(self)
			self.tree_widget.takeTopLevelItem(0)
		
	def _tree_widget_click(self, item, col, quiet=False):
		"""
		When a user clicks on the selection tree check box
		"""
		self.stacked_widget.setCurrentWidget(item.item3d().getItemInspector())
		item.setSelectionState(item.checkState(0))
		# This code is to prevent both descendants and children from being selected....
		if item.checkState(0) == QtCore.Qt.Checked: self.ensureUniqueTreeLevelSelection(item.item3d())
		if not item.item3d().isSelectedItem(): item.item3d().getItemInspector().updateItemControls() # This is too update a widget, translation and rotation may change in parent nodes change
		self.scenegraph().setCurrentSelection(item.item3d())
		if not quiet: self.updateSceneGraph()
		
	def _tree_widget_visible(self, item):
		"""
		When a user clicks on the visible icon
		"""
		item.toggleVisibleState()
		self.updateSceneGraph()
	
	def _tree_widget_edit(self):
		"""
		When a use middle clicks
		"""
		nodedialog = NodeEditDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
	
	def _on_add_button(self):
		nodedialog =  NodeDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
		
	def _tree_widget_remove(self):
		"""
		When a use wants to remove a node_name
		"""
		item = self.tree_widget.currentItem()
		if item.parent:
			self.removeTreeNode(item.parent(), item.parent().indexOfChild(item)) 
			item.parent().item3d().removeChild(item.item3d())
			# In case we delete the currently selected item, we want6 to move selected item to last selection
			if self.scenegraph().getCurrentSelection() == item.item3d():
				self.scenegraph().setCurrentSelection(self.tree_widget.currentItem().item3d())
			self.updateSceneGraph()
		else:
			print("Error cannot remove root node!!")
			
	def _get_toolbox_layout(self):
		tvbox = QtWidgets.QHBoxLayout()
		font = QtGui.QFont()
		font.setBold(True)
		toollabel = QtWidgets.QLabel("Tools")
		toollabel.setFont(font)
		self.rotatetool = EMANToolButton()
		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
		self.rotatetool.setToolTip("Rotate X/Y\nMouse: Right 'n' drag\nHot Key: R")
		self.translatetool =EMANToolButton()
		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
		self.translatetool.setToolTip("Translate X/Y\nMouse: Left 'n' drag\nHot Key: T")
		self.ztranslate = EMANToolButton()
		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
		self.ztranslate.setToolTip("Translate Z\nHot Key: Z")
		self.scaletool = EMANToolButton()
		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
		self.scaletool.setToolTip("Scale\nHot Key: S")
		self.rulertool = EMANToolButton()
		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
		self.rulertool.setToolTip("Ruler\nHot Key: U\nDoes NOT account for scaling")
		self.selectiontool = EMANToolButton()
		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
		self.selectiontool.setToolTip("Select objects\nMouse: Left 'n' drag\nMultiple = + Shift\nHot Key: Esc")
		self.multiselectiontool = EMANToolButton()
		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
		self.multiselectiontool.setToolTip("Select multiple objects\nMouse: Left 'n' drag\nHot Key: M")
		self.linetool = EMANToolButton()
		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Insert Line\nHot Key: L")
		self.cubetool = EMANToolButton()
		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
		self.cubetool.setToolTip("Insert Cube\nHot Key: C")
		self.spheretool = EMANToolButton()
		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
		self.spheretool.setToolTip("Insert Sphere\nHot Key: P")
		self.cylindertool = EMANToolButton()
		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
		self.cylindertool.setToolTip("Insert Cylinder\nHot Key: Y")
		self.conetool = EMANToolButton()
		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
		self.conetool.setToolTip("Insert Cone\nHot Key: O")
		self.texttool = EMANToolButton()
		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
		self.texttool.setToolTip("Insert Text\nHot Key: X")
		self.datatool = EMANToolButton()
		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
		self.datatool.setToolTip("Insert Data\nHot Key: E")
		self.apptool = EMANToolButton()
		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
		self.apptool.setToolTip("Application dependent\nHot Key: A")
		
		tvbox.addWidget(toollabel)
		tvbox.addWidget(self.selectiontool)
		tvbox.addWidget(self.multiselectiontool)
		tvbox.addWidget(self.translatetool)
		tvbox.addWidget(self.ztranslate)
		tvbox.addWidget(self.rotatetool)
		tvbox.addWidget(self.scaletool)
		tvbox.addWidget(self.rulertool)
		tvbox.addWidget(self.linetool)
		tvbox.addWidget(self.cubetool)
		tvbox.addWidget(self.spheretool)
		tvbox.addWidget(self.cylindertool)
		tvbox.addWidget(self.conetool)
		tvbox.addWidget(self.texttool)
		tvbox.addWidget(self.datatool)
		tvbox.addWidget(self.apptool)
		tvbox.setAlignment(QtCore.Qt.AlignLeft)
		
		self.rotatetool.clicked[int].connect(self._rotatetool_clicked)
		self.translatetool.clicked[int].connect(self._transtool_clicked)
		self.ztranslate.clicked[int].connect(self._ztranstool_clicked)
		self.scaletool.clicked[int].connect(self._scaletool_clicked)
		self.rulertool.clicked[int].connect(self._rulertool_clicked)
		self.selectiontool.clicked[int].connect(self._seltool_clicked)
		self.multiselectiontool.clicked[int].connect(self._multiseltool_clicked)
		self.linetool.clicked[int].connect(self._linetool_clicked)
		self.cubetool.clicked[int].connect(self._cubetool_clicked)
		self.spheretool.clicked[int].connect(self._spheretool_clicked)
		self.cylindertool.clicked[int].connect(self._cylindertool_clicked)
		self.conetool.clicked[int].connect(self._conetool_clicked)
		self.texttool.clicked[int].connect(self._texttool_clicked)
		self.datatool.clicked[int].connect(self._datatool_clicked)
		self.apptool.clicked[int].connect(self._apptool_clicked)
			
		return tvbox
	
	def _rotatetool_clicked(self, state):
		self.scenegraph().setMouseMode("rotate")
		
	def _transtool_clicked(self, state):
		self.scenegraph().setMouseMode("xytranslate")
		
	def _ztranstool_clicked(self, state):
		self.scenegraph().setMouseMode("ztranslate")
		
	def _scaletool_clicked(self, state):
		self.scenegraph().setMouseMode("scale")
	
	def _rulertool_clicked(self, state):
		self.scenegraph().setMouseMode("ruler")
		
	def _seltool_clicked(self, state):
		self.scenegraph().setMouseMode("selection")
		
	def _multiseltool_clicked(self, state):
		self.scenegraph().setMouseMode("multiselection")
	
	def _linetool_clicked(self, state):
		self.scenegraph().setMouseMode("line")
		
	def _cubetool_clicked(self, state):
		self.scenegraph().setMouseMode("cube")
		
	def _spheretool_clicked(self, state):
		self.scenegraph().setMouseMode("sphere")
		
	def _cylindertool_clicked(self, state):
		self.scenegraph().setMouseMode("cylinder")
	
	def _conetool_clicked(self, state):
		self.scenegraph().setMouseMode("cone")
		
	def _texttool_clicked(self, state):
		self.scenegraph().setMouseMode("text")
		
	def _datatool_clicked(self, state):
		self.scenegraph().setMouseMode("data")
				
	def _apptool_clicked(self, state):
		self.scenegraph().setMouseMode("app")
		
	def getLightsWidget(self):
		"""
		Returns the lights control widget
		"""
		self.lighttab_open = False
		lwidget = QtWidgets.QWidget()
		lvbox = QtWidgets.QVBoxLayout()
		lightslabel = QtWidgets.QLabel("Lights", lwidget)
		lightslabel.setAlignment(QtCore.Qt.AlignCenter)
		lightslabel.setMaximumHeight(30.0)
		font = QtGui.QFont()
		font.setBold(True)
		lightslabel.setFont(font)
		lvbox.addWidget(lightslabel)
		self.lightwidget = EMLightControls(GL_LIGHT1)
		positionlabel = QtWidgets.QLabel("Position", lwidget)
		positionlabel.setMaximumHeight(20.0)
		positionlabel.setAlignment(QtCore.Qt.AlignCenter)
		valslidersplitter = QtWidgets.QFrame()
		valslidersplitter.setFrameShape(QtWidgets.QFrame.StyledPanel)
		valslidersplitter.setMaximumHeight(80)
		valvbox = QtWidgets.QVBoxLayout()
		self.hvalslider = ValSlider(lwidget,(0.0,360.0),"Horizontal")
		self.vvalslider = ValSlider(lwidget,(0.0,360.0),"Vertical")
		valvbox.addWidget(self.hvalslider)
		valvbox.addWidget(self.vvalslider)
		valslidersplitter.setLayout(valvbox)
		lvbox.addWidget(self.lightwidget)
		lvbox.addWidget(positionlabel)
		lvbox.addWidget(valslidersplitter)
		self.ambientlighting = ValSlider(lwidget,(0.0,1.0),"Ambient Lighting")
		self.ambientlighting.setMaximumHeight(30)
		lvbox.addWidget(self.ambientlighting)
		lwidget.setLayout(lvbox)
		
		self.lightwidget.lightPositionMoved.connect(self._light_position_moved)
		self.hvalslider.valueChanged.connect(self._on_light_slider)
		self.vvalslider.valueChanged.connect(self._on_light_slider)
		self.ambientlighting.valueChanged.connect(self._on_light_ambient)

		return lwidget
	
	def _light_position_moved(self, position): 
		self.scenegraph().firstlight.setAngularPosition(position[0], position[1])
		self.scenegraph().updateSG()
	
	def _on_light_slider(self, value):
		self.lightwidget.setAngularPosition(self.hvalslider.getValue(), self.vvalslider.getValue())
		position = self.lightwidget.getPosition()
		self.scenegraph().firstlight.setPosition(position[0], position[1], position[2], position[3])
		self.scenegraph().update()
	
	def _on_light_ambient(self):
		ambient = self.ambientlighting.getValue()
		self.scenegraph().firstlight.setAmbient(ambient, ambient, ambient, 1.0)
		self.scenegraph().update()
	
	def getCameraWidget(self):
		"""
		Returns the camera control widget
		"""
		self.cameratab_open = False
		cwidget = QtWidgets.QWidget()
		grid = QtWidgets.QGridLayout()
		grid.setContentsMargins(11,5,11,5)
		self.camerawidget = CameraControls(scenegraph=self.scenegraph())
		grid.addWidget(self.camerawidget, 0, 0, 1, 3)
		nlabel = QtWidgets.QLabel("Near clipping plane", cwidget)
		nlabel.setMaximumHeight(30.0)
		nlabel.setAlignment(QtCore.Qt.AlignCenter)
		self.near = EMSpinWidget(self.scenegraph().camera.getClipNear(), 1.0)
		self.near.setToolTip("In the Window above:\nClick 'n' drag, near the near clipping, to move near clipping plane")
		self.near.setMaximumHeight(40.0)
		grid.addWidget(nlabel, 1, 0)
		grid.addWidget(self.near ,1, 1)
		flabel = QtWidgets.QLabel("Far clipping plane", cwidget)
		flabel.setMaximumHeight(30.0)
		flabel.setAlignment(QtCore.Qt.AlignCenter)
		self.far = EMSpinWidget(self.scenegraph().camera.getClipFar(), 1.0)
		self.far.setToolTip("In the Window above:\nClick 'n' drag, near the far clipping, to move far clipping plane")
		self.far.setMaximumHeight(40.0)
		grid.addWidget(flabel, 2, 0)
		grid.addWidget(self.far, 2, 1)
		
		# The frame for linking capping, etc
		cframe = QtWidgets.QFrame()
		cframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
		cframegrid = QtWidgets.QGridLayout()
		self.linkcb = QtWidgets.QCheckBox("Link clipping planes")
		self.linkcb.setMinimumHeight(40.0)
		cframegrid.addWidget(self.linkcb, 0, 0, 1, 2)
		self.capcb = QtWidgets.QCheckBox("Cap clipping planes")
		self.capcb.setMinimumHeight(40.0)
		cframegrid.addWidget(self.capcb, 1, 0, 1, 2)
		capcolorlabel = QtWidgets.QLabel("Cap Color")
		cframegrid.addWidget(capcolorlabel, 2, 0, 1, 1)
		self.cappingcolor = EMQTColorWidget(parent=self)
		cframegrid.addWidget(self.cappingcolor, 2, 1, 1, 1)
		cframegrid.setAlignment(QtCore.Qt.AlignCenter)
		cframe.setLayout(cframegrid)
		grid.addWidget(cframe, 1, 2, 3, 1)
		
		# The frame for project matrices
		frame = QtWidgets.QFrame()
		frame.setMaximumHeight(40.0)
		frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
		hbox = QtWidgets.QHBoxLayout()
		vvlabel = QtWidgets.QLabel("Viewing Volume")
		self.orthoradio = QtWidgets.QRadioButton("Orthographic")
		self.perspectiveradio = QtWidgets.QRadioButton("Perspective")
		hbox.addWidget(vvlabel)
		hbox.addWidget(self.orthoradio)
		hbox.addWidget(self.perspectiveradio)
		frame.setLayout(hbox)
		grid.addWidget(frame, 3, 0, 1, 2)
		cwidget.setLayout(grid)

		self.near.valueChanged[int].connect(self._on_near)
		self.far.valueChanged[int].connect(self._on_far)
		self.camerawidget.nearMoved[float].connect(self._on_near_move)
		self.camerawidget.farMoved[float].connect(self._on_far_move)
		self.orthoradio.clicked.connect(self._on_radio_click)
		self.perspectiveradio.clicked.connect(self._on_radio_click)
		self.capcb.clicked.connect(self._on_capping)
		self.linkcb.clicked.connect(self._on_linking)
		self.cappingcolor.newcolor[QtGui.QColor].connect(self._on_cap_color)
		
		return cwidget
		
	def _on_capping(self):
		self.scenegraph().camera.setCappingMode(self.capcb.isChecked())
		self.scenegraph().updateSG()
	
	def _on_linking(self):
		self.scenegraph().camera.setLinkingMode(self.linkcb.isChecked())
		
	def _on_cap_color(self, color):
		rgb = color.getRgb()
		self.scenegraph().makeCurrent()
		self.scenegraph().camera.setCapColor(old_div(float(rgb[0]),255.0), old_div(float(rgb[1]),255.0), old_div(float(rgb[2]),255.0))
		if self.scenegraph().camera.getCappingMode(): self.updateSceneGraph()
		
	def _on_near(self, value, link=True):
		if not self.scenegraph().camera.usingortho and value <= 0:
			return
		if self.scenegraph().camera.getClipFar() > value:
			if self.linkcb.isChecked() and link: self._on_far(self.scenegraph().camera.getClipFar() + (value - self.scenegraph().camera.getClipNear()), link=False)
			self.scenegraph().camera.setClipNear(value)
			self.scenegraph().updateSG()
	
	def _on_far(self, value, link=True):
		if value > self.scenegraph().camera.getClipNear():
			if self.linkcb.isChecked() and link: self._on_near(self.scenegraph().camera.getClipNear() + (value - self.scenegraph().camera.getClipFar()), link=False)
			self.scenegraph().camera.setClipFar(value)
			self.scenegraph().updateSG()
		
	def _on_near_move(self, movement, link=True):
		value = self.scenegraph().camera.getClipNear() + movement
		if not self.scenegraph().camera.usingortho and value <= 0:
			return
		if self.scenegraph().camera.getClipFar() > value:
			self.scenegraph().camera.setClipNear(value)
			if self.linkcb.isChecked() and link: self._on_far_move(movement, link=False)
			self.scenegraph().updateSG()
		
	def _on_far_move(self, movement, link=True):
		value = self.scenegraph().camera.getClipFar() + movement
		if value > self.scenegraph().camera.getClipNear():
			self.scenegraph().camera.setClipFar(value)
			if self.linkcb.isChecked() and link: self._on_near_move(movement, link=False)
			self.scenegraph().updateSG()
		
	def _on_load_camera(self, idx):
		"""
		Load the SG Z slice and only update widgets that are open
		"""
		if idx == 2: # '2' is the camera tab
			self.scenegraph().setZslice()
			self.cameratab_open = True
			self.updateInspector()
			return
		if idx == 1: # '1' is the  light tab
			self.lighttab_open = True
			self.updateInspector()
			return
			
		self.cameratab_open = False
		self.lighttab_open = False
			
	def _get_vv_state(self):
		"""
		Get the viewing volume of the camera
		"""
		if self.scenegraph().camera.usingortho: 
			self.orthoradio.setChecked(True)
		else:
			self.perspectiveradio.setChecked(True)
	
	def _on_radio_click(self):
		"""
		set the viewing volume. objectsize is the size of the bounding box of the largest displayed object
		"""
		objectsize = 50
		if self.orthoradio.isChecked(): self.scenegraph().camera.useOrtho(self.scenegraph().camera.getZclip(disregardvv=True))
		if self.perspectiveradio.isChecked(): self.scenegraph().camera.usePerspective(objectsize, 0.25, 60.0)
		self.scenegraph().updateSG()
		
	def getUtilsWidget(self):
		"""
		Return the utilities widget
		"""
		uwidget = QtWidgets.QWidget()
		uvbox = QtWidgets.QVBoxLayout()
		font = QtGui.QFont()
		font.setBold(True)
		# Controls frame
		frame = QtWidgets.QFrame()
		frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
		gridbox = QtWidgets.QGridLayout()
		backgroundcolor_label = QtWidgets.QLabel("Background Color", frame)
		backgroundcolor_label.setFont(font)
		self.backgroundcolor = EMQTColorWidget(parent=frame)
		self.hideselectionbutton = QtWidgets.QCheckBox("Hide Display Selections")
		self.hideselectionbutton.setMinimumHeight(100)
		self.hideselectionbutton.setFont(font)
		gridbox.addWidget(backgroundcolor_label, 0, 0)
		gridbox.addWidget(self.backgroundcolor, 0, 1)
		gridbox.addWidget(self.hideselectionbutton, 1, 0, 1, 2)
		gridbox.setAlignment(QtCore.Qt.AlignCenter)
		gridbox.setSpacing(10)
		frame.setLayout(gridbox)
		# Buttons frame
		uvbox.addWidget(frame)
		self.opensession_button = QtWidgets.QPushButton("Open Session")
		self.savesession_button = QtWidgets.QPushButton("Save Session")
		self.savebutton = QtWidgets.QPushButton("Save Image Snapshot")
		self.moviebutton0 = QtWidgets.QPushButton("Save GIF Movie (rotate)")
		self.moviebutton1 = QtWidgets.QPushButton("Save GIF Movie (model)")
		uvbox.addWidget(self.opensession_button)
		uvbox.addWidget(self.savesession_button)
		uvbox.addWidget(self.savebutton)
		uvbox.addWidget(self.moviebutton0)
		uvbox.addWidget(self.moviebutton1)
		uwidget.setLayout(uvbox)
		
		self.backgroundcolor.newcolor[QtGui.QColor].connect(self._on_bg_color)
		self.hideselectionbutton.clicked.connect(self._on_hide)
		self.savebutton.clicked.connect(self._on_save)
		self.moviebutton0.clicked.connect(self._on_save_movie_rotate)
		self.moviebutton1.clicked.connect(self._on_save_movie_model)
		self.savesession_button.clicked.connect(self._on_save_session)
		self.opensession_button.clicked.connect(self._on_open_session)
		
		return uwidget
	
	def _on_hide(self):
		"""
		Hide display selections
		"""
		for node in self.scenegraph().getAllNodes():
				node.setHiddenSelected(self.hideselectionbutton.isChecked())
		self.updateSceneGraph()
		
	def _on_open_session(self):
		"""
		Open a session
		"""
		# Open the file
		filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Session', os.getcwd(), "*.eman")[0]
		if filename:
			self.scenegraph().loadSession(filename)
		
	def _on_save_session(self):
		"""
		Return a list of all the child items (actually a tree of sorts)
		"""
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Session', os.getcwd(), "*.eman")[0]
		if filename: # if we cancel
			self.scenegraph().saveSession(filename)

	def _on_save(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.tiff *.jpeg *.png)")[0]
		if filename: # if we cancel
			self.scenegraph().saveSnapShot(filename)
	
	def _on_save_movie_rotate(self):
		"""
		Save a movie of the scene
		"""
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.gif)")[0]
		if filename: # if we cancel
			self.scenegraph().saveMovie_rotate(filename)
			
	def _on_save_movie_model(self):
		"""
		Save a movie of the scene
		"""
		filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.gif)")[0]
		if filename: # if we cancel
			self.scenegraph().saveMovie_model(filename)
	
	def _on_bg_color(self, color):
		rgb = color.getRgb()
		self.scenegraph().makeCurrent()
		self.scenegraph().setClearColor(old_div(float(rgb[0]),255.0), old_div(float(rgb[1]),255.0), old_div(float(rgb[2]),255.0))
		self.updateSceneGraph()
		
	def updateInspector(self):
		"""
		Update Inspector,is called whenever the scene changes
		"""
		#tool buttons
		if self.scenegraph().getMouseMode() == "selection": self.selectiontool.setDown(True)
		if self.scenegraph().getMouseMode() == "multiselection": self.multiselectiontool.setDown(True)
		if self.scenegraph().getMouseMode() == "rotate": self.rotatetool.setDown(True)
		if self.scenegraph().getMouseMode() == "xytranslate": self.translatetool.setDown(True)
		if self.scenegraph().getMouseMode() == "ztranslate": self.ztranslate.setDown(True)
		if self.scenegraph().getMouseMode() == "scale": self.scaletool.setDown(True)
		if self.scenegraph().getMouseMode() == "ruler": self.rulertool.setDown(True)
		if self.scenegraph().getMouseMode() == "cube": self.cubetool.setDown(True)
		if self.scenegraph().getMouseMode() == "sphere": self.spheretool.setDown(True)
		if self.scenegraph().getMouseMode() == "cylinder": self.cylindertool.setDown(True)
		if self.scenegraph().getMouseMode() == "cone": self.conetool.setDown(True)
		if self.scenegraph().getMouseMode() == "line": self.linetool.setDown(True)
		if self.scenegraph().getMouseMode() == "text": self.texttool.setDown(True)
		if self.scenegraph().getMouseMode() == "data": self.datatool.setDown(True)
		if self.scenegraph().getMouseMode() == "app": self.apptool.setDown(True)
		# Enable/Disable some tool buttons
		if self.scenegraph().getAPix():
			self.rulertool.setEnabled(True)
		else:
			self.rulertool.setEnabled(False)
			if self.scenegraph().getMouseMode() == "ruler":
				self.scenegraph().setMouseMode("rotate")	# Return to a default state
				self.rotatetool.setDown(True)
		# Lights
		if self.lighttab_open:
			position =  self.scenegraph().firstlight.getAngularPosition()
			self.lightwidget.setAngularPosition(position[0], position[1])
			angularposition = self.lightwidget.getAngularPosition()
			self.hvalslider.setValue(angularposition[0], quiet=1)
			self.vvalslider.setValue(angularposition[1], quiet=1)
			al = self.scenegraph().firstlight.getAmbient()
			ambient = old_div((al[0] + al[1] + al[2]),3)
			self.ambientlighting.setValue(ambient, quiet=1)
		# camera
		if self.cameratab_open:
			self.near.setValue(int(self.scenegraph().camera.getClipNear()), quiet=1)
			self.far.setValue(int(self.scenegraph().camera.getClipFar()), quiet=1)
			self.linkcb.setChecked(self.scenegraph().camera.getLinkingMode())
			self.capcb.setChecked(self.scenegraph().camera.getCappingMode())
			self.cappingcolor.setColor(QtGui.QColor(255*self.scenegraph().camera.getCapColor()[0],255*self.scenegraph().camera.getCapColor()[1],255*self.scenegraph().camera.getCapColor()[2]))
			self._get_vv_state()
			self.scenegraph().setZslice()
			self.camerawidget.updateWidget()
		# utils
		self.backgroundcolor.setColor(QtGui.QColor(255*self.scenegraph().clearcolor[0],255*self.scenegraph().clearcolor[1],255*self.scenegraph().clearcolor[2]))
		self.hideselectionbutton.setChecked(self.scenegraph().isSelectionHidded())
	
	def updateTree(self, currentnode=None):
		"""
		Update the SG tree
		"""
		self.clearTree()
		self.loadSG()
		# either set the current node to the argument or set it to the one PM has selected
		if currentnode:
			self.tree_widget.setCurrentItem(currentnode.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(currentnode.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			self.scenegraph().setCurrentSelection(currentnode)
		else:
			node = self.scenegraph().getCurrentSelection()
			self.tree_widget.setCurrentItem(node.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(node.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			
	def updateSceneGraph(self):
		""" 
		Updates SG, in the near future this will be improved to allow for slow operations
		"""
		self.scenegraph().updateSG()

class EMSGNodeInspector(EMItem3DInspector):
	"""
	Inspector for the SG node, allows special functionality for the SG node
	"""
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
		if name=="All Objects":
			tabwidget = QtWidgets.QWidget()
			gridbox = QtWidgets.QGridLayout()
			
			EMSGNodeInspector.addExtraTabAllObjects(self, gridbox)
			
			tabwidget.setLayout(gridbox)
			self.addTab(tabwidget, "misc")
		
	def addExtraTabAllObjects(self,gridbox):
		self.getthresh = QtWidgets.QLabel("Iso-threshold")
		self.isothr_box=QtWidgets.QLineEdit("0.0")
		gridbox.addWidget(self.getthresh, 1, 0, 1, 1)
		gridbox.addWidget(self.isothr_box, 1, 1, 1, 1)
		self.isothr_box.returnPressed.connect(self._on_change_threshold)
		#QtCore.QObject.connect(self.getthresh, QtCore.SIGNAL("clicked()"), self._on_get_thresh)
		
	def _on_change_threshold(self):
		try:
			value = float(self.isothr_box.text())
		except:
			return
		for c in  self.item3d().getChildren():
			for iso in c.getChildren():
				try:
					iso.setThreshold(value)
					#iso.updateSceneGraph()
				except: pass
		
	#def _on_get_thresh(self):
		##print value,
		#mx=[]
		#mn=[]
		#for c in  self.item3d().getChildren():
			#mx.append(c.data["maximum"])
			#mn.append(c.data["minimum"])
		#self.isothr_slider.setRange(min(mn),max(mx))

	def addControls(self, gridbox):
		super(EMSGNodeInspector, self).addControls(gridbox)
		buttonframe = QtWidgets.QFrame()
		buttonframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
		buttongrid = QtWidgets.QGridLayout()
		# Make buttons
		centerall = QtWidgets.QPushButton("Center All")
		distributeall =  QtWidgets.QPushButton("Distribute All")
		distributeall.setToolTip('Only works for data objects')
		selectall = QtWidgets.QPushButton("Select All")
		addaxes = QtWidgets.QPushButton("Add Axes")
		# Set layouts
		buttongrid.addWidget(centerall, 0, 0)
		buttongrid.addWidget(distributeall, 0, 1)
		buttongrid.addWidget(selectall, 1, 0)
		buttongrid.addWidget(addaxes, 1, 1)
		buttonframe.setLayout(buttongrid)
		gridbox.addWidget(buttonframe, 3, 0, 1, 1)
		# Add connections
		centerall.clicked.connect(self._on_centerall)
		distributeall.clicked.connect(self._on_distributeall)
		selectall.clicked.connect(self._on_selectall)
		addaxes.clicked.connect(self._on_addaxes)
		
	def _on_centerall(self):
		for child in self.item3d().getChildren():
			child.getTransform().set_trans(0.0, 0.0, 0.0)
		self.inspector().updateSceneGraph()
	
	def _on_distributeall(self):
		count = 0
		distcount = 0
		for child in self.item3d().getChildren():
			# Use modulo arith to distribute
			if hasattr(child, 'getBoundingBoxDimensions'):
				dims = child.getBoundingBoxDimensions()
				# distribute along X not matter what the SG root matrix
				self._set_transformSTDCorrd(child, math.pow(-1,(count%2))*distcount*dims[0], 0.0, 0.0)
				if (count + 1) % 2: distcount += 1
				count += 1
		self.inspector().updateSceneGraph()
				
	def _on_selectall(self):
		haschildren = False
		for child in self.item3d().getChildren():
			child.setSelectedItem(True)
			haschildren = True
		if haschildren:	
			self.item3d().setSelectedItem(False)
			self.inspector().updateTreeSelVis()
			self.inspector().updateSceneGraph()
	
	def _on_addaxes(self):
		# define axes size
		length = 50
		margin = 50
		linewidth = 20
		# define axes location
		xi = old_div(-self.item3d().camera.getViewPortWidthScaling()*self.item3d().camera.getWidth(),2) + margin*self.item3d().camera.getViewPortWidthScaling()
		yi = old_div(self.item3d().camera.getViewPortHeightScaling()*self.item3d().camera.getHeight(),2) - (length+margin)*self.item3d().camera.getViewPortHeightScaling()
		zi = 0
		#transform=Transform({'type':'eman','tx':xi,'ty':yi,'tz':zi,'scale':self.item3d().camera.getViewPortWidthScaling()})
		# create axes
		xaxis = EMLine(0,0,0,length,0,0,linewidth)
		xaxis.setShowLeftArrow(False)
		xaxis.setRightArrowLength(old_div(length,5))
		xaxis.setAmbientColor(1.0,0.0,0.0)
		yaxis = EMLine(0,0,0,0,length,0,linewidth)
		yaxis.setShowLeftArrow(False)
		yaxis.setRightArrowLength(old_div(length,5))
		yaxis.setAmbientColor(0.0,1.0,0.0)
		zaxis = EMLine(0,0,0,0,0,length,linewidth)
		zaxis.setShowLeftArrow(False)
		zaxis.setRightArrowLength(old_div(length,5))
		zaxis.setAmbientColor(0.0,0.0,1.0)
		# Make the axes root node
		axesnode = EMItem3D()
		self._set_transformSTDCorrd(axesnode, xi, yi, zi)
		axesnode.getTransform().set_scale(self.item3d().camera.getViewPortWidthScaling())
		self.item3d().insertNewNode("Axes", axesnode, self.item3d())
		# added the axes to root node
		self.item3d().insertNewNode("X Axis (red)", xaxis, axesnode)
		self.item3d().insertNewNode("Y Axis (green)", yaxis, axesnode)
		self.item3d().insertNewNode("Z Axis (blue)", zaxis, axesnode)
		self.inspector().updateTree()
		self.inspector().updateSceneGraph()
	
	def _set_transformSTDCorrd(self, child, tx, ty, tz):
		""" Helper function to set transform in standard coord system """
		tt = t = Transform({"tx":tx,"ty":ty,"tz":tz})
		tp = child.getParentMatrixProduct()
		if tp: tt = tp.inverse()*t
		child.getTransform().set_trans(tt.get_trans())
		
class EMQTreeWidget(QtWidgets.QTreeWidget):
	"""
	Subclassing the QTreeWidget to enable is_visible toggling
	"""
	visibleItem = QtCore.pyqtSignal(QtWidgets.QTreeWidgetItem)
	editItem = QtCore.pyqtSignal(QtWidgets.QTreeWidgetItem)

	def __init__(self, parent=None):
		QtWidgets.QTreeWidget.__init__(self, parent)
			
	def mousePressEvent(self, e):
		QtWidgets.QTreeWidget.mousePressEvent(self, e)
		if e.button()==Qt.RightButton:
			self.visibleItem.emit(self.currentItem())
		if e.button()==Qt.MidButton or (e.buttons()&Qt.LeftButton and e.modifiers()&Qt.AltModifier):
			self.editItem.emit(self.currentItem())
			
			
class EMQTreeWidgetItem(QtWidgets.QTreeWidgetItem):
	"""
	Subclass of QTreeWidgetItem
	adds functionality
	"""

	def __init__(self, qstring, item3d, parentnode):
		QtWidgets.QTreeWidgetItem.__init__(self, [qstring])
		self.name = ''.join(qstring)
		self.item3d = weakref.ref(item3d)
		if parentnode: self.parent = weakref.ref(parentnode)
		else: self.parent = None
		self.visible = QtGui.QIcon(QtGui.QPixmap(visibleicon))
		self.invisible = QtGui.QIcon(QtGui.QPixmap(invisibleicon))
		self.setCheckState(0, QtCore.Qt.Unchecked)
		self.getVisibleState()
		self.setToolTip(0, 'Click on the checkbox to select\nMiddle click to edit\nRight click to toggle visible')
	
	def setSelectionState(self, state):
		""" 
		Toggle selection state on and off
		"""
		if state == QtCore.Qt.Checked:
			self.item3d().setSelectedItem(True)
		else:
			self.item3d().setSelectedItem(False)
		self.setSelectionStateBox() # set state of TreeItemwidget
		
	def toggleVisibleState(self):
		""" 
		Toggle visible state on and off
		"""
		self.item3d().setVisibleItem(not self.item3d().isVisibleItem())
		self.getVisibleState()
		
	def getVisibleState(self):
		"""
		Set the visible icon to the item state
		"""
		if self.item3d().isVisibleItem():
			self.setIcon(0, self.visible)
		else:
			self.setIcon(0, self.invisible)
	
	def setSelectionStateBox(self):
		"""
		Set the selection state icon
		"""
		if self.item3d().isSelectedItem():
			self.setCheckState(0, QtCore.Qt.Checked)
		else:
			self.setCheckState(0, QtCore.Qt.Unchecked)
	
	def removeAllChildren(self, inspector):
		"""
		Remove all children from the SG
		"""
		for i in range(self.childCount()):
			self.child(0).removeAllChildren(inspector)
			inspector.removeTreeNode(self, 0) 

class NodeEditDialog(QtWidgets.QDialog):
	"""
	A dialog for editing the node
	"""
	def __init__(self, inspector, item):
		QtWidgets.QDialog.__init__(self)
		self.item = item
		self.inspector = weakref.ref(inspector)
		grid = QtWidgets.QGridLayout(self)
		label = QtWidgets.QLabel("Node Name")
		self.nodename = QtWidgets.QLineEdit(self.item.item3d().getLabel())
		grid.addWidget(label, 0, 0, 1, 2)
		grid.addWidget(self.nodename, 1, 0, 1, 2)
		self.ok_button = QtWidgets.QPushButton("Ok")
		self.cancel_button = QtWidgets.QPushButton("Cancel")
		grid.addWidget(self.ok_button, 2, 0, 1, 1)
		grid.addWidget(self.cancel_button, 2, 1, 1, 1)
		self.setLayout(grid)
		
		self.ok_button.clicked.connect(self._on_ok)
		self.cancel_button.clicked.connect(self._on_cancel)
	
	def _on_ok(self):
		self.item.item3d().setLabel(self.nodename.text())
		self.item.setText(0, self.nodename.text())
		self.item.item3d().setLabel(self.nodename.text())
		self.done(0)
		
	def _on_cancel(self):
		self.done(1)
	
class NodeDialog(QtWidgets.QDialog):
	"""
	Generate a dialog to add or remove node. If remove is chosen 'item' node is removed
	If add node is chosen, a node is inserted just below this node.
	"""
	def __init__(self, inspector, item):
		QtWidgets.QDialog.__init__(self)
		self.item = item
		self.inspector = weakref.ref(inspector)
		self.setWindowTitle('Node Controller')
		self.setMaximumWidth(300)
		self.transformgroup = {}
		vbox = QtWidgets.QVBoxLayout(self)
		# Stuff within the frame
		frame = QtWidgets.QFrame()
		frame.setFrameStyle(QtWidgets.QFrame.StyledPanel)
		fvbox = QtWidgets.QVBoxLayout(frame)
		label = QtWidgets.QLabel("Node Type to add")
		self.node_type_combo = QtWidgets.QComboBox() 
		self.node_stacked_widget = QtWidgets.QStackedWidget()
		self.node_stacked_widget.setFrameStyle(QtWidgets.QFrame.StyledPanel)
		self.addnode_button = QtWidgets.QPushButton("Add Node")
		fvbox.addWidget(label)
		fvbox.addWidget(self.node_type_combo)
		fvbox.addWidget(self.node_stacked_widget)
		fvbox.addWidget(self.addnode_button)
		frame.setLayout(fvbox)
		# vbox widgets
		vbox.addWidget(frame)
		self.cancel_button = QtWidgets.QPushButton("Cancel")
		vbox.addWidget(self.cancel_button)
		self.setLayout(vbox)
		
		# Populate node types		
		self.node_type_combo.addItem("Cube")
		self.cubewidgetdict = {}
		self.node_stacked_widget.addWidget(EMCube.getNodeDialogWidget(self.cubewidgetdict))
		self.node_type_combo.addItem("Sphere")
		self.spherewidgetdict = {}
		self.node_stacked_widget.addWidget(EMSphere.getNodeDialogWidget(self.spherewidgetdict))
		self.node_type_combo.addItem("Cylinder")
		self.cylinderwidgetdict = {}
		self.node_stacked_widget.addWidget(EMCylinder.getNodeDialogWidget(self.cylinderwidgetdict))
		self.node_type_combo.addItem("Cone")
		self.conewidgetdict = {}
		self.node_stacked_widget.addWidget(EMCone.getNodeDialogWidget(self.conewidgetdict))
		self.node_type_combo.addItem("Line")
		self.linewidgetdict = {}
		self.node_stacked_widget.addWidget(EMLine.getNodeDialogWidget(self.linewidgetdict))
		self.node_type_combo.addItem("Text")
		self.textwidgetdict = {}
		self.node_stacked_widget.addWidget(EM3DText.getNodeDialogWidget(self.textwidgetdict))
		
		self.node_type_combo.addItem("PDB")
		self.pdbwidgetdict = {}
		self.node_stacked_widget.addWidget(EMPDBItem3D.getNodeDialogWidget(self.pdbwidgetdict))
		if self.item and self.item.item3d().name == "PDB":
			self.node_type_combo.addItem("Ball and Stick")
			self.ballstickwidgetdict = {}
			self.node_stacked_widget.addWidget(EMBallStickModel.getNodeDialogWidget(self.ballstickwidgetdict))
			self.node_type_combo.addItem("Spheres")
			self.sphereswidgetdict = {}
			self.node_stacked_widget.addWidget(EMSphereModel.getNodeDialogWidget(self.sphereswidgetdict))
		
		self.node_type_combo.addItem("Data")
		self.datawidgetdict = {}
		self.node_stacked_widget.addWidget(EMDataItem3D.getNodeDialogWidget(self.datawidgetdict))
		if self.item and self.item.item3d().name == "Data":
			self.node_type_combo.addItem("Isosurface")
			self.isosurfacewidgetdict = {}
			self.node_stacked_widget.addWidget(EMIsosurface.getNodeDialogWidget(self.isosurfacewidgetdict))
			self.node_type_combo.addItem("Slice")
			self.slicewidgetdict = {}
			self.node_stacked_widget.addWidget(EMSliceItem3D.getNodeDialogWidget(self.slicewidgetdict))
			self.node_type_combo.addItem("Volume")
			self.volumewidgetdict = {}
			self.node_stacked_widget.addWidget(EMVolumeItem3D.getNodeDialogWidget(self.volumewidgetdict))
		
		self.addnode_button.clicked.connect(self._on_add_node)
		self.cancel_button.clicked.connect(self._on_cancel)
		self.node_type_combo.activated[int].connect(self._node_combobox_changed)
	
	def _on_add_node(self):
		insertion_node = None
		node_name = "default" 
		parentnode = None
		d = None
		reverttrans = False
		# Cube
		if self.node_type_combo.currentText() == "Cube":
			insertion_node = EMCube.getNodeForDialog(self.cubewidgetdict)
			node_name = str(self.cubewidgetdict["node_name"].text())
		# Sphere
		if self.node_type_combo.currentText() == "Sphere":
			insertion_node = EMSphere.getNodeForDialog(self.spherewidgetdict)
			node_name = str(self.spherewidgetdict["node_name"].text())
		# Cylinder
		if self.node_type_combo.currentText() == "Cylinder":
			insertion_node = EMCylinder.getNodeForDialog(self.cylinderwidgetdict)
			node_name = str(self.cylinderwidgetdict["node_name"].text())
		# Cone
		if self.node_type_combo.currentText() == "Cone":
			insertion_node = EMCone.getNodeForDialog(self.conewidgetdict)
			node_name = str(self.conewidgetdict["node_name"].text())
		# Line
		if self.node_type_combo.currentText() == "Line":
			insertion_node = EMLine.getNodeForDialog(self.linewidgetdict)
			node_name = str(self.linewidgetdict["node_name"].text())
		# Text
		if self.node_type_combo.currentText() == "Text":
			insertion_node = EM3DText.getNodeForDialog(self.textwidgetdict)
			node_name = str(self.textwidgetdict["node_name"].text())
		# PDB
		if self.node_type_combo.currentText() == "PDB":
			insertion_node = EMPDBItem3D.getNodeForDialog(self.pdbwidgetdict)
			node_name = str(self.pdbwidgetdict["node_name"].text())
		# Ball/Stick
		if self.node_type_combo.currentText() == "Ball and Stick": 
			self.ballstickwidgetdict["parent"] = parentnode = self.item.item3d()
			insertion_node = EMBallStickModel.getNodeForDialog(self.ballstickwidgetdict)
			node_name = str(self.ballstickwidgetdict["node_name"].text())
		# Spheres
		if self.node_type_combo.currentText() == "Spheres": 
			self.sphereswidgetdict["parent"] = parentnode = self.item.item3d()
			insertion_node = EMSphereModel.getNodeForDialog(self.sphereswidgetdict)
			node_name = str(self.sphereswidgetdict["node_name"].text())
		# Data
		if self.node_type_combo.currentText() == "Data": 
			insertion_node = EMDataItem3D.getNodeForDialog(self.datawidgetdict)
			node_name = str(self.datawidgetdict["node_name"].text())
		# Isosurface
		if self.node_type_combo.currentText() == "Isosurface": 
			self.isosurfacewidgetdict["parent"] = parentnode = self.item.item3d()
			insertion_node = EMIsosurface.getNodeForDialog(self.isosurfacewidgetdict)
			node_name = str(self.isosurfacewidgetdict["node_name"].text())
			reverttrans = True
		# Slice
		if self.node_type_combo.currentText() == "Slice": 
			self.slicewidgetdict["parent"] = parentnode = self.item.item3d()
			insertion_node = EMSliceItem3D.getNodeForDialog(self.slicewidgetdict)
			node_name = str(self.slicewidgetdict["node_name"].text())
			reverttrans = True
		# Volume
		if self.node_type_combo.currentText() == "Volume": 
			self.volumewidgetdict["parent"] = parentnode = self.item.item3d()
			insertion_node = EMVolumeItem3D.getNodeForDialog(self.volumewidgetdict)
			node_name = str(self.volumewidgetdict["node_name"].text())
			reverttrans = True
		
		insertion_node.setLabel(node_name)
		self.inspector().scenegraph().insertNewNode(node_name, insertion_node, parentnode=parentnode)
		insertion_node.setTransform(insertion_node.getParentMatrixProduct().inverse()*insertion_node.getTransform()) # Move to standard coordinate system
		insertion_node.getTransform().set_scale(1.0)	# The scale can be adversely affected by the above line of code. This may or may not be optimal I'll have to think about it....
		#if reverttrans: insertion_node.getTransform().set_trans(float(d[0].text()),float(d[1].text()),float(d[2].text()))
		self.inspector().updateSceneGraph()
		self.done(0)
		
	def _on_cancel(self):
		self.done(1)
	
	def _node_combobox_changed(self, nodetype):
		self.node_stacked_widget.setCurrentIndex(nodetype)
		
###################################### TEST CODE, THIS WILL NOT APPEAR IN THE WIDGET3D MODULE ##################################################
		
# All object that are rendered inherit from abstractSGnode and implement the render method

class GLdemo(QtWidgets.QWidget):
	def __init__(self):
		QtWidgets.QWidget.__init__(self)
		
		self.widget = EMScene3D()
		#self.widget.camera.usePerspective(50, 0.5)
		self.cube = EMCube(50.0)
		self.widget.addChild(self.cube)    # Something to Render something..... (this could just as well be one of Ross's SGnodes)
		self.sphere = EMSphere(50.0)
		self.widget.addChild(self.sphere)
		self.cylider = EMCylinder(50.0, 50.0)
		self.widget.addChild(self.cylider)
		
		# QT stuff to display the widget
		vbox = QtWidgets.QVBoxLayout()
		vbox.addWidget(self.widget)
		self.setLayout(vbox)
		self.setGeometry(300, 300, 600, 600)
		self.setWindowTitle('BCM EM Viewer')
		
	def closeEvent(self, event): 
		self.widget.close()


def main():
	import sys
	#from pmwidgets import PMIntEntryWidget, PMStringEntryWidget, PMBoolWidget, PMFileNameWidget, PMComboWidget
	app = QtWidgets.QApplication(sys.argv)
	window = GLdemo()
	window.show()
	app.exec_()


if __name__ == "__main__":
	main()

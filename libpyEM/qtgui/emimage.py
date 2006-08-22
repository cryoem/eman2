#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU
#from valslider import ValSlider
#from math import *
#from EMAN2 import *
#import sys
#import Numeric

class EMImage(object):
	"""This is basically a factory class that will return an instance of the appropriate EMImage* class """

	def __new__(cls,*args,*kw):
		
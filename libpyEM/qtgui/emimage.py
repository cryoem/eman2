#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emimage2d import *
from emimagemx import *
from emimage3d import *
#from OpenGL import GL,GLU
#from valslider import ValSlider
#from math import *
#from EMAN2 import *
#import sys
#import Numeric

class EMImage(object):
	"""This is basically a factory class that will return an instance of the appropriate EMImage* class """

	def __new__(cls,data=None,old=None):
		"""This will create a new EMImage* object depending on the type of 'data'. If
		old= is provided, and of the appropriate type, it will be used rather than creating
		a new instance."""
		
		if isinstance(data,EMData) and data.get_zsize()==1:
			# single 2D image
			if old:
				if isinstance(old,EMImage2D) :
					old.setData(data)
					return old
			return EMImage2D(data)
		elif isinstance(data,EMData):
			# must be a single 3D image
			if old:
				if isinstance(old,EMImage3D) :
					old.setData(data)
					return old
			return EMImage3D(data)
		elif isinstance(data,list):
			# list or tuple of images
			if old:
				if isinstance(old,EMImageMX) :
					old.setData(data)
					return old
			return EMImageMX(data)
		else:
			raise Exception,"data must be a single EMData object or a list of EMData objects"
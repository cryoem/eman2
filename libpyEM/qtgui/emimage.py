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

def get_app():
	app=QtGui.QApplication.instance()
	if not app : app = QtGui.QApplication([])
	
	try: 
		if app.updtimer : pass
	except:
		tmr=QtCore.QTimer()
		tmr.setInterval(250)
		tmr.connect(tmr,QtCore.SIGNAL("timeout()"), imageupdate)
		tmr.start()
	
		app.updtimer=tmr

	return app
		
def imageupdate():
	for i in EMImage2D.allim.keys():
		if i.data.get_attr("changecount")!=i.changec :
			i.setData(i.data)

	for i in EMImage3D.allim.keys():
		if i.data.get_attr("changecount")!=i.changec :
			i.setData(i.data)
	
	for i in EMImageMX.allim.keys():
		try:
			if len(i.data)!=i.nimg : i.setData(i.data)
		except:
			pass
		for j in i.changec.keys():
			upd=0
			if j.get_attr("changecount")!=i.changec[j] :
				upd=1
				break
		if upd : i.setData(i.data)


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
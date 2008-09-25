#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emimage2d import *
from emimagemx import *
from emimagemxrotor import *
from emimage3d import *
from emimageutil import EMParentWin
from OpenGL import GL,GLU,GLUT
from sys import argv
from copy import deepcopy
#from valslider import ValSlider
#from math import *
#from EMAN2 import *
#import sys
#import numpy

GLUT.glutInit(sys.argv )

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
		try:
			if i.data.get_attr("changecount")!=i.changec :
				i.set_data(i.data)
		except: pass

	for i in EMImage3D.allim.keys():
		try:
			if i.data.get_attr("changecount")!=i.changec :
				i.set_data(i.data)
		except: pass
		
	for i in EMImageMX.allim.keys():
		try:
			if len(i.data)!=i.nimg : i.set_data(i.data)
		except:
			pass
		upd=0
		try:
			for j in i.changec.keys():
				try:
					if j.get_attr("changecount")!=i.changec[j] :
						upd=1
						break
				except: pass
		except: pass
		if upd : i.set_data(i.data)

class EMImage(object):
	"""This is basically a factory class that will return an instance of the appropriate EMImage* class """
	def __new__(cls,data=None,old=None,parent=1,copy=True):
		"""This will create a new EMImage* object depending on the type of 'data'. If
		old= is provided, and of the appropriate type, it will be used rather than creating
		a new instance."""
		if isinstance(data,EMData) and data.get_zsize()==1:
			# single 2D image
			
			# sometimes it's necessary to copy, especially if the user is 
			# calling display from python
			if copy: local_data = data.copy()
			else: local_data = data

			if old:
				if isinstance(old,EMImage2D) :
					old.set_data(local_data)
					return old
			if parent: 
				ret=EMParentWin(EMImage2D(local_data))
			if parent : 
				ret=EMParentWin(EMImage2D(local_data))
				ret.show()
#				ret.raise()
				ret.releaseMouse()
				ret.releaseKeyboard()
				return ret
			return EMImage2D(local_data)
		elif isinstance(data,EMData):
			# data copy considerations shouldn't be necessary here 
			# seeing as the EMImage3D does internal copying of its own
			# FIXME double check this aspect of the code once things are going
			# must be a single 3D image
			if old:
				if isinstance(old,EMImage3D) :
					old.set_data(data)
					return old
			if parent : 
				ret=EMParentWin(EMImage3D(data))
				ret.show()
#				ret.raise()
				ret.releaseMouse()
				ret.releaseKeyboard()
				return ret
			return EMImage3D(data)
		elif isinstance(data,list):
			
			if ( copy ):local_data = deepcopy(data)
			else: local_data = data
			
			if old:
				if isinstance(old,EMImageMX) :
					old.set_data(local_data)
					return old
			if parent : 
				ret=EMParentWin(EMImageMX(local_data))
				ret.show()
#				ret.raise()
				ret.releaseMouse()
				ret.releaseKeyboard()
				return ret
			return EMImageMX(local_data)
		else:
			raise Exception,"data must be a single EMData object or a list of EMData objects"

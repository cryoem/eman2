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
from OpenGL import GL,GLU,GLUT
from EMAN2 import Util,EMUtil,file_exists,IMAGE_UNKNOWN,gimme_image_dimensions3D,EMData
import os
from EMAN2 import Transform
from emscene3d import EMScene3D
from emdataitem3d import EMDataItem3D, EMIsosurface

def image_update():
	from emimage2d import EMImage2DWidget
	from emimagemx import EMImageMXWidget
	from emimage3d import EMImage3DWidget
	for i in EMImage2DWidget.allim.keys():
		try:
			if i.isVisible() and i.data["changecount"] !=i.image_change_count:
				i.force_fft_redo()
				i.force_display_update()
				i.updateGL()
		except: pass
	
	for i in EMImageMXWidget.allim.keys():
		try:
			if i.isVisible() and i.data[0]["changecount"]!=i.image_change_count:
				i.force_display_update()
				i.updateGL()
		except: pass
		
	for i in EMImage3DWidget.allim.keys():
		try:
			if i.isVisible() and i.data["changecount"]!=i.image_change_count:
				i.updateGL()
		except: pass
	
	
def get_app():
	'''
	Deprecated
	But being replaced by emapplication.get_application (in progress, April 15th 2009)
	'''
	app=QtGui.QApplication.instance()
	if not app : app = QtGui.QApplication([])
	
	try: 
		if app.updtimer : pass
	except:
		tmr=QtCore.QTimer()
		tmr.setInterval(250)
		tmr.connect(tmr,QtCore.SIGNAL("timeout()"), image_update)
		tmr.start()
	
		app.updtimer=tmr

	return app


class EMImageWidget(object):
	"""This is basically a factory class that will return an instance of the appropriate EMImage* class """
	def __new__(cls,data=None,old=None,app=None,force_2d=False,force_plot=False,filename="",replace=True):
		"""This will create a new EMImage* object depending on the type of 'data'. If
		old= is provided, and of the appropriate type, it will be used rather than creating
		a new instance.
		"""
		
		if isinstance(data,EMData) and data.get_size()==0: raise RuntimeError("Can not display an EMData object that has no pixels")
		
		from EMAN2 import remove_directories_from_name
		if force_plot and force_2d:
			# ok this sucks but it suffices for the time being
			print "Error, the force_plot and force_2d options are mutually exclusive"
			return None
		
		if force_plot or (isinstance(data,EMData) and data.get_zsize()==1 and data.get_ysize()==1):
			from emplot2d import EMPlot2DWidget
			if old:
				if isinstance(old,EMPlot2DWidget) :
					old.set_data(data,remove_directories_from_name(filename),replace)
					return old
			widget = EMPlot2DWidget(application=app)
			widget.set_data(data,remove_directories_from_name(filename),replace)
			return widget	
		elif force_2d or (isinstance(data,EMData) and data.get_zsize()==1):
			from emimage2d import EMImage2DWidget
			if old:
				if isinstance(old,EMImage2DWidget) :
					old.set_data(data,filename)
					return old
			widget = EMImage2DWidget(application=app)
			widget.set_data(data,filename)
			return widget
		elif isinstance(data,EMData):
			if isinstance(old,EMScene3D): widget = old
			else: widget = EMScene3D()
			data = EMDataItem3D(data, transform=Transform())
			#data.setSelectedItem(True)
			isosurface = EMIsosurface(data, transform=Transform())
			widget.insertNewNode(os.path.basename(filename), data, parentnode=widget)
			widget.insertNewNode("Iso", isosurface, parentnode=data)
			return widget

		elif isinstance(data,list) and isinstance(data[0],EMData):
			from emimagemx import EMImageMXWidget
			if old:
				if isinstance(old,EMImageMXWidget) :
					old.set_data(data,filename)
					return old
			widget = EMImageMXWidget(application=app)
			widget.set_data(data,filename)
			return widget
		elif isinstance(data,list):
			from emplot3d import EMPlot3DWidgetNew
			if (isinstance(data[0],list) or isinstance(data[0],tuple)) and len(data) > 2:
				if old:
					if isinstance(old,EMPlot3DWidgetNew) :
						old.set_data(data,remove_directories_from_name(filename),replace)
						return old
						
				widget = EMPlot3DWidgetNew()
				widget.set_data(data,remove_directories_from_name(filename),replace)
				return widget	
			else:
				from emplot2d import EMPlot2DWidget
				if old:
					if isinstance(old,EMPlot2DWidget) :
						old.set_data(data,remove_directories_from_name(filename),replace)
						return old
				widget = EMPlot2DWidget(application=app)
				widget.set_data(data,remove_directories_from_name(filename),replace)
				return widget	
		else:
			raise Exception,"data must be a single EMData object or a list of EMData objects"


class EMWidgetFromFile(object):
	"""This is basically a factory class that will return an instance of the appropriate EMDisplay class,
	using only a file name as input. Can force plot and force 2d display, also.
	
	Used by emselector.py and e2display.py.
	
	This object was retrospectively altered to allow a file name to be passed into the EMImageMXWidget's set_data function,
	this is to facilitate viewing large images using the EMImageMXWidget's caching mechanism. The object reads the images
	from disk internally, allowing the user to view very large sets.
	
	"""
	def __new__(cls,filename,application,force_plot=False,force_2d=False,old=None):
		
		file_type = Util.get_filename_ext(filename)
		em_file_type = EMUtil.get_image_ext_type(file_type)
		if not file_exists(filename): return None
		
		if force_plot and force_2d:
			# ok this sucks but it suffices for the time being
			print "Error, the force_plot and force_2d options are mutually exclusive"
			return None
		
		if force_plot:
			from emplot2d import EMPlot2DWidget
			if isinstance(old,EMPlot2DWidget): widget = old
			else: widget = EMPlot2DWidget(application=application)
			widget.set_data_from_file(filename)
			return widget
		
		if em_file_type != IMAGE_UNKNOWN or filename[:4] == "bdb:":
			n = EMUtil.get_image_count(filename)
			nx,ny,nz = gimme_image_dimensions3D(filename)
			if n > 1 and nz == 1: 
				if force_2d:
					a = EMData()
					data=a.read_images(filename)
				else:
					data = None # This is like a flag - the ImageMXWidget only needs the file name
			elif nz == 1:
				data = [EMData(filename,0)]
			else:
				data = EMData()
				data.read_image(filename,0,not force_2d)		# This should be 3-D. We read the header-only here
				data = [data]
				
			if data != None and len(data) == 1: data = data[0]
			
			if force_2d or isinstance(data,EMData) and data.get_zsize()==1:
				if isinstance(data,list) or data.get_ysize() != 1:
					from emimage2d import EMImage2DWidget
					if isinstance(old,EMImage2DWidget): widget = old
					else: widget= EMImage2DWidget(application=application)
				else:
					from emplot2d import EMPlot2DWidget
					if isinstance(old,EMPlot2DWidget): widget = old
					else: widget = EMPlot2DWidget(application=application)
					widget.set_data_from_file(filename)
					return widget
			elif isinstance(data,EMData):
				if isinstance(old,EMScene3D): widget = old
				else: widget = EMScene3D()
#				print n,data
				for ii in xrange(n):
					data=EMData(filename,ii)
					datai = EMDataItem3D(data, transform=Transform())
					widget.insertNewNode(os.path.basename(filename), datai, parentnode=widget)
					isosurface = EMIsosurface(datai, transform=Transform())
					widget.insertNewNode("Iso", isosurface, parentnode=datai)
				return widget
				
			elif data == None or isinstance(data,list):
				from emimagemx import EMImageMXWidget
				if isinstance(old,EMImageMXWidget): widget = old
				else: widget = EMImageMXWidget(application=application)
				data = filename
			else: 
				print filename
				raise # weirdness, this should never happen
			widget.set_data(data,filename)
			return widget
		else:
			from emplot2d import EMPlot2DWidget
			if isinstance(old,EMPlot2DWidget): widget = old
			else: widget = EMPlot2DWidget(application=application)
			widget.set_data_from_file(filename)
			return widget


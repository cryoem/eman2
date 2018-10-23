#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
ploticon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccbbbbccccc',
    'ccccbbccccbbccc',
    'cccbccccccccbcc',
    'ccbccccccccccbc',
    'ccbccccccccccbc',
    'ccbccccccccccbc',
    'ccbccccccccccbc',
    'cccbccccccccbcc',
    'ccccbbccccbbccc',
    'ccccccbbbbccccc',
    'ccccccccccccccc',
    'ccccccccccccccc'
]

def safe_float(x):
	try: return float(x)
	except: return 0.0

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU
from OpenGL.GL import *
import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo
from math import *
from EMAN2 import *
import sys
from .emshape import *
import weakref
from pickle import dumps,loads
import struct, math
from numpy import *
from .valslider import *
from io import StringIO
import re
from . import emimage2d

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
#matplotlib.use('Agg')
import numpy as np

from .emapplication import EMApp, EMGLWidget
from .emglobjects import EMOpenGLFlagsAndTools

import traceback

#plt.style.use('ggplot')

linetypes=["-","--",":","-."]
symtypes=["o","s","+","2","1"]
colortypes=["k","b","r","g","y","c","m","gray"]
qt_color_map = {}
qt_color_map["k"] = QtGui.QColor(0,0,0)
qt_color_map["b"] = QtGui.QColor(0,0,255)
qt_color_map["r"] = QtGui.QColor(255,0,0)
qt_color_map["g"] = QtGui.QColor(0,255,0)
qt_color_map["y"] = QtGui.QColor(255,255,0)
qt_color_map["c"] = QtGui.QColor(0,255,255)
qt_color_map["m"] = QtGui.QColor(255,0,255)
qt_color_map["gray"] = QtGui.QColor(127,127,127)

class EMPlot2DWidget(EMGLWidget):
	"""A QT widget for drawing 2-D plots using matplotlib
	"""
	selected_sg = QtCore.pyqtSignal()
	mousedown = QtCore.pyqtSignal(QtGui.QMouseEvent,tuple)

	def __init__(self,application=None,winid=None,parent=None):

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		EMGLWidget.__init__(self, parent=parent, winid=winid)
		self.setFormat(fmt)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"plot.png"))

		self.axes={}
		self.pparm={}			# color,line,linetype,linewidth,sym,symtype,symsize
		self.inspector=None
		self.needupd=1
		self.plotimg=None
		self.shapes={}
		self.xlimits=None
		self.ylimits=None
		self.climits=None
		self.slimits=None
		self.rmousedrag=None
		self.axisparms=(None,None,"linear","linear")
		self.selected=[]
		self.comments={}			# IF reading from a file which contains per-point comments, this dictionary contains a list of comments for each point

		self.data={}				# List of Lists to plot
		self.visibility = {}  	   	# Same entries as in self.data, but entries are true or False to indicate visibility
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags

		self.tex_name = 0
		self.main_display_list = 0

		self.resize(640,480)

		self.particle_viewers = []

		self.alpha = 0.5

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
		GL.glEnable(GL_DEPTH_TEST)

	def paintGL(self):
		try: GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
		except: pass # this is a hack.

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		self.render()

	def resizeGL(self, width, height):
		#print "resize ",self.width(), self.height()
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())

		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GL.glOrtho(0.0,self.width(),0.0,self.height(),-10,10)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()

		self.resize_event(width,height)

	def closeEvent(self,event):
		for pv in self.particle_viewers: 
			if pv!=None : pv.closeEvent(event)
		self.particle_viewers=[]

		self.clear_gl_memory()
		EMGLWidget.closeEvent(self, event)

		if self.inspector :
			self.inspector.closeEvent(event)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_C:
			self.show_inspector(1)
		elif event.key() == Qt.Key_F1:
			try: from PyQt4 import QtWebKit
			except: return
			try:
				try:
					test = self.browser
				except:
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2display"))
					self.browser.resize(800,800)

				if not self.browser.isVisible(): self.browser.show()
			except: pass

	def setWindowTitle(self,filename):
		EMGLWidget.setWindowTitle(self, remove_directories_from_name(filename,1))

	def clear_gl_memory(self):
		if self.tex_name != 0:
			GL.glDeleteTextures(self.tex_name)
			self.tex_name = 0
		if self.main_display_list != 0:
			glDeleteLists(self.main_display_list,1)
			self.main_display_list = 0

	def set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,linewidth=1,linetype=-2,symtype=-2,symsize=10,comments=None):
		"""Set a keyed data set. The key should generally be a string describing the data.
		'data' is a tuple/list of tuples/list representing all values for a particular
		axis. eg - the points: 1,5; 2,7; 3,9 would be represented as ((1,2,3),(5,7,9)).
		Multiple axes may be set, and which axis represents which axis in the plot can be
		selected by the user. 'data' can also be an EMData object, in which case the entire
		data array is plotted as if it were 1-D.

		linetype and symtype are integers. -1 will disable either. Default is autoselect."""

		#print "set_data ",key
		#traceback.print_stack()

		self.del_shapes()

		self.needupd=1

		if replace:
			self.data = {}
			self.axes = {}
			self.visibility = {}

		if input_data is None :
			self.data.pop(key)
			self.visibility.pop(key)
			self.axes.pop(key)
			try: self.comments.pop(key)
			except: pass
			if self.inspector: self.inspector.datachange()
			if not quiet: self.updateGL()
			return

		if key in self.data : oldkey=True
		else: oldkey=False

		if isinstance(input_data,EMData):
			data = input_data.get_data_as_vector()
		else: data = input_data

		if not isinstance(data[0],list) and not isinstance(data[0],tuple) and not isinstance(data[0],ndarray):
			x_axis = arange(len(data))
			data = [ x_axis,array(data) ]		# replace data with our refactored version
			self.data[key]= data
			self.visibility.setdefault(key,True)
		else:
			self.data[key]=[array(i) for i in data]
			self.visibility.setdefault(key,True)

		try:
			if len(data)>1 :
				if len(data)>2:
					try:
						if data[0][2]-data[0][1]==1 : self.axes[key]=(1,2,-2,-2)	# if it looks like the first axis is a boring count
						else : self.axes[key]=(0,1,-2,-2)
					except: self.axes[key]=(0,1,-2,-2)
				else : self.axes[key]=(0,1,-2,-2)
			else : self.axes[key]=(-1,0,-2,-2)
		except:
			print("Data error:", data)
			return

		if oldkey:
			pp=self.pparm[key]
			if color<0 : color=pp[0]
			if linetype==-2: doline=pp[1]
			elif linetype==-1: doline=0
			else: doline=1
			if symtype==-2: dosym=pp[4]
			elif symtype==-1: dosym=0
			else: dosym=1
			if symtype<0: symtype=pp[5]
			if linewidth<0: linewidth=pp[3]
			if symsize<0 : symsize=pp[6]
		else:
			if symtype==-2 and linetype==-2:
				if len(data)<4 and (diff(self.data[key][0])>=0).all() : doline,linetype=1,0
				else : dosym,symtype=1,0
			if color<0 : color=len(self.data)%len(colortypes)			# Automatic color setting
			if color >len(colortypes): color = 0 # there are only a certain number of colors
			if linetype>=0 : doline=1
			else : doline,linetype=0,0
			if symtype>=0 : dosym=1
			else : dosym,symtype=0,0
		self.pparm[key]=(color,doline,linetype,linewidth,dosym,symtype,symsize)

		if comments!=None:
			self.comments[key]=comments

		self.autoscale()

		if self.inspector: self.inspector.datachange()

		if not quiet: self.updateGL()

	def get_inspector(self):
		if not self.inspector :
			self.inspector=EMPlot2DInspector(self)
			self.inspector.datachange()
		return self.inspector

	def set_data_from_file(self,filename,replace=False,quiet=False):
		"""Reads a keyed data set from a file. Automatically interpret the file contents."""
		self.del_shapes()

		if replace:
			self.data = {}
			self.axes = {}
			self.visibility = {}

		file_type = Util.get_filename_ext(filename)
		em_file_type = EMUtil.get_image_ext_type(file_type)

		if em_file_type != IMAGE_UNKNOWN or filename[:4]=="bdb:" :

			im=EMData.read_images(filename)
			if len(im) == 1:
				im = im[0]
				l = [i for i in range(im.get_size())]
				k = im.get_data_as_vector()
				if filename in self.data : filename="{}.{}".format(filename,len(self.data))
				self.set_data([l,k],filename,quiet=quiet)
			elif im[0].get_attr_default("isvector",0):
#				all=[i.get_data_as_vector() for i in im]

				all=[]
				for j in range(im[0].get_xsize()):
					r=[]
					for i in range(len(im)):
						r.append(im[i][j,0])
					all.append(r)
				self.set_data(all,vecset,quiet=quiet)
			else:
				for idx,image in enumerate(im):
					l = [i for i in range(image.get_size())]
					k = image.get_data_as_vector()
					self.set_data([l,k],filename+":"+str(idx),quiet=quiet)

		elif file_type == 'fp':
			fin=open(filename)
			fph=struct.unpack("120sII",fin.read(128))
			ny=fph[1]
			nx=fph[2]
			data=[]
			for i in range(nx):
				data.append(struct.unpack("%df"%ny,fin.read(4*ny)))

			self.set_data(data,filename,quiet=quiet)
		else:
			try:
				# this should probably be replaced with something more flexible
				fin=open(filename)
				fin.seek(0)
				rdata=fin.readlines()
				if '#' in rdata[0]:
					try: comments=[i.split("#",1)[1].strip() for i in rdata if i[0]!="#"]
					except: comments=None
				else: comments=None
				rdata=[i.split("#")[0] for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[safe_float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]

				self.set_data(data,remove_directories_from_name(filename,1),quiet=quiet,comments=comments)
			except:
				traceback.print_exc()
				print("couldn't read",filename)
				return False

		return True

	def get_data_from_file(filename):
		file_type = Util.get_filename_ext(filename)
		em_file_type = EMUtil.get_image_ext_type(file_type)
		data = None

		if em_file_type != IMAGE_UNKNOWN or filename[:4]=="bdb:" :

			im=EMData.read_images(filename)
			if len(im) == 1:
				im = im[0]
				l = [i for i in range(im.get_size())]
				k = im.get_data_as_vector()
				data = [l,k]
			elif im[0].get_attr_default("isvector",0):
#				all=[i.get_data_as_vector() for i in im]

				all=[]
				for j in range(im[0].get_xsize()):
					r=[]
					for i in range(len(im)):
						r.append(im[i][j,0])
					all.append(r)
				data = all
			else:
				for idx,image in enumerate(im):
					l = [i for i in range(image.get_size())]
					k = image.get_data_as_vector()
					data = [l,k]

		elif file_type == 'fp':
			fin=open(filename)
			fph=struct.unpack("120sII",fin.read(128))
			ny=fph[1]
			nx=fph[2]
			data=[]
			for i in range(nx):
				data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
		else:
			try:
				fin=open(filename)
				fin.seek(0)
				rdata=fin.readlines()
				rdata=[i for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				data=[[array([rdata[j][i]]) for j in range(ny)] for i in range(nx)]

			except:
				print("couldn't read",filename)

		return data

	# make the function static
	get_data_from_file = staticmethod(get_data_from_file)

	def is_file_readable(filename):
		'''
		Called by file browsing interfaces to determine if the file is plottable
		Tries to parse two lines, convert the values into floats, make sure
		the number of columns matches etc. If something goes wrong then return False.
		This function should be tightly coupled to set_data_from_file function, should
		any adaptations occur in future
		'''
		try:
			fin=open(filename)
			fin.seek(0)
			rdata = []
			while (len(rdata) < 2):
				line = fin.readline()
				if line == '': # this is equivalent to EOF
					break
				elif line[0] == '#': continue
				else: rdata.append(line)

			#print rdata
			if len(rdata) < 2: return False

			if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
			else : rdata=[[float(j) for j in i.split()] for i in rdata]

			if len(rdata[0]) == 0 or len(rdata[1]) == 0: return False
			if  len(rdata[0]) != len(rdata[1]) : return False

			return True

		except:
			#print "couldn't read",filename
			return False

	# make the function static
	is_file_readable = staticmethod(is_file_readable)

	def render(self):

		try:
			if self.data==None or len(self.data)==0 : return
			if self.xlimits==None or self.ylimits==None or self.climits==None or self.slimits==None : return
		except:
			return

		render = False

		if self.needupd or not self.plotimg:
			self.needupd=0
			if self.main_display_list != 0:
				glDeleteLists(self.main_display_list,1)
				self.main_display_list = 0

		if self.main_display_list == 0:
			self.main_display_list = glGenLists(1)
			glNewList(self.main_display_list,GL_COMPILE)
			render = True

		lighting = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)

		EMShape.font_renderer=self.font_renderer		# Important !  Each window has to have its own font_renderer. Only one context active at a time, so this is ok.
		GL.glPushMatrix()
		# overcome depth issues
		glTranslate(0,0,5)
		for k,s in list(self.shapes.items()):
#			print k,s
			s.draw(self.scr2plot)

		GL.glPopMatrix()

		if render:
			fig=Figure((old_div(self.width(),72.0),old_div(self.height(),72.0)),dpi=72.0)
			ax=fig.add_axes((.1,.1,.88,.88),autoscale_on=False,xlim=self.xlimits,ylim=self.ylimits,xscale=self.axisparms[2],yscale=self.axisparms[3])
			#else : ax=fig.add_axes((.18,.18,.9,.9),autoscale_on=True,xscale=self.axisparms[2],yscale=self.axisparms[3])
			if self.axisparms[0] and len(self.axisparms[0])>0 : ax.set_xlabel(self.axisparms[0],size="xx-large")
			if self.axisparms[1] and len(self.axisparms[1])>0 : ax.set_ylabel(self.axisparms[1],size="xx-large")
			ax.tick_params(axis='x', labelsize="x-large")
			ax.tick_params(axis='y', labelsize="x-large")
			canvas=FigureCanvasAgg(fig)

			for i in list(self.axes.keys()):
				if not self.visibility[i]: continue
				j=self.axes[i]
#				print j
				if j[0]==-1 : x=arange(len(self.data[i][0]))
				else : x=self.data[i][self.axes[i][0]]
				if j[1]==-1 : y=arange(len(self.data[i][0]))
				else : y=self.data[i][self.axes[i][1]]

				# We draw markers (if any) first
				if self.pparm[i][4]:
					mark=symtypes[self.pparm[i][5]]
					if j[2]==-2: col=colortypes[self.pparm[i][0]]
					elif j[2]==-1: col=arange(len(self.data[i][0]))*255.0/len(self.data[i][0])
					else:
						climits=self.climits
						col=old_div((self.data[i][self.axes[i][2]]-climits[0]),(climits[1]-climits[0]))*255.0

					if j[3]==-2: sz=self.pparm[i][6]
					elif j[3]==-1: sz=arange(len(self.data[i][0]))*30.0/len(self.data[i][0])
					else:
						slimits=self.slimits
						sz=(self.data[i][self.axes[i][3]]-slimits[0])*30.0/(slimits[1]-slimits[0])

					ax.scatter(x,y,sz,col,mark,linewidths=.5*(self.pparm[i][6]>4),alpha=self.alpha)

				# Then we draw the line
				if self.pparm[i][1]:
					parm=linetypes[self.pparm[i][2]]
					try: ax.plot(x,y,parm,linewidth=self.pparm[i][3],color=colortypes[self.pparm[i][0]])
					except:
						print("Error: Plot failed\n%d %s\n%d %s"%(len(x),x,len(y),y))


			canvas.draw()
			self.plotimg = canvas.tostring_rgb()  # save this and convert to bitmap as needed

			# this try except block is because the developers of matplotlib have been changing their API
			try: # this would work for matplotlib 0.98
				self.scrlim=(ax.get_window_extent().xmin,ax.get_window_extent().ymin,ax.get_window_extent().xmax-ax.get_window_extent().xmin,ax.get_window_extent().ymax-ax.get_window_extent().ymin)
			except:
				try: # this should work for matplotlib 0.91
					self.scrlim=(ax.get_window_extent().xmin(),ax.get_window_extent().ymin(),ax.get_window_extent().xmax()-ax.get_window_extent().xmin(),ax.get_window_extent().ymax()-ax.get_window_extent().ymin())
				except:
					print('there is a problem with your matplotlib')
					return
			self.plotlim=(ax.get_xlim()[0],ax.get_ylim()[0],ax.get_xlim()[1]-ax.get_xlim()[0],ax.get_ylim()[1]-ax.get_ylim()[0])

	#		print ax.get_window_extent().xmin(),ax.get_window_extent().ymin()
	#		print ax.get_window_extent().xmax(),ax.get_window_extent().ymax()
	#		print ax.get_position()
			if not self.glflags.npt_textures_unsupported():
				self.__texture_plot(self.plotimg)
			else:
				GL.glRasterPos(0,self.height()-1)
				GL.glPixelZoom(1.0,-1.0)
		#		print "paint ",self.width(),self.height(), self.width()*self.height(),len(a)
				GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1)
				GL.glDrawPixels(self.width(),self.height(),GL.GL_RGB,GL.GL_UNSIGNED_BYTE,self.plotimg)
		else:
			try:
				glCallList(self.main_display_list)
			except: pass

		if render :
			glEndList()
			try: glCallList(self.main_display_list)
			except: pass



		if lighting : glEnable(GL_LIGHTING)


	def __texture_plot(self,image_data):


		texture_2d_was_enabled = GL.glIsEnabled(GL.GL_TEXTURE_2D)
		if not texture_2d_was_enabled:GL.glEnable(GL.GL_TEXTURE_2D)

		if self.tex_name != 0: GL.glDeleteTextures(self.tex_name)
		self.tex_name = GL.glGenTextures(1)
		GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1)
		GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
		GL.glTexImage2D(GL.GL_TEXTURE_2D,0,GL.GL_RGB,self.width(),self.height(),0,GL.GL_RGB,GL.GL_UNSIGNED_BYTE, image_data)

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

		# POSITIONING POLICY - the texture occupies the entire screen area
		glBegin(GL_QUADS)

		glTexCoord2f(0,0)
		glVertex2f(0,self.height())

		glTexCoord2f(1,0)
		glVertex2f(self.width(),self.height())

		glTexCoord2f(1,1)
		glVertex2f(self.width(),0)

		glTexCoord2f(0,1)
		glVertex2f(0,0)

		glEnd()

		if not texture_2d_was_enabled: GL.glDisable(GL.GL_TEXTURE_2D)

	def scr2plot(self,x,y) :
		""" converts screen coordinates to plot coordinates """
		try:
			if self.axisparms[2]=="linear" : x2=old_div((x-self.scrlim[0]),self.scrlim[2])*self.plotlim[2]+self.plotlim[0]
			else : x2=10.0**(old_div((x-self.scrlim[0]),self.scrlim[2])*(log10(self.plotlim[2]+self.plotlim[0])-log10(self.plotlim[0]))+log10(self.plotlim[0]))
			if self.axisparms[3]=="linear" : y2=old_div((self.height()-y-self.scrlim[1]),self.scrlim[3])*self.plotlim[3]+self.plotlim[1]
			else : y2=10.0**(old_div((self.height()-y-self.scrlim[1]),self.scrlim[3])*(log10(self.plotlim[3]+self.plotlim[1])-log10(self.plotlim[1]))+log10(self.plotlim[1]))
			return (x2,y2)
		except: return (0,0)

	def plot2scr(self,x,y) :
		""" converts plot coordinates to screen coordinates """
		try:
			if self.axisparms[2]=="linear" : x2=old_div((x-self.plotlim[0]),self.plotlim[2])*self.scrlim[2]+self.scrlim[0]
			else : x2=old_div((-(self.scrlim[2]*log(x)) + (self.scrlim[0] + self.scrlim[2])*log(10)*log10(self.plotlim[0])-self.scrlim[0]*log(10)*log10(self.plotlim[0] +self.plotlim[2])),(log(10)*(log10(self.plotlim[0]) - log10(self.plotlim[0] + self.plotlim[2]))))
			if self.axisparms[3]=="linear" :y2=self.height()-(old_div((y-self.plotlim[1]),self.plotlim[3])*self.scrlim[3]+self.scrlim[1])
			else : y2=old_div((self.scrlim[3]*log(y) + self.height()*log(10.0)*log10(self.plotlim[1])-self.scrlim[1]*log(10.0)*log10(self.plotlim[1])-self.scrlim[3]*log(10.0)*log10(self.plotlim[1]) - self.height()*log(10.0)*log10(self.plotlim[1]+self.plotlim[3]) + self.scrlim[1]*log(10)*log10(self.plotlim[1]+self.plotlim[3])), (log(10)*(log10(self.plotlim[1]) - log10(self.plotlim[1]+self.plotlim[3]))))
			return (x2,y2)
		except:
			return (0,0)

	def resize_event(self,width,height):
		self.full_refresh()

	def full_refresh(self):
		'''
		This function is called from resizeGL and from the inspector when somebody toggles the display of a line
		'''
		self.needupd=1
		self.del_shapes(("xcross","ycross","lcross","Circle"))

	def setAxisParms(self,xlabel,ylabel,xlog="linear",ylog="linear"):
# skip this since we want a guaranteed redraw somewhere
#		if self.axisparms==(xlabel,ylabel,xlog,ylog): return
		self.axisparms=(xlabel,ylabel,xlog,ylog)
		self.needupd=1
		self.updateGL()

	def setAxes(self,key,xa,ya=-1,za=-2,sa=-2,quiet=False):
		if self.axes[key]==(xa,ya,za,sa) : return
		self.axes[key]=(xa,ya,za,sa)
		self.autoscale(True)
		self.needupd=1
		if not quiet : self.updateGL()

	def setPlotParms(self,key,color,line,linetype,linewidth,sym,symtype,symsize,quiet=False):
		if color==None : color=self.pparm[key][0]
		if line==None : line=self.pparm[key][1]
		if linetype==None : linetype=self.pparm[key][2]
		if linewidth==None : linewidth=self.pparm[key][3]
		if sym==None : sym=self.pparm[key][4]
		if symtype==None : symtype=self.pparm[key][5]
		if symsize==None : symsize=self.pparm[key][6]
		if self.pparm[key]==(color,line,linetype,linewidth,sym,symtype,symsize) : return
		self.pparm[key]=(color,line,linetype,linewidth,sym,symtype,symsize)
		self.needupd=1
		if not quiet: self.updateGL()

	def add_shape(self,k,s):
		"""Add a 'shape' object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed. The 'scr' shapes are in
		screen coordinates rather than data coordinates

		0            1  2  3  4  5     6     7     8
		"rect"       R  G  B  x0 y0    x1    y1    linew
		"line"       R  G  B  x0 y0    x1    y1    linew
		"label"      R  G  B  x0 y0    text  size	linew
		"circle"     R  G  B  x0 y0    r     linew
		"scrrect"    R  G  B  x0 y0    x1    y1    linew
		"scrline"    R  G  B  x0 y0    x1    y1    linew
		"scrlabel"   R  G  B  x0 y0    text  size	linew
		"scrcircle"  R  G  B  x0 y0    r     linew
		"""
		self.shapes[k]=s
		self.shapechange=1
		#self.updateGL()

	def add_shapes(self,d):
		self.shapes.update(d)
		self.shapechange=1
		#self.updateGL()

	def del_shapes(self,k=None):
		if k:
			try:
				for i in k:
					if i in self.shapes : del self.shapes[i]
			except:
				try: del self.shapes[k]
				except: return
		else:
			self.shapes={}

		self.shapechange=1
		#self.updateGL()

	def update_selected(self,evc,lc):
		"""Update the list of 'selected' points, and display coord info.
evc is the cursor selection point in screen coords
lc is the cursor selection point in plot coords"""

		j=0
		# we find the first displayed axis in the list
		for ak in list(self.axes.keys()):
			if not self.visibility[ak]: continue
			j=self.axes[ak]
			break

		if j[0]==-1 : x=arange(len(self.data[ak][0]))
		else : x=self.data[ak][self.axes[ak][0]]
		if j[1]==-1 : y=arange(len(self.data[ak][0]))
		else : y=self.data[ak][self.axes[ak][1]]

		try: comments=self.comments[ak]
		except: comments=None

		x2,y2=self.plot2scr(x,y)

		r=hypot(x2-evc[0],y2-evc[1])		# distance from cursor

		# if the closest point is more than 5 pixels away, nothing is selected
		if min(r)>5 :
			self.selected=[]
			self.del_shapes(("selp0","selp1","selp2","selp3","selp4","selpc"))
			return

		srt=argsort(r)					# indices which would produce sorted array

		# We select one point if it's within 5 pixels, then up to 4 more points within 3 pixels
		self.selected=[srt[0]]
		for i in range(1,min(5,len(srt))) :
			if r[srt[i]]<3 : self.selected.append(srt[i])

		y0=35
		if comments!=None:
			# If the comment is "#;file", we display the image, otherwise show the comment on the plot
			p=self.selected[0]
			try:
				cmts = comments[p].split(";")

				for i in range(old_div(len(cmts),2)):
					imn = int(cmts[2*i])
					imf = cmts[2*i+1]
					# The first image window will display the first 2 images when available
					if i==0 and len(cmts)>=4: 
						ptclim=[EMData(imf,imn),EMData(cmts[3],int(cmts[2]))] 
						ptclim[0].process_inplace("normalize.edgemean")
						ptclim[1].process_inplace("normalize.edgemean")
					else: ptclim=EMData(imf,imn)
					
					try: self.particle_viewers[i].set_data(ptclim)
					except: 
						self.particle_viewers.append(emimage2d.EMImage2DWidget(ptclim))
						if len(self.particle_viewers)!=i+1 : print_exc()
					self.particle_viewers[i].show()
#				if len(cmts)==2
#					p1=EMData(cmts[1],cmts[0])
#					p2=EMData(cmts[3],cmts[2])
#					p1.process_inplace("filter.highpass.tophat",{"cutoff_freq":0.01})
#					p1.process_inplace("filter.lowpass.tophat",{"cutoff_freq":1.0/30.0})
			except:
				self.add_shape("selpc",EMShape(("scrlabel",0,0,0,80,self.scrlim[3]-(35),comments[p],120.0,-1)))
				y0+=18

		for i,p in enumerate(self.selected):
			self.add_shape("selp%d"%i,EMShape(("scrlabel",0,0,0,self.scrlim[2]-220,self.scrlim[3]-(18*i+y0),"%d. %1.3g, %1.3g"%(p,x[p],y[p]),120.0,-1)))

		self.selected_sg.emit()

	def mousePressEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector(1)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.del_shapes()
			self.updateGL()
			self.rmousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
			self.add_shape("xcross",EMShape(("scrline",0,0,0,self.scrlim[0],self.height()-event.y(),self.scrlim[2]+self.scrlim[0],self.height()-event.y(),1)))
			self.add_shape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			try: recip="%1.2f"%(old_div(1.0,lc[0]))
			except: recip="-"
			self.add_shape("lcross",EMShape(("scrlabel",0,0,0,self.scrlim[2]-220,self.scrlim[3]-10,"%1.5g (%s), %1.5g"%(lc[0],recip,lc[1]),120.0,-1)))
			self.update_selected((event.x(),event.y()),lc)
			self.updateGL()
			#if self.mmode==0:
				#self.emit(QtCore.SIGNAL("mousedown"), event)
				#return
			#elif self.mmode==1 :
				#try:
					#del self.shapes["MEASL"]
				#except: pass
				#self.add_shape("MEAS",("line",.5,.1,.5,lc[0],lc[1],lc[0]+1,lc[1],2))

	def mouseMoveEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())

		if  self.rmousedrag:
			self.add_shape("zoom",EMShape(("scrrect",0,0,0,self.rmousedrag[0],self.height()-self.rmousedrag[1],event.x(),self.height()-event.y(),1)))
			self.updateGL()
		elif event.buttons()&Qt.LeftButton:
			self.add_shape("xcross",EMShape(("scrline",0,0,0,self.scrlim[0],self.height()-event.y(),self.scrlim[2]+self.scrlim[0],self.height()-event.y(),1)))
			self.add_shape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))

			try: recip="%1.2f"%(old_div(1.0,lc[0]))
			except: recip="-"
			self.add_shape("lcross",EMShape(("scrlabel",0,0,0,self.scrlim[2]-220,self.scrlim[3]-10,"%1.5g (%s), %1.5g"%(lc[0],recip,lc[1]),120.0,-1)))
			self.update_selected((event.x(),event.y()),lc)
			self.updateGL()
#			self.add_shape("mcross",EMShape(("scrlabel",0,0,0,self.scrlim[2]-80,self.scrlim[3]-20,"%1.3g, %1.3g"%(self.plot2scr(*lc)[0],self.plot2scr(*lc)[1]),1.5,1)))

	def mouseReleaseEvent(self, event):
		lc =self.scr2plot(event.x(),event.y())
		if self.rmousedrag:
			lc2=self.scr2plot(*self.rmousedrag)
			if fabs(event.x()-self.rmousedrag[0])+fabs(event.y()-self.rmousedrag[1])<3 : self.rescale(0,0,0,0)
			else : self.rescale(min(lc[0],lc2[0]),max(lc[0],lc2[0]),min(lc[1],lc2[1]),max(lc[1],lc2[1]))
			self.rmousedrag=None
		#elif event.button()==Qt.LeftButton:
			#if self.mmode==0:
				#self.emit(QtCore.SIGNAL("mouseup"), event)
				#return
			#elif self.mmode==1 :
				#self.add_shape("MEAS",("line",.5,.1,.5,self.shapes["MEAS"][4],self.shapes["MEAS"][5],lc[0],lc[1],2))

	def rescale(self,x0,x1,y0,y1,quiet=False):
		"adjusts the value range for the x/y axes"
		self.xlimits=(x0,x1)
		self.ylimits=(y0,y1)
		if x0>=x1 or y0>=y1 : self.autoscale()
		self.needupd=1
		self.del_shapes()  # also triggers an update
		self.updateGL()
		if self.inspector: self.inspector.update()

	def recolor(self,c0,c1,quiet=False):
		"adjusts the value range for the marker color display"
		if c0>=c1 : self.autoscale()
		else: self.climits=(c0,c1)
		self.needupd=1
		self.del_shapes()  # also triggers an update
		self.updateGL()
		if self.inspector: self.inspector.update()

	def remsize(self,s0,s1,quiet=False):
		"Adjusts the value range for the marker size display"
		if s0>=s1 : self.autoscale()
		else: self.slimits=(s0,s1)
		self.needupd=1
		self.del_shapes()  # also triggers an update
		self.updateGL()
		if self.inspector: self.inspector.update()

	def autoscale(self,force=False):
		"This autoscales, but only axes which currently have invalid settings"
		if force or self.xlimits==None or self.xlimits[1]<=self.xlimits[0] :
			xmin=1.0e38
			xmax=-1.0e38
			for k in list(self.axes.keys()):
				if not self.visibility[k]: continue
				xmin=min(xmin,min(self.data[k][self.axes[k][0]]))
				xmax=max(xmax,max(self.data[k][self.axes[k][0]]))

			if self.axisparms[2]!="linear" : self.xlimits=(old_div(xmin,1.1),xmax*1.1)
			else:
				margin=(xmax-xmin)*0.025
				self.xlimits=(xmin-margin,xmax+margin)

		if force or self.ylimits==None or self.ylimits[1]<=self.ylimits[0] :
			ymin=1.0e38
			ymax=-1.0e38
			for k in list(self.axes.keys()):
				if not self.visibility[k]: continue
				ymin=min(ymin,min(self.data[k][self.axes[k][1]]))
				ymax=max(ymax,max(self.data[k][self.axes[k][1]]))

			if self.axisparms[3]!="linear" : self.ylimits=(old_div(ymin,1.1),ymax*1.1)
			else:
				margin=(ymax-ymin)*0.025
				self.ylimits=(ymin-margin,ymax+margin)

		if force or self.climits==None or self.climits[1]<=self.climits[0] :
			cmin=1.0e38
			cmax=-1.0e38
			for k in list(self.axes.keys()):
				if not self.visibility[k]: continue
				cmin=min(cmin,min(self.data[k][self.axes[k][2]]))
				cmax=max(cmax,max(self.data[k][self.axes[k][2]]))
			self.climits=(cmin,cmax)

		if force or self.slimits==None or self.slimits[1]<=self.slimits[0] :
			smin=1.0e38
			smax=-1.0e38
			for k in list(self.axes.keys()):
				if not self.visibility[k]: continue
				smin=min(smin,min(self.data[k][self.axes[k][3]]))
				smax=max(smax,max(self.data[k][self.axes[k][3]]))
			self.slimits=(smin,smax)

		if self.inspector: self.inspector.update()

	def wheelEvent(self, event):
		pass
		#if event.delta() > 0:
			#self.set_scale(self.scale + MAG_INC)
		#elif event.delta() < 0:
			#if ( self.scale - MAG_INC > 0 ):
				#self.set_scale(self.scale - MAG_INC)
		## The self.scale variable is updated now, so just update with that
		#if self.inspector: self.inspector.set_scale(self.scale)

	def leaveEvent(self,event):
		pass

class EMPolarPlot2DWidget(EMGLWidget):
	"""
	A QT widget for plotting ploar plots:
	"""
	clusterStats = QtCore.pyqtSignal(list)

	def __init__(self,application=None,winid=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		EMGLWidget.__init__(self, winid=winid)
		self.setFormat(fmt)
		self.resize(640,480)
		self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(ploticon)))

		self.axes={}
		self.pparm={}			# color,line,linetype,linewidth,sym,symtype,symsize
		self.inspector=None
		self.needupd=1
		self.plotimg=None
		self.shapes={}
		self.limits=None
		self.climits=None
		self.slimits=None
		self.rmousedrag=None
		self.axisparms=(None,None,"linear","linear")
		self.selected=[]

		self.alpha = 0.75

		self.data={}				# List of Lists to plot
		self.visibility = {}  	   	# Same entries as in self.data, but entries are true or False to indicate visibility
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags

		self.tex_name = 0
		self.main_display_list = 0

		#polar plotting stuff
		self.datap = None
		self.setDataLabelsColor('#00ff00')
		self.scattercolor = None	# IF set to none default colors are used
		self.pointsizes = None		# Defalt is to use consta sizes. Overrides constant size
		self.yticklabels = True		# Default is to draw Y tick labels
		self.xticklabels = True		# Default is to draw X tick labels

	def set_yticklabels(self, boolvalue):
		self.yticklabels = boolvalue

	def set_xticklabels(self, boolvalue):
		self.xticklabels = boolvalue

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
		GL.glEnable(GL_DEPTH_TEST)

	def paintGL(self):

		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		self.render()

	def resizeGL(self, width, height):
		#print "resize ",self.width(), self.height()
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())

		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GL.glOrtho(0.0,self.width(),0.0,self.height(),-10,10)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()

		self.resize_event(width,height)

	def closeEvent(self,event):
		self.clear_gl_memory()
		EMGLWidget.closeEvent(self, event)


	def keyPressEvent(self,event):
		if event.key() == Qt.Key_C:
			self.show_inspector(1)
		elif event.key() == Qt.Key_F1:
			try: from PyQt4 import QtWebKit
			except: return
			try:
				try:
					test = self.browser
				except:
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2display"))
					self.browser.resize(800,800)

				if not self.browser.isVisible(): self.browser.show()
			except: pass

	def setWindowTitle(self,filename):
		EMGLWidget.setWindowTitle(self, remove_directories_from_name(filename,1))

	def clear_gl_memory(self):
		if self.tex_name != 0:
			GL.glDeleteTextures(self.tex_name)
			self.tex_name = 0
		if self.main_display_list != 0:
			glDeleteLists(self.main_display_list,1)
			self.main_display_list = 0

	def mousePressEvent(self, event):
		#Save a snapshot of the scene
		self.clusterorigin_rad = None
		self.clusterradius = None
		self.clusterorigin_theta = None
		lc=self.scr2plot(event.x(),event.y())
		self.lastcx = self.firstcx = event.x()
		self.lastcy = self.firstcy = event.y()
		x = self.firstcx - old_div(self.width(),2.0)
		y = self.firstcy - old_div(self.height(),2.0)
		if event.buttons()&Qt.MidButton:
			filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Publish or Perish! Save Plot', os.getcwd(), "(*.tiff *.jpeg, *.png)")
			if filename: # if we cancel
				self.saveSnapShot(filename)
		elif event.buttons()&Qt.LeftButton:
			self.clusterorigin_rad = self._computeRadius(x,y)
			self.clusterorigin_theta = self._computeTheta(x,y)
			self.valradius = 1.0
			self.add_shape("Circle",EMShape(("scrcircle",1,0,0,self.firstcx,self.height()-self.firstcy,self.valradius,2.0)))
			self.updateGL()
		elif event.buttons()&Qt.RightButton:
			best = self.find_image(self._computeTheta(x,y), self._computeRadius(x,y))
			if best == -1:
				print("No Point Selected")
			else:
				data = self.data["data"]
				self.valradius=4.0
				cxx, cyy = self._computeXY(data[0][best], data[1][best])
				self.add_shape("Circle",EMShape(("scrcircle",1,0,0, (old_div(self.width(), 2) + cxx), (old_div(self.height(), 2) + cyy), self.valradius, 2.0)))
				self.updateGL()
#				self.emit(QtCore.SIGNAL("clusterStats"), [meanAngle,meanRad,rmsdAngle,rmsdRad,pcount])
				if event.modifiers()&Qt.ShiftModifier:
						print("Shift Clicked!")
						#if self.particle_viewer == None:
							#first = True
							#self.particle_viewer = EMImage2DWidget(data=None, application=get_application())
							#self.particle_viewer.show()

		else:
			self.find_image(self._computeTheta(x,y), self._computeRadius(x,y))
	def _computeRadius(self, x, y):
		"""Return the radius of two x and y points """
		radius = math.sqrt(x**2 + y**2)
		scaling = self.width()*(self.plotdims.x1 - self.plotdims.x0)
		scaledrad = radius*(2.0*self.plotlim[3]/scaling)
		return  scaledrad

	def _computeTheta(self, x, y):
		""" Compute the theta angle for a given x and y"""
		return -math.atan2(y,x)

	def _computeXY(self, theta, rad):
		"""return x and y given theta and rad"""
		scaling = self.width()*(self.plotdims.x1 - self.plotdims.x0)
		rescaledrad = old_div(rad,(2.0*self.plotlim[3]/scaling))
		x = math.cos(theta) * rescaledrad
		y = math.sin(theta) * rescaledrad
		return x, y

	def find_image(self, theta, rad):
		data = self.data["data"]
		bestdist = float("infinity")
		for i in range(len(data[0])):
			dist = data[1][i]**2 + rad**2 - 2*data[1][i]*rad*math.cos(data[0][i] - theta)
			if dist < bestdist:
				bestdist = dist
				best = i
#				bestpoint = self.datap[i]
#		print bestdist
#		if  bestdist > 10:
#			return 0
#		else:
		return best
#		self.valradius=4.0
#		self.add_shape("Circle",EMShape(("scrcircle",0,1,0,data[0][best], self.height()-data[1][best], self.valradius,2.0)))
#		self.updateGL()
#		print data[0][best], data[1][best]
		#print "This point correpsonds to image: %s"%bestpoint
#		self.emit(QtCore.SIGNAL("pointIdentity(int)"), bestpoint)

	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton:
			lc=self.scr2plot(event.x(),event.y())
			disp = (self.lastcx - event.x()) + (self.lastcy - event.y())
			self.valradius += disp
			self.add_shape("Circle",EMShape(("scrcircle",1,0,0,self.firstcx,self.height()-self.firstcy,self.valradius,2.0)))
			self.updateGL()
			self.lastcx = event.x()
			self.lastcy = event.y()
			#If we are drawing a cluster circle, then compute its radius for use in find circumscribed particles
			if self.clusterorigin_rad:
				x = event.x() - old_div(self.width(),2.0)
				y = event.y() - old_div(self.height(),2.0)
				rad = self._computeRadius(x,y)
				self.clusterradius = self.clusterorigin_rad**2 + rad**2 - 2*self.clusterorigin_rad*rad*math.cos(self._computeTheta(x,y) - self.clusterorigin_theta)

	def mouseReleaseEvent(self, event):
		# Find all particles within the circle
		if self.clusterorigin_rad and self.clusterradius:
			pcount = 0
			sigmaAngSin = 0.0
			sigmaAngCos = 0.0
			sigmaRad = 0.0
			statsArray = []
			data = self.data["data"]
			# Compute mean angle and radius
			for i in range(len(data[0])):
				dist = data[1][i]**2 + self.clusterorigin_rad**2 - 2*data[1][i]*self.clusterorigin_rad*math.cos(data[0][i] - self.clusterorigin_theta)
				if dist < self.clusterradius:
					# We need to compute mean of angle
					sigmaAngSin += math.sin(data[0][i])
					sigmaAngCos += math.cos(data[0][i])
					sigmaRad += data[1][i]
					statsArray.append([math.sin(data[0][i]), math.cos(data[0][i]), data[1][i]])
					pcount += 1
			if pcount == 0:
				return
			# Compute stats
			meanAngle = math.degrees(math.atan2(old_div(sigmaAngSin,pcount),old_div(sigmaAngCos,pcount)))
			meanRad = old_div(sigmaRad,pcount)
			#print "Mean Angle: %3.2f, Mean Rad: %3.2f, Num particles: %d"%(meanAngle, meanRad, pcount)
			# Compute RMSD angle and radius
			varAngSin = 0.0
			varAngCos = 0.0
			varRad = 0.0
			for i,sa in enumerate(statsArray):
				varAngSin += (sa[0]-old_div(sigmaAngSin,pcount))**2
				varAngCos += (sa[1]-old_div(sigmaAngSin,pcount))**2
				varRad += (sa[2]-old_div(sigmaRad,pcount))**2
			rmsdAngle = math.degrees(math.atan2(math.sqrt(old_div(varAngSin,pcount)),math.sqrt(old_div(varAngCos,pcount))))
			rmsdRad = math.sqrt(old_div(varRad,pcount))
			#print "RMSD Angle: %3.2f, RMSD Rad %3.2f"%(rmsdAngle, rmsdRad)
			self.clusterStats.emit([meanAngle, meanRad, rmsdAngle, rmsdRad, pcount])

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

	def base_set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,linewidth=1,linetype=0,symtype=-1,symsize=4):
		"""Set a keyed data set. The key should generally be a string describing the data.
		'data' is a tuple/list of tuples/list representing all values for a particular
		axis. eg - the points: 1,5; 2,7; 3,9 would be represented as ((1,2,3),(5,7,9)).
		Multiple axes may be set, and which axis represents which axis in the plot can be
		selected by the user. 'data' can also be an EMData object, in which case the entire
		data array is plotted as if it were 1-D."""

		self.del_shapes()

	def set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,linewidth=1,linetype=-1,symtype=-1,symsize=4,radcut=-1,datapoints=None):
		"""
		Reimplemtation to set polar data
		see set_data in EMPlot2DWidget for details
		"""
		if len(input_data) != 2:
			raise ValueError("The must be Theta and R axes")

		self.del_shapes()

		self.needupd=1

		if replace:
			self.data = {}
			self.axes = {}
			self.visibility = {}

		if input_data is None :
			if not quiet: self.updateGL()
			return

		if isinstance(input_data,EMData):
			data = input_data.get_data_as_vector()
		else: data = input_data

		try:
			if len(data)>1 : self.axes[key]=(0,1,-2,-2)
			else : self.axes[key]=(-1,0,-2,-2)
		except: return

		if color<0 : color=len(self.data)%len(colortypes)			# Automatic color setting
		if color >len(colortypes): color = 0 # there are only a certain number of colors
		if linetype<0 and symtype<0 :
			if data[0]==sorted(data[0]) : linetype=0
			else: symtype=0
		if linetype>=0 : doline=1
		else : doline,linetype=0,0
		if symtype>=0 : dosym=1
		else : dosym,symtype=0,0
		self.pparm[key]=(color,doline,linetype,linewidth,dosym,symtype,symsize)

		if not isinstance(data[0],list) and not isinstance(data[0],tuple) and not isinstance(data[0],ndarray):
			x_axis = arange(len(data))
			rdata = [ x_axis,array(data) ]
			self.data[key]= rdata
			self.visibility.setdefault(key,True)
		else:
			if data :
				self.data[key]=[array(i) for i in data]
				self.visibility.setdefault(key,True)
			else :
				#del self.data[key] why del?
				self.data.pop(key)
				self.visibility.pop(key)

		lst = list(self.pparm[key])
		lst.append(radcut)
		self.pparm[key] = tuple(lst)
		self.datap = datapoints

		if not quiet: self.updateGL()
		#if self.inspector: self.inspector.datachange()

	def setDataLabelsColor(self, color):
		""" Set the color of the data labels """
		self.datalabelscolor = color

	def setScatterColor(self, color):
		""" Set a matplotlib color or list of colors. One list for each data set """
		self.scattercolor = color

	def setPointSizes(self, sizes):
		""" Sets a list of point sizes. One list for each data set """
		self.pointsizes = sizes

	def render(self):
		"""
		Reimplmentation to plot a plor plot
		"""

		if not self.data : return

		render = False

		if self.needupd or not self.plotimg:
			self.needupd=0
			if self.main_display_list != 0:
				glDeleteLists(self.main_display_list,1)
				self.main_display_list = 0

		if self.main_display_list == 0:
			self.main_display_list = glGenLists(1)
			glNewList(self.main_display_list,GL_COMPILE)
			render = True

		lighting = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)

		EMShape.font_renderer=self.font_renderer		# Important !  Each window has to have its own font_renderer. Only one context active at a time, so this is ok.
		GL.glPushMatrix()
		# overcome depth issues
		glTranslate(0,0,5)
		for k,s in list(self.shapes.items()):
#			print k,s
			s.draw(self.scr2plot)

		GL.glPopMatrix()

		if render:

			fig=Figure((old_div(self.width(),72.0),old_div(self.height(),72.0)),dpi=72.0)
			if self.limits :ax=fig.add_axes((.1,.1,.8,.8),autoscale_on=False,Polar=True,xlim=self.limits[0],ylim=self.limits[1],xscale=self.axisparms[2],yscale=self.axisparms[3])
			else : ax=fig.add_axes((.1,.1,.8,.8),autoscale_on=True,polar=True,xscale=self.axisparms[2],yscale=self.axisparms[3])
			if self.axisparms[0] and len(self.axisparms[0])>0 : ax.set_xlabel(self.axisparms[0],size="xx-large")
			if self.axisparms[1] and len(self.axisparms[1])>0 : ax.set_ylabel(self.axisparms[1],size="xx-large")
			if not self.yticklabels: ax.set_yticklabels([])
			if not self.xticklabels: ax.set_xticklabels([])
			canvas=FigureCanvasAgg(fig)

			for i in list(self.axes.keys()):
				if not self.visibility[i]: continue
				j=self.axes[i]
				theta=self.data[i][self.axes[i][0]]
				r=self.data[i][self.axes[i][1]]
				parm=""
				if self.pparm[i][1]:
					parm+=linetypes[self.pparm[i][2]]
				if self.pparm[i][4]:
					parm+=symtypes[self.pparm[i][5]]

				# Set color(s)
				if self.scattercolor:
					scattercolor = self.scattercolor[0]
				else:
					scattercolor = colortypes[self.pparm[i][0]]
				# Set size(s)
				if self.pointsizes:
					pointsizes = self.pointsizes
				else:
					pointsizes = self.pparm[i][3]

				ax.scatter(theta, r,s=pointsizes, color=scattercolor, lw=3, alpha=self.alpha)

			if len(self.pparm[i]) == 8 and self.pparm[i][7] >= 0:
				ax.set_rmax(self.pparm[i][7])

			if self.datap:
				for i in range(len(theta)):
					ax.annotate(" "+str(self.datap[i]),(theta[i],r[i]),color=self.datalabelscolor,weight='bold',horizontalalignment='left')

			canvas.draw()
			self.plotimg = canvas.tostring_rgb()  # save this and convert to bitmap as needed

			# this try except block is because the developers of matplotlib have been changing their API
			try: # this would work for matplotlib 0.98
				self.scrlim=(ax.get_window_extent().xmin,ax.get_window_extent().ymin,ax.get_window_extent().xmax-ax.get_window_extent().xmin,ax.get_window_extent().ymax-ax.get_window_extent().ymin)
			except:
				try: # this should work for matplotlib 0.91
					self.scrlim=(ax.get_window_extent().xmin(),ax.get_window_extent().ymin(),ax.get_window_extent().xmax()-ax.get_window_extent().xmin(),ax.get_window_extent().ymax()-ax.get_window_extent().ymin())
				except:
					print('there is a problem with your matplotlib')
					return
			self.plotlim=(ax.get_xlim()[0],ax.get_ylim()[0],ax.get_xlim()[1]-ax.get_xlim()[0],ax.get_ylim()[1]-ax.get_ylim()[0])
			self.plotdims = ax.get_position()

			if not self.glflags.npt_textures_unsupported():
				self.__texture_plot(self.plotimg)
			else:
				GL.glRasterPos(0,self.height()-1)
				GL.glPixelZoom(1.0,-1.0)
				GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1)
				GL.glDrawPixels(self.width(),self.height(),GL.GL_RGB,GL.GL_UNSIGNED_BYTE,self.plotimg)
		else:
			try:
				glCallList(self.main_display_list)
			except: pass

		if render :
			glEndList()
			glCallList(self.main_display_list)

		if lighting : glEnable(GL_LIGHTING)

	def __texture_plot(self,image_data):


		texture_2d_was_enabled = GL.glIsEnabled(GL.GL_TEXTURE_2D)
		if not texture_2d_was_enabled:GL.glEnable(GL.GL_TEXTURE_2D)

		if self.tex_name != 0: GL.glDeleteTextures(self.tex_name)
		self.tex_name = GL.glGenTextures(1)
		GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1)
		GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
		GL.glTexImage2D(GL.GL_TEXTURE_2D,0,GL.GL_RGB,self.width(),self.height(),0,GL.GL_RGB,GL.GL_UNSIGNED_BYTE, image_data)

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

		# POSITIONING POLICY - the texture occupies the entire screen area
		glBegin(GL_QUADS)

		glTexCoord2f(0,0)
		glVertex2f(0,self.height())

		glTexCoord2f(1,0)
		glVertex2f(self.width(),self.height())

		glTexCoord2f(1,1)
		glVertex2f(self.width(),0)

		glTexCoord2f(0,1)
		glVertex2f(0,0)

		glEnd()

		if not texture_2d_was_enabled: GL.glDisable(GL.GL_TEXTURE_2D)

	def scr2plot(self,x,y) :
		""" converts screen coordinates to plot coordinates """
		try:
			if self.axisparms[2]=="linear" : x2=old_div((x-self.scrlim[0]),self.scrlim[2])*self.plotlim[2]+self.plotlim[0]
			else : x2=10.0**(old_div((x-self.scrlim[0]),self.scrlim[2])*(log10(self.plotlim[2]+self.plotlim[0])-log10(self.plotlim[0]))+log10(self.plotlim[0]))
			if self.axisparms[3]=="linear" : y2=old_div((self.height()-y-self.scrlim[1]),self.scrlim[3])*self.plotlim[3]+self.plotlim[1]
			else : y2=10.0**(old_div((self.height()-y-self.scrlim[1]),self.scrlim[3])*(log10(self.plotlim[3]+self.plotlim[1])-log10(self.plotlim[1]))+log10(self.plotlim[1]))
			return (x2,y2)
		except: return (0,0)

	def plot2scr(self,x,y) :
		""" converts plot coordinates to screen coordinates """
		try:
			if self.axisparms[2]=="linear" : x2=old_div((x-self.plotlim[0]),self.plotlim[2])*self.scrlim[2]+self.scrlim[0]
			else : x2=old_div((-(self.scrlim[2]*log(x)) + (self.scrlim[0] + self.scrlim[2])*log(10)*log10(self.plotlim[0])-self.scrlim[0]*log(10)*log10(self.plotlim[0] +self.plotlim[2])),(log(10)*(log10(self.plotlim[0]) - log10(self.plotlim[0] + self.plotlim[2]))))
			if self.axisparms[3]=="linear" :y2=self.height()-(old_div((y-self.plotlim[1]),self.plotlim[3])*self.scrlim[3]+self.scrlim[1])
			else : y2=old_div((self.scrlim[3]*log(y) + self.height()*log(10.0)*log10(self.plotlim[1])-self.scrlim[1]*log(10.0)*log10(self.plotlim[1])-self.scrlim[3]*log(10.0)*log10(self.plotlim[1]) - self.height()*log(10.0)*log10(self.plotlim[1]+self.plotlim[3]) + self.scrlim[1]*log(10)*log10(self.plotlim[1]+self.plotlim[3])), (log(10)*(log10(self.plotlim[1]) - log10(self.plotlim[1]+self.plotlim[3]))))
			return (x2,y2)
		except:
			return (0,0)

	def resize_event(self,width,height):
		self.full_refresh()

	def full_refresh(self):
		'''
		This function is called from resizeGL and from the inspector when somebody toggles the display of a line
		'''
		self.needupd=1
		self.del_shapes(("xcross","ycross","lcross","Circle"))

	def setAxisParms(self,xlabel,ylabel,xlog="linear",ylog="linear"):
# skip this since we want a guaranteed redraw somewhere
#		if self.axisparms==(xlabel,ylabel,xlog,ylog): return
		self.axisparms=(xlabel,ylabel,xlog,ylog)
		self.needupd=1
		self.updateGL()

	def add_shape(self,k,s):
		"""Add a 'shape' object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed. The 'scr' shapes are in
		screen coordinates rather than data coordinates

		0            1  2  3  4  5     6     7     8
		"rect"       R  G  B  x0 y0    x1    y1    linew
		"line"       R  G  B  x0 y0    x1    y1    linew
		"label"      R  G  B  x0 y0    text  size	linew
		"circle"     R  G  B  x0 y0    r     linew
		"scrrect"    R  G  B  x0 y0    x1    y1    linew
		"scrline"    R  G  B  x0 y0    x1    y1    linew
		"scrlabel"   R  G  B  x0 y0    text  size	linew
		"scrcircle"  R  G  B  x0 y0    r     linew
		"""
		self.shapes[k]=s
		self.shapechange=1
		#self.updateGL()

	def add_shapes(self,d):
		self.shapes.update(d)
		self.shapechange=1
		#self.updateGL()

	def del_shapes(self,k=None):
		if k:
			try:
				for i in k:
					if i in self.shapes : del self.shapes[i]
			except:
				try: del self.shapes[k]
				except: return
		else:
			self.shapes={}

		self.shapechange=1
		#self.updateGL()

class EMPlot2DStatsInsp(QtGui.QWidget):

	"""This class implements the statistics pop-up from the EMPlot2DInspector"""

	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		gbl0=QtWidgets.QGridLayout(self)

		self.summary=QtGui.QPushButton(self)
		self.summary.setText("Summary Table")
		gbl0.addWidget(self.summary,2,0,1,2)

		hl1 = QtWidgets.QFrame()
		hl1.setFrameStyle(QtWidgets.QFrame.HLine)
		hl1.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
		gbl0.addWidget(hl1,3,0,1,2)

		self.wlnorm=QtWidgets.QLabel(self)
		self.wlnorm.setText("Test:")
		gbl0.addWidget(self.wlnorm,4,0)

		self.wcomb_test=QtWidgets.QComboBox(self)
		#self.wcomb_test.addItem("Welch's t-test")
		#self.wcomb_test.addItem("Student's t-test")
		#self.wcomb_test.addItem("Hotelling's T-squared test")
		#self.wcomb_test.addItem("F-test")
		#self.wcomb_test.addItem("Z-test")
		#self.wcomb_test.addItem("Chi-Squared test")
		#self.wcomb_test.addItem("Sharpio-Wilk test")
		#self.wcomb_test.addItem("Wilcoxon signed-rank test")
		#self.wcomb_test.addItem("Kolomogorov-Smirnov test")
		#self.wcomb_test.addItem("Mann-Whitney U-test")
		#self.wcomb_test.addItem("Covariance")
		#self.wcomb_test.addItem("Pearson Correlation")
		gbl0.addWidget(self.wcomb_test,4,1)

		self.wsbcols=StringBox(label="Col(s)",value="0,1")
		gbl0.addWidget(self.wsbcols,6,0)

		self.wnround=ValBox(rng=(0,10),label="Round:",value=2)
		self.wnround.intonly=1
		gbl0.addWidget(self.wnround,6,1)

		self.run=QtGui.QPushButton(self)
		self.run.setText("Compute")
		gbl0.addWidget(self.run,8,0,1,2)

		hl2 = QtWidgets.QFrame()
		hl2.setFrameStyle(QtWidgets.QFrame.HLine)
		hl2.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
		gbl0.addWidget(hl2,9,0,1,2)

		self.table = QtGui.QTableWidget() #QtGui.QTextEdit()
		self.table.setRowCount(1)
		self.table.setColumnCount(1)
		self.table.setSortingEnabled(True)
		gbl0.addWidget(self.table,10,0,2,2)

		self.summary.clicked.connect(self.printSummary)
		self.run.clicked.connect(self.runTest)

		self.imgwin=None

	def printSummary(self):
		"""Computes and plots a polynomial fit (of order N) for the current x and y axes"""
		insp=self.target().get_inspector()				# inspector
		name=str(insp.setlist.currentItem().text())		# currently hilighted item
		data=self.target().data[name]					# data set we will operate on
		rnd=self.wnround.getValue()

		d = np.asarray(self.target().data[name])

		self.table.clearContents()
		column_labels = ["Axis","Mean","Median","Std","Var","Max","Min","Range","Q1","Q3","IQR","IQM","MAD","Skewness"]
		row_labels = [str(i) for i in np.arange(len(d))]
		self.table.setRowCount(len(d[0]))
		self.table.setColumnCount(len(column_labels))
		self.replaceRowLabels(row_labels)
		self.replaceColumnLabels(column_labels)

		for c in range(len(d[0])):
			col = d[:,c]
			n = len(col)
			try: mean = np.mean(col)
			except: mean = ""
			try: q3, q2, q1 = np.percentile(col, [75, 50 ,25])
			except: q3,q2,q1 = "","",""
			try: std = np.std(col)
			except: std = ""
			try: var = np.var(col)
			except: var = ""
			try: min = np.min(col)
			except: min = ""
			try: max = np.max(col)
			except: max = ""
			try: rng = max-min
			except: rng = ""
			try: iqr = q3-q1
			except: iqr = ""
			try: mad = old_div((q3+q1),2)
			except: mad = ""
			try: skew = old_div((mad - q2), mad)
			except: skew = ""
			try: iq = np.where(np.logical_and(col>=q1, col<=q3))[0]
			except: iq = ""
			try: iqm = np.mean(col[iq])
			except: iqm = ""
			stats = [c,mean,q2,std,var,max,min,rng,q1,q3,iqr,iqm,mad,skew]
			for s,stat in enumerate(stats):
				if s == 0: item = str(int(stat))
				else: item = str(round(stat,rnd))
				self.table.setItem( c, s, QtGui.QTableWidgetItem(item) )

	def runTest(self):
		stat = str(self.wcomb_test.currentText())
		cols = [str(i) for i in self.wsbcols.getValue().split(",")]
		insp = self.target().get_inspector() # inspector
		name=str(insp.setlist.currentItem().text())
		rnd = 2

		if len(cols) <= 0:
			self.textout.setText("Please specify the columns on which you wish to compute this test or statistic")
			return

		data = np.asarray(self.target().data[name]).copy()

		self.table.clearContents()

		x = np.asarray(data).T[:,cols]
		#if stat == "Covariance":
		#	self.replaceTableLabels(cols)
		#	result = np.cov(x,rowvar=False) #result = ["\t".join([str(round(j,rnd)) for j in i]) for i in cov]
		#if stat == "Pearson Correlation":
		#	self.replaceTableLabels(cols)
		#	result = np.corrcoef(x,rowvar=False) #result = ["\t".join([str(round(j,2)) for j in i]) for i in corrcoef]
		#else:
		#	print(("{} not yet implemented!".format(stat)))
		#	return

		self.table.setRowCount(result.shape[0])
		self.table.setColumnCount(result.shape[1])

		for i, r in enumerate(result):
			for j, c in enumerate(r):
				item = str(c)
				self.table.setItem( j, i, QtGui.QTableWidgetItem(item) )

	def replaceRowLabels(self,rows):
		self.table.setVerticalHeaderLabels(rows)

	def replaceColumnLabels(self,cols):
		self.table.setHorizontalHeaderLabels(cols)

	def replaceTableLabels(self,cols):
		self.table.setHorizontalHeaderLabels(cols)
		self.table.setVerticalHeaderLabels(cols)

class EMPlot2DRegrInsp(QtGui.QWidget):
	"""This class implements the regression pop-up from the EMPlot2DInspector"""

	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		gbl0=QtWidgets.QGridLayout(self)

		insp = self.target().get_inspector()

		cx = str(insp.slidex.value())
		self.wsbnax=StringBox(label="X Cols:",value=cx)
		gbl0.addWidget(self.wsbnax,2,0)

		cy = str(insp.slidey.value())
		self.wsbnay=StringBox(label="Y Cols:",value=cy)
		gbl0.addWidget(self.wsbnay,2,1)

		self.wnord=ValBox(rng=(1,25),label="Degree:",value=2)
		self.wnord.intonly=1
		gbl0.addWidget(self.wnord,4,0)

		self.wnpts=ValBox(rng=(1,10000),label="N Interp:",value=100)
		self.wnpts.intonly=1
		gbl0.addWidget(self.wnpts,4,1)

		self.wlnorm=QtWidgets.QLabel(self)
		self.wlnorm.setText("Normalization:")
		gbl0.addWidget(self.wlnorm,6,0)

		self.wcomb_norm=QtWidgets.QComboBox(self)
		self.wcomb_norm.addItem("None")
		self.wcomb_norm.addItem("Standardize")
		self.wcomb_norm.addItem("Maxmin")
		gbl0.addWidget(self.wcomb_norm,6,1)

		self.regrb=QtGui.QPushButton(self)
		self.regrb.setText("Regress")
		gbl0.addWidget(self.regrb,8,0,1,2)

		self.regrb.clicked.connect(self.doRegression)

		self.imgwin=None

	def doRegression(self):
		"""Computes and plots a polynomial fit (of order N) for the current x and y axes"""

		degree=self.wnord.getValue()	# selected order
		xaxes = self.wsbnax.getValue()
		yaxes = self.wsbnay.getValue()
		norm = str(self.wcomb_norm.currentText())
		npts = self.wnpts.getValue()
		insp = self.target().get_inspector() # inspector
		name=str(insp.setlist.currentItem().text())	# currently hilighted items

		try:
			xaxes=[int(i) for i in xaxes.split(",")]
			if max(xaxes)>=ncol : raise Exception
			yaxes=[int(i) for i in yaxes.split(",")]
			if max(yaxes)>=ncol : raise Exception
		except:
			pass #QtGui.QMessageBox.warning(self, "Axes must be a comma separated list of column numbers")
			#return

		xs = ",".join([str(i) for i in xaxes])
		ys = ",".join([str(i) for i in yaxes])
		result_name = "Regression (Degree {}; X: {}; Y: {})".format(degree,xs,ys)

		data = np.asarray(self.target().data[name])

		x = data[xaxes].T
		y = data[yaxes].T

		x = self.normalize(x,norm) # normalize features

		if norm != "none":
			normed = np.c_[x,y].T
			color = insp.color.currentIndex() + 1
			self.target().set_data(normed,"{}_{}".format(norm,xs),replace=False,quiet=False,color=color,linewidth=0,linetype=-1,symtype=0,symsize=10,comments=None)

		A = self.vandermonde(x,degree) # polynomial features

		# perform actual regression
		coefs, _, _, _ = np.linalg.lstsq(A, y) #coefs, residuals, rank, svals

		print(("Polynomial Regression (Degree {})".format(degree)))
		print(("X: {}\tY: {}".format(xs,ys)))
		print("Coefficients:")
		for i,c in enumerate(coefs):
			print(("{}:\t{}".format(i,c)))

		# construct interpolated polynomial
		xmin = np.min(x)
		xmax = np.max(x)
		xx = np.linspace(xmin,xmax,npts)
		yy = np.polyval(coefs.flatten(),xx)

		result = np.c_[xx,yy].T

		color = insp.color.currentIndex() + 1
		self.target().set_data(result,result_name,replace=False,quiet=False,color=color,linewidth=2,linetype=0,symtype=-1,symsize=0,comments=None)
		self.target().render()

	def vandermonde(self, x, order): # Source: http://stackoverflow.com/questions/11723779/2d-numpy-power-for-polynomial-expansion
		x = np.asarray(x).T[np.newaxis]
		n = x.shape[1]
		power_matrix = np.tile(np.arange(order + 1), (n, 1)).T[..., np.newaxis]
		X = np.power(x, power_matrix)
		I = np.indices((order + 1, ) * n).reshape((n, (order + 1) ** n)).T
		F = np.product(np.diagonal(X[I], 0, 1, 2), axis=2)
		return np.fliplr(F.T) # this matrix should have the form [x**n x**n-1 ... x**2 x 1]

	def normalize(self, x,norm="none"):
		if norm == "Standardize":
			mu = np.mean(x,axis=0)
			sigma = np.std(x,axis=0)
			return old_div((x-mu),sigma)
		elif norm == "Maxmin":
			xmin = np.min(x,axis=0)
			xmax = np.max(x,axis=0)
			return old_div((x-xmin),(xmax-xmin))
		elif norm == "None":
			return x

class EMPlot2DClassInsp(QtGui.QWidget):
	"""This class implements the classification pop-up from the EMPlot2DInspector"""
	def __init__(self,target) :

		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		gbl0=QtWidgets.QGridLayout(self)

		#self.wimgfile=StringBox(label="Images:")
		#gbl0.addWidget(self.wimgfile,0,0)

		#self.wimgfilebut=QtGui.QPushButton(self)
		#self.wimgfilebut.setText("Browse")
		#gbl0.addWidget(self.wimgfilebut,0,1)

		self.kmeansb=QtGui.QPushButton(self)
		self.kmeansb.setText("K-means")
		gbl0.addWidget(self.kmeansb,2,0,1,2)

		self.wnseg=ValBox(rng=(2,32),label="Nseg:",value=2)
		self.wnseg.intonly=1
		gbl0.addWidget(self.wnseg,4,0)

		self.wnax=StringBox(label="Axes:",value="all")
		gbl0.addWidget(self.wnax,4,1)

		self.wcbaxnorm=CheckBox(label="Eq Wt Axes:",value=0)
		gbl0.addWidget(self.wcbaxnorm,6,0)

		hl1 = QtWidgets.QFrame()
		hl1.setFrameStyle(QtWidgets.QFrame.HLine)
		hl1.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
		gbl0.addWidget(hl1,7,0,1,2)

		self.threshb=QtGui.QPushButton(self)
		self.threshb.setText("Threshold")
		gbl0.addWidget(self.threshb,8,0,1,2)

		self.wnval=StringBox(label="Value(s):",value="0.50")
		self.wnval.intonly=0
		gbl0.addWidget(self.wnval,10,0)

		self.wnax_thresh=StringBox(label="Axes:",value="0")
		gbl0.addWidget(self.wnax_thresh,10,1)

		self.wcomb_threshtype=QtWidgets.QComboBox(self)
		self.wcomb_threshtype.addItem("value")
		self.wcomb_threshtype.addItem("sigma")
		self.wcomb_threshtype.addItem("median")
		self.wcomb_threshtype.addItem("mean")
		self.wcomb_threshtype.addItem("percentile")
		gbl0.addWidget(self.wcomb_threshtype,12,0)

		hl2 = QtWidgets.QFrame()
		hl2.setFrameStyle(QtWidgets.QFrame.HLine)
		hl2.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
		gbl0.addWidget(hl2,13,0,1,2)

		self.wspfix=StringBox(label="Prefix:",value="split")
		gbl0.addWidget(self.wspfix,14,0)

		self.wbmakeset=QtGui.QPushButton()
		self.wbmakeset.setText("New Sets")
		gbl0.addWidget(self.wbmakeset,14,1)

		self.kmeansb.clicked.connect(self.doKMeans)
		self.threshb.clicked.connect(self.doThresh)
		self.wbmakeset.clicked.connect(self.doMakeSet)

		#QtCore.QObject.connect(self.wimgfilebut,QtCore.SIGNAL("clicked()"),self.selectImgFile)
		#QtCore.QObject.connect(self.target(),QtCore.SIGNAL("selected"),self.imgSelect)

		self.imgwin=None

	#def disableValue(self,event):
	#	EMGLWidget.closeEvent(self, event)

	#def enableValue(self,event):
	#	EMGLWidget.closeEvent(self, event)

	def doMakeSet(self):
		"""Saves selected plots as new .lst files in sets/ if 'comment' field contains image specifiers"""
		insp=self.target().get_inspector()				# inspector
		names=[str(i.text()) for i in insp.setlist.selectedItems()]		# currently hilighted items
		lsx={}

		nums=set()
		outs=[]
		for name in names:
			try: num=int(name.rsplit("_",1)[1])
			except:
				QtGui.QMessageBox.warning(self,"Error","Please hilight sets with names ending in _# !")
				return
			if num in nums:
				QtGui.QMessageBox.warning(self, "Error","Please select only one group of sets at a time !")
				return
			nums.add(num)

			outs.append("sets/{}_{}.lst".format(self.wspfix.getValue(),num))
			out=LSXFile(outs[-1])

			try: comments=self.target().comments[name]
			except:
				QtGui.QMessageBox.warning(self,"Error", "No filenames stored in {}".format(name))
				return

			for r in range(len(comments)):
				try: imn,imf=comments[r].split(";")[:2]
				except:
					QtGui.QMessageBox.warning(self,"Error", "Invalid filename {} in {}, line {}".format(comments[r],name,r))
					return

				imn=int(imn)
				if imf not in lsx : lsx[imf]=LSXFile(imf,True)	# open the LSX file for reading
				val=lsx[imf][imn]
				out[r]=val

		QtGui.QMessageBox.information(None,"Finished","New sets created: "+", ".join(outs))

	def doKMeans(self):
		"""Performs K-means classification, and produces nseg new data sets"""
		nseg=self.wnseg.getValue()		# number of segments
		axes=self.wnax.getValue()		# selected axes
		axnorm=self.wcbaxnorm.getValue()
		insp=self.target().get_inspector()				# inspector
		name=str(insp.setlist.currentItem().text())		# currently hilighted item
		data=self.target().data[name]					# data set we will operate on
		try: comments=self.target().comments[name]		# comments for copying to results
		except: comments=None
		ncol=len(data)
		nrow=len(data[0])

		if axes=="all":
			axes=list(range(ncol))
		else:
			try:
				axes=[int(i) for i in axes.split(",")]
				if max(axes)>=ncol : raise Exception
			except:
				QtGui.QMessageBox.warning(self, "Axes must be 'all' or a comma separated list of column numbers")
				return

		# Sometimes one axis dominates the classification improperly, this makes each axis equally weighted
		if axnorm:
			print("Normalize Axes")
			datafix=[old_div(i.copy(),std(i)) for i in data]
		else: datafix=data

		# build our array data into images for analysis ... this may not be the most efficient approach
		imdata=[]
		for r in range(nrow):
			imdata.append(EMData(len(axes),1,1))
			for c,ax in enumerate(axes):
				imdata[r][c]=datafix[ax][r]

		an=Analyzers.get("kmeans")
		an.set_params({"ncls":nseg,"minchange":nrow//100,"verbose":1,"slowseed":0,"mininclass":5})
		an.insert_images_list(imdata)
		centers=an.analyze()

		#print "centers ",centers

		# extract classified results as new sets
		results=[[[] for i in range(ncol)] for j in range(nseg)]	# nseg lists of data sets. Each data set is a list of ncol lists
		resultc=[[] for j in range(nseg)]							# nseg lists of comments
		for r in range(nrow):
			s=imdata[r]["class_id"]
			for c in range(ncol):
				results[s][c].append(data[c][r])
			if comments!=None: resultc[s].append(comments[r])

		for s in range(nseg) :
			if comments!=None: self.target().set_data(results[s],key="{}_{}".format(name,s),comments=resultc[s])
			else: self.target().set_data(results[s],key="{}_{}".format(name,s))

	def doThresh(self):
		"""Performs threshold-based classification, and produces nseg new data sets"""
		axes=self.wnax_thresh.getValue()						# selected axes
		vals=self.wnval.getValue()								# values
		thresh_type=str(self.wcomb_threshtype.currentText())	# type of threshold to perform
		insp=self.target().get_inspector()						# inspector
		name=str(insp.setlist.currentItem().text())				# currently hilighted item
		data=self.target().data[name]							# data set we will operate on
		try: comments=self.target().comments[name]				# comments for copying to results
		except: comments=None

		nseg = 2
		ncol=len(data)
		nrow=len(data[0])

		if axes == "all":
			axes=list(range(ncol))
		else:
			try:
				axes=[int(i) for i in axes.split(",")]
				if max(axes)>=ncol : raise Exception
			except:
				QtGui.QMessageBox.warning(self, "Axes must be 'all' or a comma separated list of column numbers")
				return

		try:
			vals=[float(i) for i in vals.split(",")]
			if len(vals) != len(axes): raise Exception
		except:
			QtGui.QMessageBox.warning(self, "You must specify one (comma separated) value for each axis.")
			return

		if thresh_type == "value":
			axvals = {a:v for a,v in zip(axes,vals)}

		elif thresh_type == "sigma":
			tmp = np.asarray(data)[:,axes]
			mu = np.mean(tmp,axis=0)
			sig = np.std(tmp,axis=0)
			axvals = {a:(m-v*s,m+v*s) for a,v,m,s in zip(axes,vals,mu,sig)}

		elif thresh_type == "median":
			tmp = np.asarray(data)[:,axes]
			vals = np.median(tmp,axis=0)
			axvals = {a:v for a,v in zip(axes,vals)}

		elif thresh_type == "mean":
			tmp = np.asarray(data)[:,axes]
			vals = np.mean(tmp,axis=0)
			axvals = {a:v for a,v in zip(axes,vals)}

		elif thresh_type == "percentile":
			for i,(a,v) in enumerate(zip(axes,vals)):
				tmp = np.asarray(data)[:,a]
				vals[i] = np.percentile(tmp,v)
			axvals = {a:v for a,v in zip(axes,vals)}

		# build our array data into images for analysis
		imdata=[]
		for r in range(nrow):
			imdata.append(EMData(len(axes),1,1))
			for c,ax in enumerate(axes):
				imdata[r][c]=data[ax][r]
				# Make sure this class doesn't already have an assigned class id.
				try: cid = imdata[r]["class_id"]
				except: imdata[r]["class_id"] = 1
				# We only want points satisfying ALL threshold values.
				if thresh_type == "sigma":
					if data[ax][r]>axvals[ax][0] and data[ax][r]<axvals[ax][1]:
						imdata[r]["class_id"] *= 1
					else: imdata[r]["class_id"] *= 0
				else:
					if data[ax][r] < axvals[ax]:
						imdata[r]["class_id"] *= 0
					else:
						imdata[r]["class_id"] *= 1

		# extract classified results as new sets
		results=[[[] for i in range(ncol)] for j in range(nseg)]	# nseg lists of data sets. Each data set is a list of ncol lists
		resultc=[[] for j in range(nseg)]							# nseg lists of comments
		for r in range(nrow):
			s=imdata[r]["class_id"]
			for c in range(ncol):
				results[s][c].append(data[c][r])
			if comments!=None: resultc[s].append(comments[r])

		for s in range(nseg):
			if comments!=None: self.target().set_data(results[s],key="{}_{}".format(name,s),comments=resultc[s])
			else: self.target().set_data(results[s],key="{}_{}".format(name,s))

	def imgSelect(self,sel=None):
		if self.imgwin==None :
			from .emimagemx import EMImageMXWidget
			self.imgwin=EMImageMXWidget()

		try:
			if sel==None: sel=self.target().selected
		except:
			print("imgSelect with no selection")
			return

		try:
			data=[EMData(self.wimgfile.getValue(),int(i)) for i in sel]
		except:
			traceback.print_exc()
			data=None

		self.imgwin.set_data(data)
		self.imgwin.show()

	def selectImgFile(self):
		from .embrowser import EMBrowserWidget
		self.browse = EMBrowserWidget(withmodal=True,multiselect=False)
		self.browse.show()
		self.browse.ok.connect(self.setImgFile)
		self.browse.cancel.connect(self.canImgFile)

	def canImgFile(self,file=None):
		return

	def setImgFile(self,file=None):

		if file==None :
			try:
				self.wimgfile.setValue(self.browse.getResult()[0])
			except:
				traceback.print_exc()
				pass
			return

	def closeEvent(self, event):
		try: self.imgwin.close()
		except: pass

class DragListWidget(QtWidgets.QListWidget):
	"This is a minor modification of the QListWidget to support drag-drop of data sets"
	def setDataSource(self,trg):
		"""We keep a weak reference to our data source so we can pull the data only when dragging actually starts"""
		self.datasource=weakref.ref(trg)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_Backspace:
			name=str(self.currentItem().text())		# currently hilighted item
			self.datasource().target().set_data(None,key=name)
		else: QtWidgets.QListWidget.keyPressEvent(self,event)

	def dragEnterEvent(self,e):
		if e.mimeData().hasText() : e.acceptProposedAction()

	def dragMoveEvent(self,e):
		if e.mimeData().hasText() : e.acceptProposedAction()

	def dropEvent(self,e):
		if e.mimeData().hasText() :
			sdata=str(e.mimeData().text()).split("\n")

			rex=re.compile("\s*([0-9Ee\-\+\.]+)(?:[\s,;:]*)")		# regular expression for parsing text with any of these separators: <space>,;:

			# parse the data
			data=None
			for s in sdata:
				if len(s.strip())==0 or s[0]=="#" : continue

				if data==None:					# first good line
					n=len(rex.findall(s))		# count numbers on the line
					data=[ [] for i in range(n)]		# initialize empty data arrays

				# parses out each number from each line and puts it in our list of lists
				for i,f in enumerate(rex.findall(s)):
					try: data[i].append(float(f))
					except: print("Error (%d): %s"%(i,f))

			# Find an unused name for the data set
			trgplot=self.datasource().target()
			name="Dropped"
			nn=1
			while name in trgplot.data :
				name="Dropped_%d"%nn
				nn+=1

			trgplot.set_data(data,name,quiet=True)
			if n==1: trgplot.setAxes(name,0)
			elif n==2: trgplot.setAxes(name,0,1)
			elif n==3: trgplot.setAxes(name,0,1,2)
			elif n==4: trgplot.setAxes(name,0,1,2,3)

			e.acceptProposedAction()

	def supportedDropActions(self):
		return Qt.DropActions(Qt.CopyAction)

	def setMovement(self,x):
		"""The ListView and ListWidget unfortunately make use of drag-drop for internal rearrangement, but we need to use it for widget->widget copy. This prevents the parent from disabling drag/drop."""
		QtWidgets.QListWidget.setMovement(self,x)
		self.setlist.setDragEnabled(True)
		self.setlist.setAcceptDrops(True)

	def setViewMode(self,x):
		"""The ListView and ListWidget unfortunately make use of drag-drop for internal rearrangement, but we need to use it for widget->widget copy. This prevents the parent from disabling drag/drop."""
		QtWidgets.QListWidget.setViewMode(self,x)
		self.setlist.setDragEnabled(True)
		self.setlist.setAcceptDrops(True)

	def startDrag(self,actions):

		data,axes,pparm=self.datasource().getCurrentData()		# get the "current" data object
		if data==None : return						# don't start a drag if nothing is selected

		# we only copy the specific axes that are used in the current plot settings !
		if axes[1]<0: axes=[axes[0]]
		elif axes[2]<0: axes=axes[:2]
		elif axes[3]<0: axes=axes[:3]

		# create the string representation of the data set
		sdata=StringIO()		# easier to write as if to a file
		for y in range(len(data[0])):
			sdata.write("%1.8g"%data[axes[0]][y])
			for x in range(1,len(axes)):
				sdata.write("\t%1.8g"%data[axes[x]][y])
			sdata.write("\n")

		# start the drag operation
		drag = QtGui.QDrag(self)
		mimeData = QtCore.QMimeData()

		mimeData.setText(sdata.getvalue())
		drag.setMimeData(mimeData)
#		drag.setPixmap(iconPixmap);

		dropact = drag.exec_(Qt.CopyAction)
#		print "Dropped ",dropact


class EMPlot2DInspector(QtGui.QWidget):

	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"plot.png"))
		self.target=weakref.ref(target)
		vbl0=QtGui.QVBoxLayout(self)

		hbl = QtWidgets.QHBoxLayout()
		hbl.setMargin(2)
		hbl.setSpacing(6)
		hbl.setObjectName("hbl")

		gbx = QtWidgets.QGroupBox("Data sets")

		vbl3 = QtGui.QVBoxLayout()
		vbl3.setMargin(4)
		vbl3.setSpacing(6)
		vbl3.setObjectName("vbl3")
		gbx.setLayout(vbl3)
		hbl.addWidget(gbx)

		# plot list
		self.setlist=DragListWidget(self)
		self.setlist.setDataSource(self)
		self.setlist.setSelectionMode(3)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.setlist.setDragEnabled(True)
		self.setlist.setAcceptDrops(True)
		vbl3.addWidget(self.setlist)

		# none and all buttons for turning plot display on and off
		hbl6 = QtWidgets.QHBoxLayout()
		hbl.setObjectName("hbl6")
		vbl3.addLayout(hbl6)

		self.nonebut=QtGui.QPushButton(self)
		self.nonebut.setText("None")
		hbl6.addWidget(self.nonebut)

		self.allbut=QtGui.QPushButton(self)
		self.allbut.setText("All")
		hbl6.addWidget(self.allbut)

		# Slider for moving within the range
		self.showslide=ValSlider(self,(0,5),"Sel1:",0,30)
		self.showslide.setIntonly(1)
		vbl3.addWidget(self.showslide)

		# number and step for the slider
		hbl7 = QtWidgets.QHBoxLayout()
		hbl.setObjectName("hbl7")
		vbl3.addLayout(hbl7)

		self.nbox=ValBox(label="ns:",value=1)
		hbl7.addWidget(self.nbox)

		self.stepbox=ValBox(label="stp:",value=1)
		hbl7.addWidget(self.stepbox)

		vbl = QtGui.QVBoxLayout()
		vbl.setMargin(0)
		vbl.setSpacing(6)
		vbl.setObjectName("vbl")
		hbl.addLayout(vbl)

		hbl0=QtWidgets.QHBoxLayout()
		hbl0.setMargin(0)
		hbl0.setSpacing(6)
		vbl.addLayout(hbl0)

		self.saveb=QtGui.QPushButton(self)
		self.saveb.setText("Save")
		hbl0.addWidget(self.saveb)

		self.concatb=QtGui.QPushButton(self)
		self.concatb.setText("Concat")
		hbl0.addWidget(self.concatb)

		self.pdfb=QtGui.QPushButton(self)
		self.pdfb.setText("PDF")
#		self.pdfb.setEnabled(False)
		hbl0.addWidget(self.pdfb)

		hbl01=QtWidgets.QHBoxLayout()
		hbl01.setMargin(0)
		hbl01.setSpacing(6)
		vbl.addLayout(hbl01)

		self.stats=QtGui.QPushButton(self)
		self.stats.setText("Statistics")
		hbl01.addWidget(self.stats)

		self.regress=QtGui.QPushButton(self)
		self.regress.setText("Regression")
		hbl01.addWidget(self.regress)

		self.classb=QtGui.QPushButton(self)
		self.classb.setText("Classification")
		hbl01.addWidget(self.classb)

		hbl1 = QtWidgets.QHBoxLayout()
		hbl1.setMargin(0)
		hbl1.setSpacing(6)

		self.color=QtWidgets.QComboBox(self)
		self.color.addItem("black")
		self.color.addItem("blue")
		self.color.addItem("red")
		self.color.addItem("green")
		self.color.addItem("yellow")
		self.color.addItem("cyan")
		self.color.addItem("magenta")
		self.color.addItem("grey")
		hbl1.addWidget(self.color)

		vbl.addLayout(hbl1)

		hbl2 = QtWidgets.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(6)
		vbl.addLayout(hbl2)

		# This is for line parms
		vbl2a = QtGui.QVBoxLayout()
		vbl2a.setMargin(0)
		vbl2a.setSpacing(6)
		hbl2.addLayout(vbl2a)

		self.lintog=QtGui.QPushButton(self)
		self.lintog.setText("Line")
		self.lintog.setCheckable(1)
		vbl2a.addWidget(self.lintog)

		self.linsel=QtWidgets.QComboBox(self)
		self.linsel.addItem("------")
		self.linsel.addItem("- - - -")
		self.linsel.addItem(".......")
		self.linsel.addItem("-.-.-.-")
		vbl2a.addWidget(self.linsel)

		self.linwid=QtGui.QSpinBox(self)
		self.linwid.setRange(1,10)
		vbl2a.addWidget(self.linwid)

		# This is for point parms
		vbl2b = QtGui.QVBoxLayout()
		vbl2b.setMargin(0)
		vbl2b.setSpacing(6)
		hbl2.addLayout(vbl2b)

		self.symtog=QtGui.QPushButton(self)
		self.symtog.setText("Symbol")
		self.symtog.setCheckable(1)
		vbl2b.addWidget(self.symtog)

		self.symsel=QtWidgets.QComboBox(self)
		self.symsel.addItem("circle")
		self.symsel.addItem("square")
		self.symsel.addItem("plus")
		self.symsel.addItem("triup")
		self.symsel.addItem("tridown")
		vbl2b.addWidget(self.symsel)

		self.symsize=QtGui.QSpinBox(self)
		self.symsize.setRange(0,25)
		vbl2b.addWidget(self.symsize)

		# This is for "heatmap"/2D hexbin parms
		#vbl2c = QtGui.QVBoxLayout()
		#vbl2c.setMargin(0)
		#vbl2c.setSpacing(6)
		#hbl2.addLayout(vbl2c)

		#self.hmtog=QtGui.QPushButton(self)
		#self.hmtog.setText("Heatmap")
		#self.hmtog.setCheckable(1)
		#vbl2c.addWidget(self.hmtog)

		#self.hmsel=QtWidgets.QComboBox(self)
		#self.hmsel.addItem("hex")
		#self.hmsel.addItem("square")
		#vbl2c.addWidget(self.hmsel)

		#self.hmbins=QtGui.QSpinBox(self)
		#self.hmbins.setRange(1,50)
		#vbl2c.addWidget(self.hmbins)

		# per plot column selectors
		gl=QtWidgets.QGridLayout()
		gl.addWidget(QtWidgets.QLabel("X Col:",self),0,0,Qt.AlignRight)
		self.slidex=QtGui.QSpinBox(self)
		self.slidex.setRange(-1,1)
		gl.addWidget(self.slidex,0,1,Qt.AlignLeft)

		gl.addWidget(QtWidgets.QLabel("Y Col:",self),1,0,Qt.AlignRight)
		self.slidey=QtGui.QSpinBox(self)
		self.slidey.setRange(-1,1)
		gl.addWidget(self.slidey,1,1,Qt.AlignLeft)

		gl.addWidget(QtWidgets.QLabel("C Col:",self),0,2,Qt.AlignRight)
		self.slidec=QtGui.QSpinBox(self)
		self.slidec.setRange(-2,1)
		gl.addWidget(self.slidec,0,3,Qt.AlignLeft)

		gl.addWidget(QtWidgets.QLabel("S Col:",self),1,2,Qt.AlignRight)
		self.slides=QtGui.QSpinBox(self)
		self.slides.setRange(-2,1)
		gl.addWidget(self.slides,1,3,Qt.AlignLeft)
		vbl.addLayout(gl)

		hbl2 = QtWidgets.QHBoxLayout()

		self.xlogtog=QtGui.QPushButton(self)
		self.xlogtog.setText("X Log")
		self.xlogtog.setCheckable(1)
		hbl2.addWidget(self.xlogtog)

		self.ylogtog=QtGui.QPushButton(self)
		self.ylogtog.setText("Y Log")
		self.ylogtog.setCheckable(1)
		hbl2.addWidget(self.ylogtog)

		#self.zlogtog=QtGui.QPushButton(self)
		#self.zlogtog.setText("Z Log")
		#self.zlogtog.setCheckable(1)
		#hbl2.addWidget(self.zlogtog)

		vbl.addLayout(hbl2)

		self.wrescale=QtGui.QPushButton(self)
		self.wrescale.setText("Rescale")
		vbl.addWidget(self.wrescale)

		vbl0.addLayout(hbl)

		hbl2a=QtWidgets.QHBoxLayout()

		self.wl1=QtWidgets.QLabel("Min")
		self.wl1.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl1)
		self.wl2=QtWidgets.QLabel("Max")
		self.wl2.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl2)
		self.wl3=QtWidgets.QLabel("Min")
		self.wl3.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl3)
		self.wl4=QtWidgets.QLabel("Max")
		self.wl4.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl4)
		vbl0.addLayout(hbl2a)

		hbl2=QtWidgets.QHBoxLayout()

		#hbl2.addWidget(QtWidgets.QLabel("X:",self))
		#self.wxmin=QtWidgets.QLineEdit(self)
		self.wxmin=ValBox(label="X:")
		hbl2.addWidget(self.wxmin)
		#hbl2.addWidget(QtWidgets.QLabel("-",self))
		#self.wxmax=QtWidgets.QLineEdit(self)
		self.wxmax=ValBox(label="  ")
		hbl2.addWidget(self.wxmax)

		self.wymin=ValBox(label="Y:")
		hbl2.addWidget(self.wymin)
		self.wymax=ValBox(label="  ")
		hbl2.addWidget(self.wymax)

		vbl0.addLayout(hbl2)

		hbl3=QtWidgets.QHBoxLayout()

		self.wcmin=ValBox(label="C:")
		hbl3.addWidget(self.wcmin)
		self.wcmax=ValBox(label="  ")
		hbl3.addWidget(self.wcmax)

		self.wsmin=ValBox(label="S:")
		hbl3.addWidget(self.wsmin)
		self.wsmax=ValBox(label="  ")
		hbl3.addWidget(self.wsmax)
		vbl0.addLayout(hbl3)


		hbl4 = QtWidgets.QHBoxLayout()
		hbl4.addWidget(QtWidgets.QLabel("X Label:",self))
		self.xlabel=QtWidgets.QLineEdit(self)
		hbl4.addWidget(self.xlabel)
		vbl0.addLayout(hbl4)

		hbl5 = QtWidgets.QHBoxLayout()
		hbl5.addWidget(QtWidgets.QLabel("Y Label:",self))
		self.ylabel=QtWidgets.QLineEdit(self)
		hbl5.addWidget(self.ylabel)
		vbl0.addLayout(hbl5)

		hbl6 = QtWidgets.QHBoxLayout()
		#hbl6.addWidget(QtWidgets.QLabel("Transparency:",self))
		self.alphaslider=ValSlider(self,(0,1),"Transparency:",0.5,50)
		hbl6.addWidget(self.alphaslider)
		vbl0.addLayout(hbl6)

#		self.setLayout(vbl0)

		self.quiet=0
		self.busy=0
		self.classwin=None
		self.statswin=None
		self.regresswin=None

		self.showslide.valueChanged.connect(self.selSlide)
		self.allbut.clicked.connect(self.selAll)
		self.nonebut.clicked.connect(self.selNone)

		self.slidex.valueChanged[int].connect(self.newCols)
		self.slidey.valueChanged[int].connect(self.newCols)
		self.slidec.valueChanged[int].connect(self.newCols)
		self.slides.valueChanged[int].connect(self.newCols)
		self.setlist.currentRowChanged[int].connect(self.newSet)
		self.setlist.itemChanged[QtGui.QListWidgetItem].connect(self.list_item_changed)
		self.color.currentIndexChanged[str].connect(self.updPlotColor)
		self.classb.clicked.connect(self.openClassWin)
		#QtCore.QObject.connect(self.hmsel,QtCore.SIGNAL("clicked()"),self.updPlot)
		self.symtog.clicked.connect(self.updPlot)
		#QtCore.QObject.connect(self.hmsel,QtCore.SIGNAL("clicked()"),self.updPlotHmsel)
		#QtCore.QObject.connect(self.hmbins,QtCore.SIGNAL("clicked()"),self.updPlotHmbins)
		self.symsel.currentIndexChanged[str].connect(self.updPlotSymsel)
		self.symsize.valueChanged[int].connect(self.updPlotSymsize)
		self.xlogtog.clicked.connect(self.updPlot)
		self.ylogtog.clicked.connect(self.updPlot)
		#QtCore.QObject.connect(self.zlogtog,QtCore.SIGNAL("clicked()"),self.updPlot)
		self.lintog.clicked.connect(self.updPlot)
		#QtCore.QObject.connect(self.hmtog,QtCore.SIGNAL("clicked()"),self.updPlot)
		self.linsel.currentIndexChanged[str].connect(self.updPlotLinsel)
		self.linwid.valueChanged[int].connect(self.updPlotLinwid)
		self.xlabel.textChanged[str].connect(self.updPlot)
		self.ylabel.textChanged[str].connect(self.updPlot)
		self.stats.clicked.connect(self.openStatsWin)
		self.regress.clicked.connect(self.openRegrWin)
		self.saveb.clicked.connect(self.savePlot)
		self.pdfb.clicked.connect(self.savePdf)
		self.concatb.clicked.connect(self.saveConcatPlot)
		self.wxmin.valueChanged.connect(self.newLimits)
		self.wxmax.valueChanged.connect(self.newLimits)
		self.wymin.valueChanged.connect(self.newLimits)
		self.wymax.valueChanged.connect(self.newLimits)
		self.wcmin.valueChanged.connect(self.newCLimits)
		self.wcmax.valueChanged.connect(self.newCLimits)
		self.wsmin.valueChanged.connect(self.newSLimits)
		self.wsmax.valueChanged.connect(self.newSLimits)
		self.wrescale.clicked.connect(self.autoScale)
		self.alphaslider.valueChanged.connect(self.updAlpha)

		self.newSet(0)
		self.datachange()
		self.resize(500,540)

	def selSlide(self,val):
		rngn0=int(val)
		rngn1=int(self.nbox.getValue())
		rngstp=int(self.stepbox.getValue())
		rng=list(range(rngn0,rngn0+rngstp*rngn1,rngstp))
		for i,k in enumerate(sorted(self.target().visibility.keys())) :
			self.target().visibility[k]=i in rng
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def newAlpha(self,val):
		self.target().updateGL()

	def selAll(self):
		for k in list(self.target().visibility.keys()) : self.target().visibility[k]=True
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def selNone(self):
		for k in list(self.target().visibility.keys()) : self.target().visibility[k]=False
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def openStatsWin(self):
		"""This launches a separate window for regressing points in a 2-D plot"""

		if self.statswin==None : self.statswin=EMPlot2DStatsInsp(self.target())
		self.statswin.show()

	def openRegrWin(self):
		"""This launches a separate window for regressing points in a 2-D plot"""

		if self.regresswin==None : self.regresswin=EMPlot2DRegrInsp(self.target())
		self.regresswin.show()

	def openClassWin(self):
		"""This launches a separate window for classifying points in a 2-D plot"""

		if self.classwin==None : self.classwin=EMPlot2DClassInsp(self.target())
		self.classwin.show()

	def on_bad_button(self,unused=False):
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names) == 0: return

		self.__remove_from_file("xian_good.txt",names)
		self.__at_to_file("xian_bad.txt",names)

	def __remove_from_file(self,fname,names):

		if os.path.exists(fname):
			inf = open(fname,"r")
			lines = inf.readlines()
			inf.close()
		else:
			return

		f = open(fname,"w")

		for i in range(0,len(lines)):
			lines[i] = lines[i].strip()

		for i in range(len(lines)-1,-1,-1):
			if lines[i] in names:
				lines.pop(i)


		for line in lines:
			f.write(line+"\n")

		f.close()

	def __at_to_file(self,fname,names):

		if os.path.exists(fname):
			inf = open(fname,"r")
			lines = inf.readlines()
			inf.close()
		else:
			lines = []
		f = open(fname,"w")

		for i in range(0,len(lines)):
			lines[i] = lines[i].strip()

		for name in names:
			if name not in lines:
				lines.append(name)


		for line in lines:
			f.write(line+"\n")

		f.close()


	def on_good_button(self,unused=False):
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names) == 0: return


		self.__remove_from_file("xian_bad.txt",names)
		self.__at_to_file("xian_good.txt",names)

	def getCurrentData(self):
		"""Returns (data,axes,pparm) for the single 'current' plot"""
		try:
			name=str(self.setlist.currentItem().text())
			data=self.target().data[name]
			pparm=self.target().pparm[name]
			axes=self.target().axes[name]
		except: return None,None,None

		return data,axes,pparm

	def saveConcatPlot(self):
		"""Saves the contents of the current plot to a text file. All sets are concatenated together into a single file"""
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names) == 0: return

		name2="plt_concat.txt"
		i=0
		while os.path.exists(name2):
			name2="plt_concat_%02d.txt"%(i)
			i+=1
		out=open(name2,"a")

		xcol=self.slidex.value()
		ycol=self.slidey.value()
		for name in names :
			data=self.target().data[name]

			for i in range(len(data[0])):
				out.write("%g\t%g\n"%(data[xcol][i],data[ycol][i]))

		out=None
		print("Wrote ",name2)


	def savePlot(self):
		"""Saves the contents of the current plot to a text file"""
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names) == 0: return

		for name in names :
			sname=name.replace(" ","_")
			sname=sname.replace("(","_")
			sname=sname.replace(")","")
			sname=sname.replace("/","_")
			data=self.target().data[name]
			name2="plt_%s.txt"%(sname)
			i=0
			while os.path.exists(name2):
				name2="plt_%s_%02d.txt"%(sname,i)
				i+=1

			out=open(name2,"w")
			xcol=self.slidex.value()
			ycol=self.slidey.value()
			for i in range(len(data[0])):
				out.write("%g\t%g\n"%(data[xcol][i],data[ycol][i]))

			print("Wrote ",name2)

	def savePdf(self):
		"""Saves the contents of the current plot to a pdf"""
		plt.savefig("plot.pdf")

	def updPlot(self,s=None):
		if self.quiet : return
		if self.xlogtog.isChecked() : xl="log"
		else : xl="linear"
		if self.ylogtog.isChecked() : yl="log"
		else : yl="linear"
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		self.target().setAxisParms(self.xlabel.text(),self.ylabel.text(),xl,yl)
		self.target().autoscale(True)
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,None,self.lintog.isChecked(),None,None,self.symtog.isChecked(),None,None,True)
			self.target().updateGL()

	#

	def updPlotColor(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,self.color.currentIndex(),None,None,None,None,None,None,True)
			self.target().updateGL()

	#def updPlotHmsel(self,s=None):
		#if self.quiet : return
		#names = [str(item.text()) for item in self.setlist.selectedItems()]
		#if len(names)==1:
			#self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				#self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		#else:
			#for name in names:
				#self.target().setPlotParms(name,None,None,None,None,None,self.symsel.currentIndex(),None,True)
			#self.target().updateGL()

	def updPlotSymsel(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,None,None,None,None,None,self.symsel.currentIndex(),None,True)
			self.target().updateGL()

	def updPlotSymsize(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,None,None,None,None,None,None,self.symsize.value(),True)
			self.target().updateGL()

	#def updPlotHmbins(self,s=None):
		#if self.quiet : return
		#names = [str(item.text()) for item in self.setlist.selectedItems()]
		#if len(names)==1:
			#self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				#self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		#else:
			#for name in names:
				#self.target().setPlotParms(name,None,None,None,None,None,self.symsel.currentIndex(),None,True)
			#self.target().updateGL()

	def updPlotLinsel(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,None,None,self.linsel.currentIndex(),None,None,None,None,True)
			self.target().updateGL()

	def updPlotLinwid(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],self.color.currentIndex(),self.lintog.isChecked(),
				self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),self.symsel.currentIndex(),self.symsize.value())
		else:
			for name in names:
				self.target().setPlotParms(name,None,None,None,self.linwid.value(),None,None,None,True)
			self.target().updateGL()

	def updAlpha(self,val):
		if self.quiet : return
		self.target().alpha = val
		self.target().needupd=1
		self.target().updateGL()

	def newSet(self,row):
		self.quiet=1
		try:
			i=str(self.setlist.item(row).text())
		except:
#			print "plot error"
			return

		self.slidex.setRange(-1,len(self.target().data[i])-1)
		self.slidey.setRange(-1,len(self.target().data[i])-1)
		self.slidec.setRange(-2,len(self.target().data[i])-1)
		self.slides.setRange(-2,len(self.target().data[i])-1)
		self.slidex.setValue(self.target().axes[i][0])
		self.slidey.setValue(self.target().axes[i][1])
		self.slidec.setValue(self.target().axes[i][2])
		self.slides.setValue(self.target().axes[i][3])

		pp=self.target().pparm[i]
		self.color.setCurrentIndex(pp[0])

		self.lintog.setChecked(pp[1])
		self.linsel.setCurrentIndex(pp[2])
		self.linwid.setValue(pp[3])

		self.symtog.setChecked(pp[4])
		self.symsel.setCurrentIndex(pp[5])
		self.symsize.setValue(pp[6])
		self.quiet=0

	def newCols(self,val):
		if self.target and not self.quiet:
#			print "newcols",self.slidex.value(),self.slidey.value(),self.slidec.value(),self.slides.value()
			names = [str(item.text()) for item in self.setlist.selectedItems()]
			for name in names:
				self.target().setAxes(name,self.slidex.value(),self.slidey.value(),self.slidec.value(),self.slides.value(),True)
			self.target().updateGL()

#			self.target().setAxes(str(self.setlist.currentItem().text()),self.slidex.value(),self.slidey.value(),self.slidec.value(),self.slides.value())

	def newLimits(self,val=None):
		if self.busy: return
		try:
			xmin=self.wxmin.getValue()
			xmax=self.wxmax.getValue()
			ymin=self.wymin.getValue()
			ymax=self.wymax.getValue()
			self.target().rescale(xmin,xmax,ymin,ymax,True)
		except:
			self.target().rescale(0,0,0,0)

	def newCLimits(self,val=None):
		if self.busy: return
		try:
			cmin=self.wcmin.getValue()
			cmax=self.wcmax.getValue()
			self.target().recolor(cmin,cmax,True)
		except:
			self.target().recolor(0,0)

	def newSLimits(self,val=None):
		if self.busy: return
		try:
			smin=self.wsmin.getValue()
			smax=self.wsmax.getValue()
			self.target().remsize(smin,smax,True)
		except:
			self.target().remsize(0,0)

	def update(self):
		self.busy=1
		try:
			self.wxmin.setValue(self.target().xlimits[0])
			self.wxmax.setValue(self.target().xlimits[1])
			self.wymin.setValue(self.target().ylimits[0])
			self.wymax.setValue(self.target().ylimits[1])
		except:
			self.wxmin.setValue(None)
			self.wxmax.setValue(None)
			self.wymin.setValue(None)
			self.wymax.setValue(None)
		self.busy=0

		try:
			self.wcmin.setValue(self.target().climits[0])
			self.wcmax.setValue(self.target().climits[1])
		except:
			self.wcmin.setValue(None)
			self.wcmax.setValue(None)

		try:
			self.wsmin.setValue(self.target().slimits[0])
			self.wsmax.setValue(self.target().slimits[1])
		except:
			self.wsmin.setValue(None)
			self.wsmax.setValue(None)

	def autoScale(self):
		self.target().autoscale(True)
		self.target().rescale(0,0,0,0)
#		self.update()

	def datachange(self):

		self.setlist.clear()

		flags= Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)|Qt.ItemFlags(Qt.ItemIsDragEnabled)

		keys=list(self.target().data.keys())
		visible = self.target().visibility
		keys.sort()
		parms = self.target().pparm # get the colors from this

		for i,j in enumerate(keys) :
			a = QtGui.QListWidgetItem(j)
			a.setFlags(flags)
			try: a.setTextColor(qt_color_map[colortypes[parms[j][0]]])
			except:
				print("Color error")
				print(list(sorted(parms.keys())))
				print(parms[j][0])
				print(colortypes[parms[j][0]])
				print(qt_color_map[colortypes[parms[j][0]]])
			if visible[j]: a.setCheckState(Qt.Checked)
			else: a.setCheckState(Qt.Unchecked)

			self.setlist.addItem(a)

		if len(keys) > 0 : self.setlist.setCurrentRow(0)
		self.showslide.setRange(0,len(keys))

	def list_item_changed(self,item):
		checked = False
		if item.checkState() == Qt.Checked: checked = True

		name = str(item.text())
		if self.target().visibility[name] != checked:
			self.target().visibility[name] = checked
			self.target().full_refresh()
			self.target().updateGL()

	def closeEvent(self, event):
		try: self.classwin.close()
		except: pass
		try: self.statswin.close()
		except: pass
		try: self.regresswin.close()
		except: pass

class EMDataFnPlotter(QtGui.QWidget):

	def __init__(self, parent = None, data=None):
		QtGui.QWidget.__init__(self, parent)

		self.setWindowTitle("Plotter")

		self.resize(780, 580)
		self.gbl = QtWidgets.QGridLayout(self)

		self.plot = EMPlot2DWidget(parent=self)
		self.gbl.addWidget(self.plot,0,0,1,1)

		self.lplot = QtWidgets.QLabel("Plot")
		self.gbl.addWidget(self.plot,1,0)

		if data!=None :
			self.plot.set_data(data,"init")

	def closeEvent(self,event):
		self.plot.closeEvent(event)

	def set_label(self,lbl):
		self.lplot.setText(lbl)

	def set_data(self, data, key):
		self.plot.set_data(data,key)

#
# class EMGLPlot2DItem(QGLWidget):
#
# 	# default window size
# 	width, height = 640, 480
#
# 	def set_data(self, data):
# 		self.data = data
# 		self.count = data.shape[0]
#
# 	def initializeGL(self):
# 		gl.glClearColor(0,0,0,0)
# 		self.vbo = glvbo.VBO(self.data)
#
# 	def paintGL(self):
# 		# clear the buffer
# 		gl.glClear(gl.GL_COLOR_BUFFER_BIT)
# 		# set yellow color for subsequent drawing rendering calls
# 		gl.glColor(1,1,0)
# 		# bind the VBO
# 		self.vbo.bind()
# 		# tell OpenGL that the VBO contains an array of vertices
# 		gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
# 		# these vertices contain 2 single precision coordinates
# 		gl.glVertexPointer(2, gl.GL_FLOAT, 0, self.vbo)
# 		# draw "count" points from the VBO
# 		gl.glDrawArrays(gl.GL_POINTS, 0, self.count)
#
# 	def resizeGL(self, width, height):
# 		# update the window size
# 		self.width, self.height = width, height
# 		# paint within the whole window
# 		gl.glViewport(0, 0, width, height)
# 		# set orthographic projection (2D only)
# 		gl.glMatrixMode(gl.GL_PROJECTION)
# 		gl.glLoadIdentity()
# 		# the window corner OpenGL coordinates are (-+1, -+1)
# 		gl.glOrtho(-1, 1, 1, -1, -1, 1)
#
# class EMGLPlot2DWidget(QtGui.QMainWindow):
#
# 	def __init__(self):
# 		super(EMGLPlot2DWidget, self).__init__()
# 		self.data = None
#
# 	def render(self):
# 		self.plot = EMGLPlot2DItem()
# 		self.plot.set_data(self.data)
#
# 		self.setGeometry(400, 400, self.plot.width, self.plot.height)
# 		self.setCentralWidget(self.plot)
# 		self.show()
#
# 	def set_data_from_file(self, f):
# 		try: self.data = np.loadtxt(f,dtype=np.float32)
# 		except: print("Could not read {}".format(f))
# 		self.render()

# This is just for testing, of course
def main():
	
	app = EMApp()
	window = EMPlot2DWidget(app)
	if len(sys.argv)==1 :
		l=[i/30.*pi for i in range(30)]
		window.set_data([[1,2,3,4],[2,3,4,3]],"test")
		window.set_data([l,[sin(2*i) for i in l]],"test2")
	else:
		for i in range(1,len(sys.argv)):
			window.set_data_from_file(sys.argv[i])
	
	app.show()
	app.execute()


if __name__ == '__main__':
	main()

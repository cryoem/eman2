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

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtOpenGL import QGLWidget
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU
from OpenGL.GL import *
import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo
from math import *
from EMAN2 import *
import sys
from emshape import *
import weakref
from cPickle import dumps,loads
import struct, math
from numpy import *
from valslider import *
from cStringIO import StringIO
import re
#import emimage2d

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
#matplotlib.use('Agg')

import numpy as np

from emapplication import EMApp, EMGLWidget
from emglobjects import EMOpenGLFlagsAndTools

import traceback

#plt.style.use('fivethirtyeight')

alignments=["mid","right","left"]
histtypes=["bar","step","stepfilled"]
orientations=["vertical"] #,"horizontal"] # need to implement
colortypes=["k","b","r","g","y","c","m","gray"]
matplotlib.rcParams['axes.color_cycle'] = colortypes[1:-1] # temporary hack

qt_color_map = {}
qt_color_map["k"] = QtGui.QColor(0,0,0)
qt_color_map["b"] = QtGui.QColor(0,0,255)
qt_color_map["r"] = QtGui.QColor(255,0,0)
qt_color_map["g"] = QtGui.QColor(0,255,0)
qt_color_map["y"] = QtGui.QColor(255,255,0)
qt_color_map["c"] = QtGui.QColor(0,255,255)
qt_color_map["m"] = QtGui.QColor(255,0,255)
qt_color_map["gray"] = QtGui.QColor(127,127,127)

class EMHistogramWidget(EMGLWidget):
	"""A QT widget for drawing 2-D plots using matplotlib
	"""

	def __init__(self,application=None,winid=None,parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		EMGLWidget.__init__(self, parent=parent, winid=winid)
		self.setFormat(fmt)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"plot.png"))
		self.axes={}
		self.pparm={}			# nbins,color,histtype,orient,align,alpha,width,norm,cumul,logy,stacked
		self.inspector=None
		self.needupd=1
		self.plotimg=None
		self.shapes={}
		self.xlimits=None
		self.ylimits=None
		self.rmousedrag=None
		self.axisparms=(None,None,"linear","linear")
		self.selected=[]
		self.comments={}			# IF reading from a file which contains per-point comments, this dictionary contains a list of comments for each point
		self.data={}				# List of Lists to plot
		self.visibility = {}		# Same entries as in self.data, but entries are true or False to indicate visibility
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags
		self.tex_name = 0
		self.main_display_list = 0
		self.resize(640,480)
	
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
		if self.inspector : self.inspector.closeEvent(self, event)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_C:
			self.show_inspector(1)
		elif event.key() == Qt.Key_F1:
			try: from PyQt4 import QtWebKit
			except: return
			try:
				try: test = self.browser
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

	def set_data(self,input_data,key="data",replace=False,quiet=False,nbins=10,color=-1,htype="bar",orient="vertical",align="mid",alpha=0.8,width=0.8,norm=False,cumul=False,logy=False,stacked=False):
		"""Set a keyed data set. The key should generally be a string describing the data.
		'data' is a tuple/list of tuples/list representing all values for a particular
		axis. eg - the points: 1,5; 2,7; 3,9 would be represented as ((1,2,3),(5,7,9)).
		Multiple axes may be set, and which axis represents which axis in the plot can be
		selected by the user. 'data' can also be an EMData object, in which case the entire
		data array is plotted as if it were 1-D.

		linetype and symtype are integers. -1 will disable either. Default is autoselect."""
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
		if self.data.has_key(key) : oldkey=True
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
						if data[0][2]-data[0][1]==1 : self.axes[key]=(1,)#,2,-2,-2)	# if it looks like the first axis is a boring count
						else : self.axes[key]=(0,)#,1,-2,-2)
					except: self.axes[key]=(0,)#,1,-2,-2)
				else : self.axes[key]=(0,)#,1,-2,-2)
			else : self.axes[key]=(-1,)#,0,-2,-2)
		except:
			print "Data error:", data
			return
		
# 		self.inspector.slidecol.setRange(-1,len(self.data[key])-1)
# 		self.inspector.slidecol.setValue(self.axes[key][0])
		
		if oldkey:
			pp=self.pparm[key]
			if nbins < 0: nbins = pp[0]
			if color<0 : color= pp[1]
			if color > len(colortypes): color = pp[1]
			if alpha > 1.0: alpha = pp[5]
			if alpha < 0.0: alpha = pp[5]
		else:
			if nbins < 1: nbins = 1
			if color < 0: color=len(self.data)%len(colortypes)
			if color > len(colortypes): color = 0
			if alpha > 1.0: alpha = 1.0
			if alpha < 0.0: alpha = 0.8
			if width < 0.0: width = 0.8
		self.pparm[key]=(nbins,color,htype,orient,align,alpha,width,norm,cumul,logy,stacked)
		
		self.autoscale()
		if self.inspector: self.inspector.datachange()
		if not quiet: self.updateGL()

	def get_inspector(self):
		if not self.inspector :
			self.inspector=EMHistogramInspector(self)
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
				if self.data.has_key(filename) : filename="{}.{}".format(filename,len(self.data))
				self.set_data([l,k],filename,quiet=quiet)
			elif im[0].get_attr_default("isvector",0):
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
			fin=file(filename)
			fph=struct.unpack("120sII",fin.read(128))
			ny=fph[1]
			nx=fph[2]
			data=[]
			for i in range(nx):
				data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
			self.set_data(data,filename,quiet=quiet)
		else:
			try:
				fin=file(filename)
				fin.seek(0)
				rdata=fin.readlines()
				if '#' in rdata[0]:
					try: comments=[i.split("#",1)[1].strip() for i in rdata if i[0]!="#"]
					except: comments=None
				else: comments=None
				rdata=[i.split("#")[0] for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]
				self.set_data(data,remove_directories_from_name(filename,1),quiet=quiet)#,comments=comments)
			except:
				traceback.print_exc()
				print "couldn't read",filename
				return False
		return True

	@staticmethod
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
			fin=file(filename)
			fph=struct.unpack("120sII",fin.read(128))
			ny=fph[1]
			nx=fph[2]
			data=[]
			for i in range(nx):
				data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
		else:
			try:
				fin=file(filename)
				fin.seek(0)
				rdata=fin.readlines()
				rdata=[i for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				data=[[array([rdata[j][i]]) for j in range(ny)] for i in range(nx)]
			except:
				print "couldn't read",filename
		return data

	@staticmethod
	def is_file_readable(filename):
		'''
		Called by file browsing interfaces to determine if the file is plottable
		Tries to parse two lines, convert the values into floats, make sure
		the number of columns matches etc. If something goes wrong then return False.
		This function should be tightly coupled to set_data_from_file function, should
		any adaptations occur in future
		'''
		try:
			fin=file(filename)
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
			return False

	def render(self):
		try:
			if self.data==None or len(self.data)==0: return
			if self.xlimits==None or self.ylimits==None: return
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
		for k,s in self.shapes.items():
			s.draw(self.scr2plot)
		GL.glPopMatrix()
		if render:
			fig=Figure((self.width()/72.0,self.height()/72.0),dpi=72.0)
			ax=fig.add_axes((.1,.1,.88,.88),autoscale_on=False,xlim=self.xlimits,ylim=self.ylimits,xscale=self.axisparms[2],yscale=self.axisparms[3])
			#if self.axisparms[0] and len(self.axisparms[0])>0 : ax.set_xlabel(self.axisparms[0],size="xx-large")
			#if self.axisparms[1] and len(self.axisparms[1])>0 : ax.set_ylabel(self.axisparms[1],size="xx-large")
			ax.tick_params(axis='x', labelsize="x-large")
			ax.tick_params(axis='y', labelsize="x-large")
			canvas=FigureCanvasAgg(fig)
			
			if self.inspector == None: 
				self.inspector = self.get_inspector()
			
			x = []
			colors = []
			for i in self.axes.keys():
				if not self.visibility[i]: continue
				x.append(self.data[i][self.axes[i][0]])
				colors.append(colortypes[self.inspector.color.currentIndex()])
			
			stacked = self.inspector.stacked.isChecked()
			bins=int(self.inspector.slidenbs.value())
			rng = self.xlimits
			normed=self.inspector.normed.isChecked()
			histtype=histtypes[self.inspector.histtype.currentIndex()]
			orient=orientations[self.inspector.orient.currentIndex()]
			align=alignments[self.inspector.align.currentIndex()]
			alpha = self.inspector.alpha.getValue()
			rwidth = self.inspector.rwidth.getValue()
			cumulative = self.inspector.cumulative.isChecked()
			logtog = self.inspector.logtogy.isChecked()
			
			self.hist,self.hist_edges = np.histogram(x,bins,range=rng,density=normed)
			ax.hist(x,bins=self.hist_edges,color=None,range=rng,histtype=histtype,orientation=orient,align=align,alpha=alpha,rwidth=rwidth,normed=normed,cumulative=cumulative,log=logtog,stacked=stacked)
			
			if logtog:
				if cumulative: 
					ymax = max(np.cumsum(np.log10(self.hist)))
				else: 
					ymax = max(self.hist)
			elif cumulative: 
				ymax = max(np.cumsum(self.hist))
			else: 
				ymax = max(self.hist)
		
			if not stacked: ymax /= len(x)
			
			margin = ymax*0.025
			ax.set_ylim([0-margin,ymax+margin])
			self.ylimits=(0-margin,ymax+margin)
				
			# else:
# 				for i in self.axes.keys():
# 					if not self.visibility[i]: continue
# 				
# 					x=self.data[i][self.axes[i][0]]
# 				
# 					bins=int(self.pparm[i][0])
# 					color=colortypes[self.pparm[i][1]]
# 					histtype=self.pparm[i][2]
# 					orient=self.pparm[i][3]
# 					align=self.pparm[i][4]
# 					alpha=self.pparm[i][5]
# 					rwidth=self.pparm[i][6]
# 					normed=self.pparm[i][7]
# 					cumulative=self.pparm[i][8]
# 					log=self.pparm[i][9]
# 					#stacked=self.pparm[i][10]
# 					rng = self.xlimits
# 				
# 					self.hist,self.hist_edges = np.histogram(x,bins,range=rng,density=normed)
# 					ax.hist(x,bins=self.hist_edges,color=color,range=rng,histtype=histtype,orientation=orient,align=align,alpha=alpha,rwidth=rwidth,normed=normed,cumulative=cumulative,log=log,stacked=False)
# 				
# 					if log: 
# 						if cumulative: 
# 							ymax = max(np.cumsum(np.log10(self.hist)))
# 						else: 
# 							ymax = max(self.hist)
# 					elif cumulative: 
# 						ymax = max(np.cumsum(self.hist))
# 					else: 
# 						ymax = max(self.hist)
# 				
# 					margin = ymax*0.025
# 					ax.set_ylim([0-margin,ymax+margin])
# 					self.ylimits=(0-margin,ymax+margin)
# 				
			canvas.draw()
			self.plotimg = canvas.tostring_rgb()  # save this and convert to bitmap as needed

			# this try except block is because the developers of matplotlib have been changing their API
			try: # this would work for matplotlib 0.98
				self.scrlim=(ax.get_window_extent().xmin,ax.get_window_extent().ymin,ax.get_window_extent().xmax-ax.get_window_extent().xmin,ax.get_window_extent().ymax-ax.get_window_extent().ymin)
			except:
				try: # this should work for matplotlib 0.91
					self.scrlim=(ax.get_window_extent().xmin(),ax.get_window_extent().ymin(),ax.get_window_extent().xmax()-ax.get_window_extent().xmin(),ax.get_window_extent().ymax()-ax.get_window_extent().ymin())
				except:
					print 'there is a problem with your matplotlib'
					return
			self.plotlim=(ax.get_xlim()[0],ax.get_ylim()[0],ax.get_xlim()[1]-ax.get_xlim()[0],ax.get_ylim()[1]-ax.get_ylim()[0])
			
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
			if self.axisparms[2]=="linear" : x2=(x-self.scrlim[0])/self.scrlim[2]*self.plotlim[2]+self.plotlim[0]
			else : x2=10.0**((x-self.scrlim[0])/self.scrlim[2]*(log10(self.plotlim[2]+self.plotlim[0])-log10(self.plotlim[0]))+log10(self.plotlim[0]))
			if self.axisparms[3]=="linear" : y2=(self.height()-y-self.scrlim[1])/self.scrlim[3]*self.plotlim[3]+self.plotlim[1]
			else : y2=10.0**((self.height()-y-self.scrlim[1])/self.scrlim[3]*(log10(self.plotlim[3]+self.plotlim[1])-log10(self.plotlim[1]))+log10(self.plotlim[1]))
			return (x2,y2)
		except: return (0,0)

	def plot2scr(self,x,y) :
		""" converts plot coordinates to screen coordinates """
		try:
			if self.axisparms[2]=="linear" : x2=(x-self.plotlim[0])/self.plotlim[2]*self.scrlim[2]+self.scrlim[0]
			else : x2=(-(self.scrlim[2]*log(x)) + (self.scrlim[0] + self.scrlim[2])*log(10)*log10(self.plotlim[0])-self.scrlim[0]*log(10)*log10(self.plotlim[0] +self.plotlim[2])) /(log(10)*(log10(self.plotlim[0]) - log10(self.plotlim[0] + self.plotlim[2])))
			if self.axisparms[3]=="linear" :y2=self.height()-((y-self.plotlim[1])/self.plotlim[3]*self.scrlim[3]+self.scrlim[1])
			else : y2=(self.scrlim[3]*log(y) + self.height()*log(10.0)*log10(self.plotlim[1])-self.scrlim[1]*log(10.0)*log10(self.plotlim[1])-self.scrlim[3]*log(10.0)*log10(self.plotlim[1]) - self.height()*log(10.0)*log10(self.plotlim[1]+self.plotlim[3]) + self.scrlim[1]*log(10)*log10(self.plotlim[1]+self.plotlim[3])) / (log(10)*(log10(self.plotlim[1]) - log10(self.plotlim[1]+self.plotlim[3])))
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
		#self.axisparms=(xlabel,ylabel,xlog,ylog)
		self.needupd=1
		self.updateGL()
	
	def setAxes(self,key,xa,quiet=False):
		if self.axes[key]==(xa,) : return
		self.axes[key]=(xa,)
		self.autoscale(True)
		self.needupd=1
		if not quiet : self.updateGL()
	
	def setPlotParms(self,key,nbins=None,color=None,htype=None,orient=None,align=None,alpha=None,width=None,norm=None,cumul=None,logy=None,stacked=None,quiet=False):
		if nbins==None : nbins=self.pparm[key][0]
		if color==None : color=self.pparm[key][1]
		if htype==None : htype=self.pparm[key][2]
		if orient==None : orient=self.pparm[key][3]
		if align==None : align=self.pparm[key][4]
		if alpha==None : alpha=self.pparm[key][5]
		if width==None : width=self.pparm[key][6]
		if norm==None : norm=self.pparm[key][7]
		if cumul==None : cumul=self.pparm[key][8]
		if logy==None: logy=self.pparm[key][9]
		if stacked==None: stacked=self.pparm[key][10]
		
		if self.pparm[key]==(nbins,color,htype,orient,align,alpha,width,norm,cumul,logy,stacked) : return
		self.pparm[key]=(nbins,color,htype,orient,align,alpha,width,norm,cumul,logy,stacked)
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

	def add_shapes(self,d):
		self.shapes.update(d)
		self.shapechange=1

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

	def update_selected(self,evc,lc): # could use this to highlight selected bin in histogram
		"""Update the list of 'selected' points, and display coord info.
evc is the cursor selection point in screen coords
lc is the cursor selection point in plot coords"""
		return

	def mousePressEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector(1)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.del_shapes()
			self.updateGL()
			self.rmousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
			self.add_shape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			histlabel = ""
			for i in range(len(self.hist_edges)-1):
				if lc[0] > self.hist_edges[i] and lc[0] < self.hist_edges[i+1]:
					histlabel = "{}; ({:1.5g}, {:1.5g})".format(self.hist[i],self.hist_edges[i],self.hist_edges[i+1])
					break
			self.add_shape("lcrosshist",EMShape(("scrlabel",0,0,0,self.scrlim[2]-175,self.scrlim[3]-10,histlabel,120.0,-1)))
			self.update_selected((event.x(),event.y()),lc)
			self.updateGL()
	
	def mouseMoveEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())
		if  self.rmousedrag:
			self.add_shape("xzoom1",EMShape(("scrline",0,0,0,self.rmousedrag[0],self.scrlim[1],self.rmousedrag[0],self.scrlim[3]+self.scrlim[1],1)))
			self.add_shape("xzoom2",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			zm = self.scr2plot(self.rmousedrag[0],self.rmousedrag[1])
			zoomlabel = "{:1.5g}".format(lc[0]-zm[0])
			self.add_shape("lzoom",EMShape(("scrlabel",0,0,0,self.scrlim[2]-175,self.scrlim[3]-10,zoomlabel,120.0,-1)))
			self.updateGL()
		elif event.buttons()&Qt.LeftButton:
			self.add_shape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			histlabel = ""
			for i in range(len(self.hist_edges)-1):
				if lc[0] > self.hist_edges[i] and lc[0] < self.hist_edges[i+1]:
					histlabel = "{}; ({:1.5g}, {:1.5g})".format(self.hist[i],self.hist_edges[i],self.hist_edges[i+1])
					break
			self.add_shape("lcrosshist",EMShape(("scrlabel",0,0,0,self.scrlim[2]-175,self.scrlim[3]-10,histlabel,120.0,-1)))
			self.update_selected((event.x(),event.y()),lc)
			self.updateGL()
	
	def mouseReleaseEvent(self, event):
		lc =self.scr2plot(event.x(),event.y())
		if self.rmousedrag:
			lc2=self.scr2plot(*self.rmousedrag)
			if fabs(event.x()-self.rmousedrag[0])+fabs(event.y()-self.rmousedrag[1])<3 : self.rescale(0,0,0,0)
			else :
				ymax = max(self.hist)
				ymin = 0
				margin=(ymax-ymin)*0.025
				self.rescale(min(lc[0],lc2[0]),max(lc[0],lc2[0]),ymin-margin,ymax+margin)
			self.rmousedrag=None
	
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
		self.needupd=1
		self.del_shapes()  # also triggers an update
		self.updateGL()
		if self.inspector: self.inspector.update()

	def autoscale(self,force=False):
		"This autoscales, but only axes which currently have invalid settings"
		
		if force or self.xlimits==None or self.xlimits[1]<=self.xlimits[0] :
			xmin=1.0e38
			xmax=-1.0e38
			for k in self.axes.keys():
				if not self.visibility[k]: continue
				xmin=min(xmin,min(self.data[k][self.axes[k][0]]))
				xmax=max(xmax,max(self.data[k][self.axes[k][0]]))

			if self.axisparms[2]!="linear" : self.xlimits=(xmin/1.1,xmax*1.1)
			else:
				margin=(xmax-xmin)*0.025
				self.xlimits=(xmin-margin,xmax+margin)

		if force or self.ylimits==None or self.ylimits[1]<=self.ylimits[0] :
			ymin=0
			ymax=-1.0e38
			for k in self.axes.keys():
				if not self.visibility[k]: continue
				self.hist,self.hist_edges = np.histogram(self.data[k][self.axes[k][0]],self.pparm[k][1])
				ymax=max(self.hist)
			if self.axisparms[3]!="linear" : self.ylimits=(ymin/1.1,ymax*1.1)
			else:
				margin=(ymax-ymin)*0.025
				self.ylimits=(ymin-margin,ymax+margin)

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


class EMHistogramInspector(QtGui.QWidget):

	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"plot.png"))
		self.target=weakref.ref(target)
		
		vbl0=QtGui.QVBoxLayout(self)

		hbl = QtGui.QHBoxLayout()
		hbl.setMargin(2)
		hbl.setSpacing(6)
		hbl.setObjectName("hbl")

		gbx = QtGui.QGroupBox("Data sets")

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
		hbl6 = QtGui.QHBoxLayout()
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
		hbl7 = QtGui.QHBoxLayout()
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

		hbl0=QtGui.QHBoxLayout()
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
		hbl0.addWidget(self.pdfb)

		hbl01=QtGui.QHBoxLayout()
		hbl01.setMargin(0)
		hbl01.setSpacing(6)
		vbl.addLayout(hbl01)

		self.histtype=QtGui.QComboBox(self)
		self.histtype.addItem("bar")
		#self.histtype.addItem("barstacked")
		self.histtype.addItem("step")
		self.histtype.addItem("stepfilled")
		hbl01.addWidget(self.histtype)

		self.orient=QtGui.QComboBox(self)
		self.orient.addItem("vertical")
		#self.orient.addItem("horizontal")
		hbl01.addWidget(self.orient)
		
		self.align=QtGui.QComboBox(self)
		self.align.addItem("mid")
		self.align.addItem("left")
		self.align.addItem("right")
		hbl01.addWidget(self.align)
		
		hbl1 = QtGui.QHBoxLayout()
		hbl1.setMargin(0)
		hbl1.setSpacing(6)
		
		self.color=QtGui.QComboBox(self)
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
		
		hbl000 = QtGui.QHBoxLayout()
		hbl000.setMargin(0)
		hbl000.setSpacing(6)
		
		self.alpha=ValSlider(self,(0,1),"Alpha:",0,25)
		self.alpha.setValue(0.8)
		hbl000.addWidget(self.alpha)
		vbl.addLayout(hbl000)
		
		hbl001 = QtGui.QHBoxLayout()
		hbl001.setMargin(0)
		hbl001.setSpacing(6)
		
		self.rwidth=ValSlider(self,(0,1),"Width:",0,25)
		self.rwidth.setValue(0.8)
		hbl001.addWidget(self.rwidth)
		
		vbl.addLayout(hbl001)
		
		hbl2 = QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(6)
		
		vbl.addLayout(hbl2)
		
		hbl2 = QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(6)
		vbl.addLayout(hbl2)
		
		# per plot column selectors
		gl=QtGui.QGridLayout()
		gl.addWidget(QtGui.QLabel("Column:",self),0,0,Qt.AlignRight)
		self.slidecol=QtGui.QSpinBox(self)
		self.slidecol.setRange(0,1)
		self.slidecol.setValue(1)
		gl.addWidget(self.slidecol,0,1,Qt.AlignLeft)

		gl.addWidget(QtGui.QLabel("N Bins:",self),0,2,Qt.AlignRight)
		self.slidenbs=QtGui.QSpinBox(self)
		self.slidenbs.setRange(1,10000)
		self.slidenbs.setValue(10)
		gl.addWidget(self.slidenbs,0,3,Qt.AlignLeft)
		
		vbl.addLayout(gl)
		
		hbl02=QtGui.QHBoxLayout()
		hbl02.setMargin(0)
		hbl02.setSpacing(6)
		vbl.addLayout(hbl02)

		self.normed=QtGui.QCheckBox(self)
		self.normed.setText("Norm")
		hbl02.addWidget(self.normed)

		self.cumulative=QtGui.QCheckBox(self)
		self.cumulative.setText("Cumulative")
		hbl02.addWidget(self.cumulative)
		
		hbl03=QtGui.QHBoxLayout()
		hbl03.setMargin(0)
		hbl03.setSpacing(6)
		vbl.addLayout(hbl03)
		
		self.logtogy=QtGui.QCheckBox(self)
		self.logtogy.setText("Log Y")
		hbl03.addWidget(self.logtogy)
		
		self.stacked=QtGui.QCheckBox(self)
		self.stacked.setText("Stacked")
		hbl03.addWidget(self.stacked)

		self.wrescale=QtGui.QPushButton(self)
		self.wrescale.setText("Rescale")
		vbl.addWidget(self.wrescale)
		
		vbl0.addLayout(hbl)

		hbl2a=QtGui.QHBoxLayout()

		self.wl1=QtGui.QLabel("Min")
		self.wl1.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl1)
		self.wl2=QtGui.QLabel("Max")
		self.wl2.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl2)
		self.wl3=QtGui.QLabel("Min")
		self.wl3.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl3)
		self.wl4=QtGui.QLabel("Max")
		self.wl4.setAlignment(Qt.AlignHCenter)
		hbl2a.addWidget(self.wl4)
		vbl0.addLayout(hbl2a)

		hbl2=QtGui.QHBoxLayout()

		self.wxmin=ValBox(label="X:")
		hbl2.addWidget(self.wxmin)
		self.wxmax=ValBox(label="  ")
		hbl2.addWidget(self.wxmax)

		self.wymin=ValBox(label="Y:")
		hbl2.addWidget(self.wymin)
		self.wymax=ValBox(label="  ")
		hbl2.addWidget(self.wymax)

		vbl0.addLayout(hbl2)
# 
# 		hbl4 = QtGui.QHBoxLayout()
# 		hbl4.addWidget(QtGui.QLabel("X Label:",self))
# 		self.xlabel=QtGui.QLineEdit(self)
# 		hbl4.addWidget(self.xlabel)
# 		vbl0.addLayout(hbl4)
# 
# 		hbl5 = QtGui.QHBoxLayout()
# 		hbl5.addWidget(QtGui.QLabel("Y Label:",self))
# 		self.ylabel=QtGui.QLineEdit(self)
# 		hbl5.addWidget(self.ylabel)
# 		vbl0.addLayout(hbl5)

		self.quiet=0
		self.busy=0

		QtCore.QObject.connect(self.showslide, QtCore.SIGNAL("valueChanged"), self.selSlide)
		QtCore.QObject.connect(self.allbut, QtCore.SIGNAL("clicked()"), self.selAll)
		QtCore.QObject.connect(self.nonebut, QtCore.SIGNAL("clicked()"), self.selNone)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.list_item_changed)
		QtCore.QObject.connect(self.saveb,QtCore.SIGNAL("clicked()"),self.savePlot)
		QtCore.QObject.connect(self.pdfb,QtCore.SIGNAL("clicked()"),self.savePdf)
		QtCore.QObject.connect(self.concatb,QtCore.SIGNAL("clicked()"),self.saveConcatPlot)
		QtCore.QObject.connect(self.normed,QtCore.SIGNAL("stateChanged(int)"),self.updPlot)
		QtCore.QObject.connect(self.logtogy,QtCore.SIGNAL("stateChanged(int)"),self.updPlot)
		QtCore.QObject.connect(self.cumulative,QtCore.SIGNAL("stateChanged(int)"),self.updPlot)
		QtCore.QObject.connect(self.stacked,QtCore.SIGNAL("stateChanged(int)"),self.updPlot)
		QtCore.QObject.connect(self.slidecol, QtCore.SIGNAL("valueChanged(int)"), self.newCols)
		QtCore.QObject.connect(self.slidenbs, QtCore.SIGNAL("valueChanged(int)"), self.updPlot)
		QtCore.QObject.connect(self.rwidth,QtCore.SIGNAL("valueChanged"),self.updPlot)
		QtCore.QObject.connect(self.alpha,QtCore.SIGNAL("valueChanged"),self.updPlot)
		QtCore.QObject.connect(self.color,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlotColor)
		QtCore.QObject.connect(self.histtype,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.orient,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.align,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		#QtCore.QObject.connect(self.xlabel,QtCore.SIGNAL("textChanged(QString)"),self.updPlot)
		#QtCore.QObject.connect(self.ylabel,QtCore.SIGNAL("textChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.wxmin,QtCore.SIGNAL("valueChanged"),self.newLimits)
		QtCore.QObject.connect(self.wxmax,QtCore.SIGNAL("valueChanged"),self.newLimits)
		QtCore.QObject.connect(self.wymin,QtCore.SIGNAL("valueChanged"),self.newLimits)
		QtCore.QObject.connect(self.wymax,QtCore.SIGNAL("valueChanged"),self.newLimits)
		QtCore.QObject.connect(self.wrescale,QtCore.SIGNAL("clicked()"),self.autoScale)

		self.newSet(0)
		self.datachange()
		#self.resize(500,540)

	def selSlide(self,val):
		rngn0=int(val)
		rngn1=int(self.nbox.getValue())
		rngstp=int(self.stepbox.getValue())
		rng=range(rngn0,rngn0+rngstp*rngn1,rngstp)
		for i,k in enumerate(sorted(self.target().visibility.keys())) :
			self.target().visibility[k]=i in rng
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def selAll(self):
		for k in self.target().visibility.keys() : self.target().visibility[k]=True
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def selNone(self):
		for k in self.target().visibility.keys() : self.target().visibility[k]=False
		self.target().full_refresh()
		self.target().updateGL()
		self.datachange()

	def __remove_from_file(self,fname,names):
		if os.path.exists(fname):
			inf = open(fname,"r")
			lines = inf.readlines()
			inf.close()
		else:
			return
		f = open(fname,"w")
		for i in xrange(0,len(lines)):
			lines[i] = lines[i].strip()
		for i in xrange(len(lines)-1,-1,-1):
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
		for i in xrange(0,len(lines)):
			lines[i] = lines[i].strip()
		for name in names:
			if name not in lines:
				lines.append(name)
		for line in lines:
			f.write(line+"\n")
		f.close()

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
		out=file(name2,"a")
		for name in names :
			data=self.target().data[name]
			for i in xrange(len(data[0])):
				out.write("%g\t%g\n"%(data[0][i],data[1][i]))
		out=None
		print "Wrote ",name2

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
			out=file(name2,"w")
			for i in xrange(len(data[0])):
				out.write("%g\t%g\n"%(data[0][i],data[1][i]))
			print "Wrote ",name2

	def savePdf(self):
		"""Saves the contents of the current plot to a pdf"""
		plt.savefig("plot.pdf")

	def updPlot(self,s=None):
		"""
		self.setlist=DragListWidget(self)
		self.nonebut=QtGui.QPushButton(self)
		self.allbut=QtGui.QPushButton(self)
		self.showslide=ValSlider(self,(0,5),"Sel1:",0,30)
		self.nbox=ValBox(label="ns:",value=1)
		self.stepbox=ValBox(label="stp:",value=1)
		self.saveb=QtGui.QPushButton(self)
		self.concatb=QtGui.QPushButton(self)
		self.pdfb=QtGui.QPushButton(self)
		self.histtype=QtGui.QComboBox(self)
		self.orient=QtGui.QComboBox(self)
		self.align=QtGui.QComboBox(self)
		self.color=QtGui.QComboBox(self)
		self.alpha=ValSlider(self,(0,1),"Opacity:",0,25)
		self.rwidth=ValSlider(self,(0,1),"Bin Width:",0,25)
		self.slidecol=QtGui.QSpinBox(self)
		self.slidenbs=QtGui.QSpinBox(self)
		self.normed=QtGui.QCheckBox(self)
		self.cumulative=QtGui.QCheckBox(self)
		self.logtogy=QtGui.QCheckBox(self)
		self.stacked=QtGui.QCheckBox(self)
		self.wrescale=QtGui.QPushButton(self)
		self.wxmin=ValBox(label="X:")
		self.wxmax=ValBox(label="  ")
		self.wymin=ValBox(label="Y:")
		self.wymax=ValBox(label="  ")
		"""
		if self.quiet : return
		#if self.xlogtog.isChecked() : xl="log"
		xl="linear"
		if self.logtogy.isChecked() : yl="log"
		else : yl="linear"
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		#self.target().setAxisParms(self.xlabel.text(),self.ylabel.text(),xl,yl)
		self.target().autoscale(True)
		if len(names)==1:
			ht = histtypes[self.histtype.currentIndex()]
			ori = orientations[self.orient.currentIndex()]
			ali = alignments[self.align.currentIndex()]
			self.target().setPlotParms(names[0],self.slidenbs.value(),self.color.currentIndex(),ht,ori,ali,self.alpha.getValue(),self.rwidth.getValue(),self.normed.isChecked(),self.cumulative.isChecked(),self.logtogy.isChecked(),self.stacked.isChecked(),False)
		else:
			for name in names:
				ht = histtypes[self.histtype.currentIndex()]
				ori = orientations[self.orient.currentIndex()]
				ali = alignments[self.align.currentIndex()]
				self.target().setPlotParms(name,self.slidenbs.value(),self.color.currentIndex(),ht,ori,ali,self.alpha.getValue(),self.rwidth.getValue(),self.normed.isChecked(),self.cumulative.isChecked(),self.logtogy.isChecked(),self.stacked.isChecked(),False)
			self.target().updateGL()
	
	def updPlotColor(self,s=None):
		if self.quiet : return
		names = [str(item.text()) for item in self.setlist.selectedItems()]
		if len(names)==1:
			self.target().setPlotParms(names[0],None,self.color.currentIndex(),None,None,None,self.alpha.getValue(),None,None,None,None,None,False)
		else:
			for name in names:
				self.target().setPlotParms(name,None,self.color.currentIndex(),None,None,None,self.alpha.getValue(),None,None,None,None,None,True)
			self.target().updateGL()

	def newSet(self,row):
		self.quiet=1
		try: i=str(self.setlist.item(row).text())
		except: return
		self.slidecol.setRange(-1,len(self.target().data[i])-1)
		self.slidecol.setValue(self.target().axes[i][0])
		pp=self.target().pparm[i]
		self.color.setCurrentIndex(pp[1])
		self.quiet=0

	def newCols(self,val):
		if self.target and not self.quiet:
			names = [str(item.text()) for item in self.setlist.selectedItems()]
			for name in names:
				self.target().setAxes(name,self.slidecol.value(),True)
			self.target().updateGL()

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

	def autoScale(self):
		self.target().autoscale(True)
		self.target().rescale(0,0,0,0)

	def datachange(self):
		self.setlist.clear()
		flags= Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)|Qt.ItemFlags(Qt.ItemIsDragEnabled)
		keys=self.target().data.keys()
		visible = self.target().visibility
		keys.sort()
		parms = self.target().pparm # get the colors from this
		for i,j in enumerate(keys) :
			a = QtGui.QListWidgetItem(j)
			a.setFlags(flags)
			try: a.setTextColor(qt_color_map[colortypes[parms[j][1]]])
			except:
				print "Color error"
				print list(sorted(parms.keys()))
				print parms[j][1]
				print colortypes[parms[j][1]]
				print qt_color_map[colortypes[parms[j][1]]]
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
		pass

class DragListWidget(QtGui.QListWidget):
	"This is a minor modification of the QListWidget to support drag-drop of data sets"
	def setDataSource(self,trg):
		"""We keep a weak reference to our data source so we can pull the data only when dragging actually starts"""
		self.datasource=weakref.ref(trg)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_Backspace:
			name=str(self.currentItem().text())		# currently hilighted item
			nbins=len(self.datasource().target().hist_edges)
			self.datasource().target().set_data(None,key=name,alpha=0.8,width=0.8,nbins=nbins)
		else: QtGui.QListWidget.keyPressEvent(self,event)

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
					data=[ [] for i in xrange(n)]		# initialize empty data arrays
				# parses out each number from each line and puts it in our list of lists
				for i,f in enumerate(rex.findall(s)):
					try: data[i].append(float(f))
					except: print "Error (%d): %s"%(i,f)
			# Find an unused name for the data set
			trgplot=self.datasource().target()
			name="Dropped"
			nn=1
			while trgplot.data.has_key(name) :
				name="Dropped_%d"%nn
				nn+=1
			trgplot.set_data(data,name,quiet=True)
			if n==1: trgplot.setAxes(name,0)
			e.acceptProposedAction()

	def supportedDropActions(self):
		return Qt.DropActions(Qt.CopyAction)

	def setMovement(self,x):
		"""The ListView and ListWidget unfortunately make use of drag-drop for internal rearrangement, but we need to use it for widget->widget copy. This prevents the parent from disabling drag/drop."""
		QtGui.QListWidget.setMovement(self,x)
		self.setlist.setDragEnabled(True)
		self.setlist.setAcceptDrops(True)

	def setViewMode(self,x):
		"""The ListView and ListWidget unfortunately make use of drag-drop for internal rearrangement, but we need to use it for widget->widget copy. This prevents the parent from disabling drag/drop."""
		QtGui.QListWidget.setViewMode(self,x)
		self.setlist.setDragEnabled(True)
		self.setlist.setAcceptDrops(True)

	def startDrag(self,actions):
		data,axes,pparm=self.datasource().getCurrentData()		# get the "current" data object
		if data==None : return						# don't start a drag if nothing is selected
		# we only copy the specific axes that are used in the current plot settings !
		if axes[0]<0: axes=[axes[0]]
		## create the string representation of the data set
		sdata=StringIO()		# easier to write as if to a file
		for y in xrange(len(data[0])):
			sdata.write("%1.8g"%data[axes[0]][y])
			for x in xrange(1,len(axes)):
				sdata.write("\t%1.8g"%data[axes[x]][y])
			sdata.write("\n")
		# start the drag operation
		drag = QtGui.QDrag(self)
		mimeData = QtCore.QMimeData()
		mimeData.setText(sdata.getvalue())
		drag.setMimeData(mimeData)
		dropact = drag.exec_(Qt.CopyAction)


# if __name__ == '__main__': # This is just for testing, of course
# 	app = EMApp()
# 	window = EMHistogramWidget(app)
# 	if len(sys.argv)==1 :
# 		l=[i/30.*pi for i in range(30)]
# 		window.set_data([[1,2,3,4],[2,3,4,3]],"test")
# 		window.set_data([l,[sin(2*i) for i in l]],"test2")
# 	else:
# 		for i in range(1,len(sys.argv)):
# 			window.set_data_from_file(sys.argv[i])
# 
# 	app.show()
# 	app.execute()

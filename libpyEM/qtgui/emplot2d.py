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

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from valslider import ValSlider
from math import *
from EMAN2 import *
import EMAN2
import sys
import numpy
from emshape import *
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
import struct

import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
#matplotlib.use('Agg')

linetypes=["-","--",":","-."]
symtypes=["o","s","+","2","1"]
colortypes=["k","b","r","g"]

def NewPlot2DWin():
	return EMParentWin(EMPlot2D())


class EMPlot2D(QtOpenGL.QGLWidget):
	"""A QT widget for drawing 2-D plots using matplotlib
	"""
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		self.axes={}
		self.pparm={}
		self.inspector=None
		self.needupd=1
		self.plotimg=None
		self.shapes={}
		self.limits=None
		self.rmousedrag=None

		self.data={}				# List of Lists to plot 
		
	def setData(self,key,data):
		"""Set a keyed data set. The key should generally be a string describing the data.
		'data' is a tuple/list of tuples/list representing all values for a particular
		axis. eg - the points: 1,5; 2,7; 3,9 would be represented as ((1,2,3),(5,7,9)).
		Multiple axes may be set, and which axis represents which axis in the plot can be
		selected by the user. 'data' can also be an EMData object, in which case the entire
		data array is plotted as if it were 1-D."""
		
		self.needupd=1
		
		if type(data)==EMData :
			data=(data.get_data_pickle(),)
		
		try:
			if len(data)>1 : self.axes[key]=(0,1,-1)
			else : self.axes[key]=(-1,0,-1)
		except: return
		
		self.pparm[key]=(0,1,0,1,0,0,5)
		
		if data : self.data[key]=data
		else : del self.data[key]
		
		
		if self.inspector: self.inspector.datachange()
		
		self.updateGL()
		
	def setDataFromFile(self,key,filename):
		"""Reads a keyed data set from a file. Automatically interpret the file contents."""
		
		data=None
		try:
			# if it's an image file
			im=EMData(filename)
			# This could be faster...
			data=[[im.get_value_at(i,j) for j in range(im.get_ysize())] for j in range(im.get_xsize())]
		except:
			try:
				# Maybe an EMAN1 .fp file
				fin=file(filename)
				fph=struct.unpack("120sII",fin.read(128))
				ny=fph[1]
				nx=fph[2]
				data=[]
				for i in range(nx):
					data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
			except:
				# Probably a text file
				fin=file(filename)
				fin.seek(0)
				rdata=fin.readlines()
				rdata=[i for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				
				data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]
				
		if data : self.setData(filename.split("/")[-1],data)
		
	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
		
	def paintGL(self):
		if not self.parentWidget() : return
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		
		if self.needupd or not self.plotimg:
			self.needupd=0
			fig=Figure((self.width()/72.0,self.height()/72.0),dpi=72.0)
			if self.limits :ax=fig.add_axes((.1,.05,.85,.9),autoscale_on=False,xlim=self.limits[0],ylim=self.limits[1])
			else : ax=fig.add_axes((.1,.05,.85,.9),autoscale_on=True)
			canvas=FigureCanvasAgg(fig)
			
			for i in self.axes.keys():
				j=self.axes[i]
				if j[0]==-1 : x=range(len(self.data[i][0]))
				else : x=self.data[i][self.axes[i][0]]
				if j[1]==-1 : y=range(len(self.data[i][0]))
				else : y=self.data[i][self.axes[i][1]]
				parm=""
				parm+=colortypes[self.pparm[i][0]]
				if self.pparm[i][1]: 
					parm+=linetypes[self.pparm[i][2]]
				if self.pparm[i][4]:
					parm+=symtypes[self.pparm[i][5]]
					
				ax.plot(x,y,parm,linewidth=self.pparm[i][3],markersize=self.pparm[i][6])
			
			
			canvas.draw()
			self.plotimg = canvas.tostring_rgb()  # save this and convert to bitmap as needed
			
			self.scrlim=(ax.get_window_extent().xmin(),ax.get_window_extent().ymin(),ax.get_window_extent().xmax()-ax.get_window_extent().xmin(),ax.get_window_extent().ymax()-ax.get_window_extent().ymin())
			self.plotlim=(ax.get_xlim()[0],ax.get_ylim()[0],ax.get_xlim()[1]-ax.get_xlim()[0],ax.get_ylim()[1]-ax.get_ylim()[0])
		
#		print ax.get_window_extent().xmin(),ax.get_window_extent().ymin()
#		print ax.get_window_extent().xmax(),ax.get_window_extent().ymax()
#		print ax.get_position()
		
		GL.glRasterPos(0,self.height()-1)
		GL.glPixelZoom(1.0,-1.0)
#		print "paint ",self.width(),self.height(), self.width()*self.height(),len(a)
		GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1)
		GL.glDrawPixels(self.width(),self.height(),GL.GL_RGB,GL.GL_UNSIGNED_BYTE,self.plotimg)
		
		GL.glPushMatrix()
		for k,s in self.shapes.items():
			s.draw(self.scr2plot)
		GL.glPopMatrix()


	def scr2plot(self,x,y) :
		""" converts screen coordinates to plot coordinates """
		try: 
			return ((x-self.scrlim[0])/self.scrlim[2]*self.plotlim[2]+self.plotlim[0],(self.height()-y-self.scrlim[1])/self.scrlim[3]*self.plotlim[3]+self.plotlim[1])
		except: return (0,0)
		
	def plot2scr(self,x,y) :
		""" converts plot coordinates to screen coordinates """
		try: 
			return ((x-self.plotlim[0])/self.plotlim[2]*self.scrlim[2]+self.scrlim[0],(self.height()-y-self.plotlim[1])/self.plotlim[3]*self.scrlim[3]+self.scrlim[1])
		except: return (0,0)


	def resizeGL(self, width, height):
#		print "resize ",self.width()
		self.needupd=1
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())
	
		self.delShapes(("xcross","ycross","lcross"))
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMPlot2DInspector(self)
		self.inspector.show()
		self.inspector.datachange()
		self.needupd=1
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
	
	def setAxes(self,key,xa,ya,za):
		if self.axes[key]==(xa,ya,za) : return
		self.axes[key]=(xa,ya,za)
		self.needupd=1
		self.updateGL()
		
	def setPlotParms(self,key,color,line,linetype,linewidth,sym,symtype,symsize):
		self.pparm[key]=(color,line,linetype,linewidth,sym,symtype,symsize)
		self.needupd=1
		self.updateGL()
	
	def addShape(self,k,s):
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
		self.updateGL()
	
	def addShapes(self,d):
		self.shapes.update(d)
		self.shapechange=1
		self.updateGL()
	
	def delShapes(self,k=None):
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
		self.updateGL()
		
	#def dragEnterEvent(self,event):
		
		#if event.provides("application/x-eman"):
			#event.setDropAction(Qt.CopyAction)
			#event.accept()

	
	#def dropEvent(self,event):
#=		if EMAN2.GUIbeingdragged:
			#self.setData(EMAN2.GUIbeingdragged)
			#EMAN2.GUIbeingdragged=None
		#elif event.provides("application/x-eman"):
			#x=loads(event.mimeData().data("application/x-eman"))
			#self.setData(x)
			#event.acceptProposedAction()

	
	def mousePressEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):
			self.showInspector(1)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.delShapes()
			self.rmousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
			self.addShape("xcross",EMShape(("scrline",0,0,0,self.scrlim[0],self.height()-event.y(),self.scrlim[2]+self.scrlim[0],self.height()-event.y(),1)))
			self.addShape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			#if self.mmode==0:
				#self.emit(QtCore.SIGNAL("mousedown"), event)
				#return
			#elif self.mmode==1 :
				#try: 
					#del self.shapes["MEASL"]
				#except: pass
				#self.addShape("MEAS",("line",.5,.1,.5,lc[0],lc[1],lc[0]+1,lc[1],2))
	
	def mouseMoveEvent(self, event):
		lc=self.scr2plot(event.x(),event.y())
		
		if  self.rmousedrag:
			self.addShape("zoom",EMShape(("scrrect",0,0,0,self.rmousedrag[0],self.height()-self.rmousedrag[1],event.x(),self.height()-event.y(),1)))
		elif event.buttons()&Qt.LeftButton:
			self.addShape("xcross",EMShape(("scrline",0,0,0,self.scrlim[0],self.height()-event.y(),self.scrlim[2]+self.scrlim[0],self.height()-event.y(),1)))
			self.addShape("ycross",EMShape(("scrline",0,0,0,event.x(),self.scrlim[1],event.x(),self.scrlim[3]+self.scrlim[1],1)))
			self.addShape("lcross",EMShape(("scrlabel",0,0,0,self.scrlim[2]-100,self.scrlim[3]-10,"%1.4g, %1.4g"%(lc[0],lc[1]),1.5,1)))
#			self.addShape("mcross",EMShape(("scrlabel",0,0,0,self.scrlim[2]-80,self.scrlim[3]-20,"%1.3g, %1.3g"%(self.plot2scr(*lc)[0],self.plot2scr(*lc)[1]),1.5,1)))
	
	def mouseReleaseEvent(self, event):
		lc =self.scr2plot(event.x(),event.y())
		if self.rmousedrag:
			lc2=self.scr2plot(*self.rmousedrag)
			if fabs(event.x()-self.rmousedrag[0])+fabs(event.y()-self.rmousedrag[1])<3 : self.limits=None
			else : self.limits=((min(lc[0],lc2[0]),max(lc[0],lc2[0])),(min(lc[1],lc2[1]),max(lc[1],lc2[1])))
			self.rmousedrag=None
			self.needupd=1
			self.delShapes()  # also triggers an update
		#elif event.button()==Qt.LeftButton:
			#if self.mmode==0:
				#self.emit(QtCore.SIGNAL("mouseup"), event)
				#return
			#elif self.mmode==1 :
				#self.addShape("MEAS",("line",.5,.1,.5,self.shapes["MEAS"][4],self.shapes["MEAS"][5],lc[0],lc[1],2))

	def wheelEvent(self, event):
		pass
		#if event.delta() > 0:
			#self.setScale(self.scale + MAG_INC)	
		#elif event.delta() < 0:
			#if ( self.scale - MAG_INC > 0 ):
				#self.setScale(self.scale - MAG_INC)
		## The self.scale variable is updated now, so just update with that
		#if self.inspector: self.inspector.setScale(self.scale)

class EMPlot2DInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(2)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		# plot list
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.hbl.addWidget(self.setlist)
		
		self.vbl = QtGui.QVBoxLayout()
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		self.hbl.addLayout(self.vbl)
		
		self.color=QtGui.QComboBox(self)
		self.color.addItem("black")
		self.color.addItem("blue")
		self.color.addItem("red")
		self.color.addItem("green")
		self.vbl.addWidget(self.color)

		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.vbl.addLayout(self.hbl2)
				
		# This is for line parms
		self.vbl2b = QtGui.QVBoxLayout()
		self.vbl2b.setMargin(0)
		self.vbl2b.setSpacing(6)
		self.hbl2.addLayout(self.vbl2b)
				
		self.lintog=QtGui.QPushButton(self)
		self.lintog.setText("Line")
		self.lintog.setCheckable(1)
		self.vbl2b.addWidget(self.lintog)
				
		self.linsel=QtGui.QComboBox(self)
		self.linsel.addItem("------")
		self.linsel.addItem("- - - -")
		self.linsel.addItem(".......")
		self.linsel.addItem("-.-.-.-")
		self.vbl2b.addWidget(self.linsel)
		
		self.linwid=QtGui.QSpinBox(self)
		self.linwid.setRange(1,10)
		self.vbl2b.addWidget(self.linwid)
		
		# This is for point parms
		self.vbl2a = QtGui.QVBoxLayout()
		self.vbl2a.setMargin(0)
		self.vbl2a.setSpacing(6)
		self.hbl2.addLayout(self.vbl2a)
				
		self.symtog=QtGui.QPushButton(self)
		self.symtog.setText("Symbol")
		self.symtog.setCheckable(1)
		self.vbl2a.addWidget(self.symtog)

		self.symsel=QtGui.QComboBox(self)
		self.symsel.addItem("circle")
		self.symsel.addItem("square")
		self.symsel.addItem("plus")
		self.symsel.addItem("triup")
		self.symsel.addItem("tridown")
		self.vbl2a.addWidget(self.symsel)
		
		self.symsize=QtGui.QSpinBox(self)
		self.symsize.setRange(0,25)
		self.vbl2a.addWidget(self.symsize)
		
		# per plot column selectors
		self.slidex=QtGui.QSpinBox(self)
		self.slidex.setRange(-1,1)
		self.vbl.addWidget(self.slidex)
		#self.slidex=ValSlider(self,(-1,1),"X col:",0)
		#self.slidex.setIntonly(1)
		#self.vbl.addWidget(self.slidex)
		
		self.slidey=QtGui.QSpinBox(self)
		self.slidey.setRange(-1,1)
		self.vbl.addWidget(self.slidey)
		#self.slidey=ValSlider(self,(-1,1),"Y col:",1)
		#self.slidey.setIntonly(1)
		#self.vbl.addWidget(self.slidey)
		
		self.slidec=QtGui.QSpinBox(self)
		self.slidec.setRange(-1,1)
		self.vbl.addWidget(self.slidec)
		#self.slidec=ValSlider(self,(-1,1),"C col:",-1)
		#self.slidec.setIntonly(1)
		#self.vbl.addWidget(self.slidec)
		
		self.setLayout(self.hbl)

		self.quiet=0
		
		QtCore.QObject.connect(self.slidex, QtCore.SIGNAL("valueChanged(int)"), self.newCols)
		QtCore.QObject.connect(self.slidey, QtCore.SIGNAL("valueChanged(int)"), self.newCols)
		QtCore.QObject.connect(self.slidec, QtCore.SIGNAL("valueChanged(int)"), self.newCols)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.color,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.symtog,QtCore.SIGNAL("clicked()"),self.updPlot)
		QtCore.QObject.connect(self.symsel,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.symsize,QtCore.SIGNAL("valueChanged(int)"),self.updPlot)
		QtCore.QObject.connect(self.lintog,QtCore.SIGNAL("clicked()"),self.updPlot)
		QtCore.QObject.connect(self.linsel,QtCore.SIGNAL("currentIndexChanged(QString)"),self.updPlot)
		QtCore.QObject.connect(self.linwid,QtCore.SIGNAL("valueChanged(int)"),self.updPlot)
		self.datachange()
		
		
		#QtCore.QObject.connect(self.gammas, QtCore.SIGNAL("valueChanged"), self.newGamma)
		#QtCore.QObject.connect(self.invtog, QtCore.SIGNAL("toggled(bool)"), target.setInvert)
		#QtCore.QObject.connect(self.mmode, QtCore.SIGNAL("buttonClicked(int)"), target.setMMode)

	def updPlot(self,s=None):
		if self.quiet : return
		self.target.setPlotParms(str(self.setlist.currentItem().text()),self.color.currentIndex(),self.lintog.isChecked(),
			self.linsel.currentIndex(),self.linwid.value(),self.symtog.isChecked(),
			self.symsel.currentIndex(),self.symsize.value())

	def newSet(self,row):
		self.quiet=1
		i=str(self.setlist.item(row).text())
		self.slidex.setRange(-1,len(self.target.data[i])-1)
		self.slidey.setRange(-1,len(self.target.data[i])-1)
		self.slidec.setRange(-1,len(self.target.data[i])-1)
		self.slidex.setValue(self.target.axes[i][0])
		self.slidey.setValue(self.target.axes[i][1])
		self.slidec.setValue(self.target.axes[i][2])
		
		pp=self.target.pparm[i]
		self.color.setCurrentIndex(pp[0])
		
		self.lintog.setChecked(pp[1])
		self.linsel.setCurrentIndex(pp[2])
		self.linwid.setValue(pp[3])
		
		self.symtog.setChecked(pp[4])
		self.symsel.setCurrentIndex(pp[5])
		self.symsize.setValue(pp[6])
		self.quiet=0

	def newCols(self,val):
		if self.target: 
			self.target.setAxes(str(self.setlist.currentItem().text()),self.slidex.value(),self.slidey.value(),self.slidec.value())
	
	def datachange(self):
		
		self.setlist.clear()
		
		keys=self.target.data.keys()
		keys.sort()
		
		for i,j in enumerate(keys) :
			self.setlist.addItem(j)
		
		
# This is just for testing, of course
if __name__ == '__main__':
	GLUT.glutInit(sys.argv )
	app = QtGui.QApplication(sys.argv)
	window = EMPlot2D()
	if len(sys.argv)==1 : 
		l=[i/30.*pi for i in range(30)]
		window.setData("test",[[1,2,3,4],[2,3,4,3],[3,4,5,2]])
		window.setData("test2",[l,[sin(i) for i in l],[cos(i) for i in l]])
	else:
		try:
			# if it's an image file
			im=EMData(sys.argv[1])
			# This could be faster...
			data=[[im.get_value_at(i,j) for j in range(im.get_ysize())] for j in range(im.get_xsize())]
		except:
			try:
				# Maybe a .fp file
				fin=file(sys.argv[1])
				fph=struct.unpack("120sII",fin.read(128))
				ny=fph[1]
				nx=fph[2]
				data=[]
				for i in range(nx):
					data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
			except:
				# Probably a text file
				fin=file(sys.argv[1])
				fin.seek(0)
				rdata=fin.readlines()
				rdata=[i for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				
				print nx,ny
				data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]
				
		window.setData(sys.argv[1].split("/")[-1],data)
			
	window2=EMParentWin(window)
	window2.show()
		
	sys.exit(app.exec_())

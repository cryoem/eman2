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
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

t3d_stack = []

class EMImage3D(QtOpenGL.QGLWidget):
	""" This class is not yet complete.
	A QT widget for rendering 3D EMData objects.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage3D.allim[self]=0
# 		try: 
# 			if EMImage2D.gq : pass
# 		except:
# 			EMImage2D.gq=GLU.gluNewQuadric()
# 			GLU.gluQuadricDrawStyle(EMImage2D.gq,GLU.GLU_FILL)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_SMOOTH)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_NONE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_OUTSIDE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_INSIDE)
# 			GLU.gluQuadricTexture(EMImage2D.gq,GL.GL_FALSE)
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		
		self.data=None
		
		self.aspect=1.0
		self.gq=0
		self.mmode=0
		self.isothr=0.5
		self.isorender=None
		self.isodl = 0
		
		self.scale = 1.0
		
		self.inspector=None
		
		if image :
			print "We have an image"
			self.setData(image)
			self.show()
		
	
	def timeout(self):
		self.updateGL()
		
	
	def setData(self,data):
		"""Pass in a 3D EMData object"""
#		if not self.data and data: self.resize(data.get_xsize(),data.get_ysize())
		
		self.data=data
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			self.updateGL()
			return
		
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma = data.get_attr("sigma")
		
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		
		
		self.start_z = -1.25*data.get_zsize()
		
		#a=self.data.render_amp8()
		#hist=numpy.fromstring(a[-1024:],'i')
		
		#self.inspector.setHist(hist,self.minden,self.maxden) 
		
		self.inspector.setThrs(m0,m1,mean+sigma)
		self.isothr = mean+sigma
		
		self.isorender=MarchingCubes(data,1)
		self.updateGL()
		#self.timer.start(25)
		
	def initializeGL(self):
		glEnable(GL_NORMALIZE)
		
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])

		GL.glClearColor(0,0,0,0)
		
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		# Precompile a displaylist for the display volume border
		self.volcubedl=glGenLists(1)
		glNewList(self.volcubedl,GL_COMPILE)
		glPushMatrix()
		glColor(.7,.7,1.0)
#		glRotate(90.,1.,0.,0.)
		glTranslate(-self.aspect-.01,1.01,-4.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(self.aspect*2.0+.02,0.0,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(0.0,-2.02,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(-self.aspect*2.0-.02,0.0,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)		
		glPopMatrix()
		glEndList()
		
		cl = 50
		self.cylinderdl = glGenLists(1)
		glNewList(self.cylinderdl,GL_COMPILE)
		glColor(.7,.7,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, [.2,.1,0.4,1.0])
		glMaterial(GL_FRONT, GL_SPECULAR, [.2,.2,0.1,1.0])
		glMaterial(GL_FRONT, GL_SHININESS, 32)
		gluCylinder(self.gq,5,5,cl,16,16)
		glEndList()
		
		self.xshapedl = glGenLists(1)
		glNewList(self.xshapedl,GL_COMPILE)
		glPushMatrix()
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glPushMatrix()
		glRotate(90,0,1,0)
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glPushMatrix()
		glRotate(90,1,0,0)
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glEndList()
		
		t3d = Transform3D()
		rot = {}
		rot["az"] = 0.01
		rot["alt"] = 0.01
		rot["phi"] = 0.01
		
		t3d.set_rotation( EULER_EMAN, rot )
		self.t3d_stack = []
		self.t3d_stack.append(t3d)
	
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		# because the viewing volume extends from -0.001 to -10000 in the z direction,
		# translate the scene back into it so it is visible...
		glTranslated(0, 0, self.start_z)

		# get the current rotation from the rotation stack and apply
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		self.updateInspector(rot)
				
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		# here is where zoom is effectively applied
		glScalef(self.scale,self.scale,self.scale)
		
		if self.isodl == 0:
			self.createIsoDL()
		
		glCallList(self.isodl)
		
	def createIsoDL(self):
		# create the isosurface display list
		
		if ( self.isodl != 0 ): glDeleteLists(self.isodl,1)
		
		self.isorender.set_surface_value(self.isothr)
		a=self.isorender.get_isosurface(True)
		f=a["faces"]
		n=a["normals"]
		p=a["points"]
		
		f=[i/3 for i in f]
		n=[i/3 for i in n]
		
		self.isodl = glGenLists(1)
		glNewList(self.isodl,GL_COMPILE)
		glPushMatrix()
		glScalef(self.data.get_xsize(), self.data.get_ysize(), self.data.get_zsize())
		glTranslate(-0.5, -0.5, -0.5)
		glColor(.2,.1,0.4,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, [.2,.1,0.4,1.0])
		glMaterial(GL_FRONT, GL_SPECULAR, [.2,.2,0.1,1.0])
		glMaterial(GL_FRONT, GL_SHININESS, 32)
		glBegin(GL_TRIANGLES)
		for i in f:
			glVertex(p[i*3],p[i*3+1],p[i*3+2])
			glNormal(n[i*3],n[i*3+1],n[i*3+2])
			#print p[i*3],p[i*3+1],p[i*3+2]
		glEnd()
		glPopMatrix()
		glEndList()
		
		print "created isosurface"
		
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(50,aspect,0.001,1000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
	def setupShapes(self):
		# make our own cirle rather than use gluDisk or somesuch
		pass
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector3D(self)
		self.inspector.show()
	
	def updateInspector(self,rot):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		self.inspector.newAz(float(rot["az"])+180)
		self.inspector.newAlt(rot["alt"])
		self.inspector.newPhi(float(rot["phi"])+180)
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
			
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.rpressx = event.x()
				self.rpressy = event.y()
				
				tmp = self.t3d_stack.pop()
				t3d = Transform3D(tmp)
				self.t3d_stack.append(tmp)
				self.t3d_stack.append(t3d)
				
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
		
	def mouseMoveEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
# 		if self.rmousedrag and event.buttons()&Qt.RightButton:
# 			self.origin=(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
# 			self.rmousedrag=(event.x(),event.y())
# 			self.update()
		if event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
		if event.buttons()&Qt.RightButton:
			if self.mmode==0:
				self.motionRotate(self.rpressx - event.x(), self.rpressy - event.y())
				self.updateGL()
				self.rpressx = event.x()
				self.rpressy = event.y()
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
	
	def mouseReleaseEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
# 		if event.button()==Qt.RightButton:
# 			self.rmousedrag=None
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mouseup"), event)
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.rreleasex = event.x()
				self.rreleasey = event.y()
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return

	def motionRotate(self,x,y):
		if ( x == 0 and y == 0): return
		
		theta = atan2(-y,x)

		rotaxis_x = sin(theta)
		rotaxis_y = cos(theta)
		rotaxis_z = 0
		
		length = sqrt(x*x + y*y)
		angle = length/8.0*pi
		
		t3d = Transform3D()
		quaternion = {}
		quaternion["Omega"] = angle
		quaternion["n1"] = rotaxis_x
		quaternion["n2"] = rotaxis_y
		quaternion["n3"] = rotaxis_z
		
		t3d.set_rotation( EULER_SPIN, quaternion )
		
		size = len(self.t3d_stack)
		self.t3d_stack[size-1] = t3d*self.t3d_stack[size-1]
		
	def setScale(self,val):
		self.scale = val
		self.updateGL()
		
	def setAz(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["az"] = val - 180.0
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
		
	def setAlt(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["alt"] = val
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
		
	def setPhi(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["phi"] = val - 180.0
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
	
	def setThr(self,val):
		if (self.isothr != val):
			self.isothr = val
			self.createIsoDL()
			self.updateGL()

class EMImageInspector3D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist)
		
		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
		
		self.invtog = QtGui.QPushButton("Invert")
		self.invtog.setCheckable(1)
		self.vbl2.addWidget(self.invtog)
		
		self.ffttog = QtGui.QPushButton("FFT")
		self.ffttog.setCheckable(1)
		self.vbl2.addWidget(self.ffttog)

		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)
		
		self.mapp = QtGui.QPushButton("App")
		self.mapp.setCheckable(1)
		self.hbl2.addWidget(self.mapp)

		self.mmeas = QtGui.QPushButton("Meas")
		self.mmeas.setCheckable(1)
		self.hbl2.addWidget(self.mmeas)

		self.mmode=QtGui.QButtonGroup()
		self.mmode.setExclusive(1)
		self.mmode.addButton(self.mapp,0)
		self.mmode.addButton(self.mmeas,1)
		
		self.scale = ValSlider(self,(0.0,15.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0

		self.thr = ValSlider(self,(0.0,2.0),"Thr:")
		self.thr.setObjectName("az")
		self.thr.setValue(0.5)
		self.vbl.addWidget(self.thr)
			
		self.az = ValSlider(self,(0.0,360.0),"Az:")
		self.az.setObjectName("az")
		self.az.setValue(0.0)
		self.vbl.addWidget(self.az)
		
		self.alt = ValSlider(self,(0.01,180.0),"Alt:")
		self.alt.setObjectName("alt")
		self.alt.setValue(0.0)
		self.vbl.addWidget(self.alt)
		
		self.phi = ValSlider(self,(0.0,360.0),"Phi:")
		self.phi.setObjectName("phi")
		self.phi.setValue(0.0)
		self.vbl.addWidget(self.phi)

		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), target.setAz)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), target.setAlt)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), target.setPhi)
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), target.setThr)
		
	
	def newAz(self,val):
		if self.busy : return
		self.busy=1
		self.az.setValue(val, True)
		self.busy=0

	def newAlt(self,val):
		if self.busy : return
		self.busy=1
		self.alt.setValue(val, True)
		self.busy=0
	
	def newPhi(self,val):
		if self.busy : return
		self.busy=1
		self.phi.setValue(val, True)
		self.busy=0

	def setThrs(self,low,high,val):
		self.thr.setRange(low,high)
		self.thr.setValue(val, True)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage3D()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(33,33,33)
		e.process_inplace('testimage.x')
 		window.setData(e)

		# these lines are for testing shape rendering
# 		window.addShape("a",["rect",.2,.8,.2,20,20,80,80,2])
# 		window.addShape("b",["circle",.5,.8,.2,120,50,30.0,2])
# 		window.addShape("c",["line",.2,.8,.5,20,120,100,200,2])
# 		window.addShape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())

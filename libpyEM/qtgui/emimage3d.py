#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
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

from emimage3diso import EMIsosurface
from emimage3dvol import EMVolume

from time import *

MAG_INCREMENT_FACTOR = 1.1

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
		
		self.image = image
		self.currentselection = -1
		self.inspector = None
		#self.isosurface = EMIsosurface(image,self)
		#self.volume = EMVolume(image,self)
		self.viewables = []
		self.num_iso = 0
		self.num_vol = 0
		self.num_sli = 0
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		self.aspect=1.0
		self.fov = 50 # field of view angle used by gluPerspective
		
		self.addIsosurface()
		
	def timeout(self):
		self.updateGL()
		
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

		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		
		GL.glClearColor(0,0,0,0)
		
		# For the time being
		glEnable(GL_CULL_FACE);
		
		glPolygonMode(GL_FRONT,GL_FILL);

		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		for i in self.viewables:
			if (i.getType() == "volume"):
				glPushMatrix()
				i.render()
				glPopMatrix()
			
		for i in self.viewables:
			if (i.getType() != "volume"):
				glPushMatrix()
				i.render()
				glPopMatrix()


	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,0.001,1000000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		for i in self.viewables:
			i.resizeEvent()
			
	def setData(self,data):
		self.image = data
		for i in self.viewables:
			i.setData(data)
		#self.volume.setData(data)
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector3D(self)
		self.inspector.show()
	
	def closeEvent(self,event) :
		for i in self.viewables:
			i.closeEvent(event)
		if self.inspector: self.inspector.close()
		
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				self.inspector=EMImageInspector3D(self)
			self.showInspector(1)
		else:
			for i in self.viewables:
				i.mousePressEvent(event)

	def mouseMoveEvent(self, event):
		for i in self.viewables:
			i.mouseMoveEvent(event)
	
	def mouseReleaseEvent(self, event):
		for i in self.viewables:
			i.mouseReleaseEvent(event)
			
	def wheelEvent(self, event):
		for i in self.viewables:
			i.wheelEvent(event)

	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]

	def getSundryInspector(self):
		return self.viewables[self.currentselection].getInspector()
	
	def addIsosurface(self):
		self.viewables.append(EMIsosurface(self.image,self))
		self.num_iso += 1
		name = "Isosurface " + str(self.num_iso)
		self.viewables[len(self.viewables)-1].setName(name)
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def addVolume(self):
		self.viewables.append(EMVolume(self.image,self))
		self.num_vol += 1
		name = "Volume " + str(self.num_vol)
		self.viewables[len(self.viewables)-1].setName(name)
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def rowChanged(self,row):
		if ( row == self.currentselection ): return
		self.currentselection=row
		self.updateGL()
		
	def getCurrentName(self):
		return self.viewables[self.currentselection].getName()
	
	def getCurrentInspector(self):
		return self.viewables[self.currentselection].getInspector()
	
	def deleteCurrent(self):
		if ( len(self.viewables) == 0 ): return

		self.viewables.pop(self.currentselection)
		if (len(self.viewables) == 0 ) :  self.currentselection = -1
		elif ( self.currentselection == 0): pass
		else : self.currentselection -= 1

		self.updateGL()
	
class EMImageInspector3D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.listwidget = QtGui.QListWidget(self)
		self.listwidget.addItem("Isosurface")
		self.listwidget.setCurrentRow(0)
		self.vbl.addWidget(self.listwidget)
		
		self.hbl_buttons = QtGui.QHBoxLayout(self)
		self.hbl_buttons.setMargin(0)
		self.hbl_buttons.setSpacing(6)
		self.hbl_buttons.setObjectName("hbl_buttons")
		
		self.addIso = QtGui.QPushButton("Isosurface")
		self.hbl_buttons.addWidget(self.addIso)
		
		self.addVol = QtGui.QPushButton("Volume")
		self.hbl_buttons.addWidget(self.addVol)
		
		self.addSli = QtGui.QPushButton("Slices")
		self.hbl_buttons.addWidget(self.addSli)

		self.vbl.addLayout(self.hbl_buttons)
		
		self.hbl_buttons2 = QtGui.QHBoxLayout(self)
		self.delete = QtGui.QPushButton("Delete")
		self.hbl_buttons2.addWidget(self.delete)
		self.vbl.addLayout(self.hbl_buttons2)
		
		self.tabwidget = QtGui.QTabWidget(self)
		self.tabwidget.addTab(self.target.getSundryInspector(), "Isosurface")
		self.hbl.addLayout(self.vbl)
		self.hbl.addWidget(self.tabwidget)
		
		QtCore.QObject.connect(self.addIso, QtCore.SIGNAL("clicked()"), self.addIsosurface)
		QtCore.QObject.connect(self.addVol, QtCore.SIGNAL("clicked()"), self.addVolume)
		QtCore.QObject.connect(self.addSli, QtCore.SIGNAL("clicked()"), self.addSlices)
		QtCore.QObject.connect(self.delete, QtCore.SIGNAL("clicked()"), self.deleteSelection)
		
		QtCore.QObject.connect(self.listwidget, QtCore.SIGNAL("currentRowChanged(int)"), self.rowChanged)
		QtCore.QObject.connect(self.tabwidget, QtCore.SIGNAL("currentChanged(int)"), self.tabChanged)
		
	def tabChanged(self,tab):
		self.target.rowChanged(tab)
		self.listwidget.setCurrentRow(self.target.currentselection)
		
	def rowChanged(self,row):
		self.target.rowChanged(row)
		self.tabwidget.setCurrentIndex(self.target.currentselection)

	def addIsosurface(self):
		self.target.addIsosurface()
		self.updateSelection()
	
	def addVolume(self):
		self.target.addVolume()
		self.updateSelection()
	
	def updateSelection(self):
		self.listwidget.addItem(self.target.getCurrentName())
		self.listwidget.setCurrentRow(self.target.currentselection)
		self.tabwidget.addTab(self.target.getCurrentInspector(), self.target.getCurrentName())
		self.tabwidget.setCurrentIndex(self.target.currentselection)
	
	
	def addSlices(self):
		pass
	
	def deleteSelection(self):
		tmp = self.target.currentselection
		self.target.deleteCurrent()
		self.tabwidget.removeTab(tmp)
		self.listwidget.takeItem(tmp)
		
		
		#self.x_trans = QtGui.QDoubleSpinBox(self)
		#self.x_trans.setMinimum(-10000)
		#self.x_trans.setMaximum(10000)
		#self.x_trans.setValue(0.0)
		#self.hbl_trans.addWidget(self.x_trans)
		
		#self.y_label = QtGui.QLabel()
		#self.y_label.setText('y')
		#self.hbl_trans.addWidget(self.y_label)
		
		#self.y_trans = QtGui.QDoubleSpinBox(self)
		#self.y_trans.setMinimum(-10000)
		#self.y_trans.setMaximum(10000)
		#self.y_trans.setValue(0.0)
		#self.hbl_trans.addWidget(self.y_trans)
		
		
		#self.z_label = QtGui.QLabel()
		#self.z_label.setText('z')
		#self.hbl_trans.addWidget(self.z_label)
		
		#self.z_trans = QtGui.QDoubleSpinBox(self)
		#self.z_trans.setMinimum(-10000)
		#self.z_trans.setMaximum(10000)
		#self.z_trans.setValue(0.0)
		#self.hbl_trans.addWidget(self.z_trans)
		
		#self.hbl_src = QtGui.QHBoxLayout()
		#self.hbl_src.setMargin(0)
		#self.hbl_src.setSpacing(6)
		#self.hbl_src.setObjectName("hbl")
		#self.vbl.addLayout(self.hbl_src)
		
		#self.label_src = QtGui.QLabel()
		#self.label_src.setText('Rotation Convention')
		#self.hbl_src.addWidget(self.label_src)
		
		#self.src = QtGui.QComboBox(self)
		#self.load_src_options(self.src)
		#self.hbl_src.addWidget(self.src)
		
		## set default value -1 ensures that the val slider is updated the first time it is created
		#self.az = ValSlider(self,(-360.0,360.0),"az",-1)
		#self.az.setObjectName("az")
		#self.vbl.addWidget(self.az)
		
		#self.alt = ValSlider(self,(-180.0,180.0),"alt",-1)
		#self.alt.setObjectName("alt")
		#self.vbl.addWidget(self.alt)
		
		#self.phi = ValSlider(self,(-360.0,360.0),"phi",-1)
		#self.phi.setObjectName("phi")
		#self.vbl.addWidget(self.phi)
		
		#self.n3_showing = False
		
		#self.current_src = EULER_EMAN
		
		#QtCore.QObject.connect(self.zoom, QtCore.SIGNAL("valueChanged"), target.setZoom)
		#QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		#QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		#QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		#QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		#QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		#QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		#QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)

	#def setXYTrans(self, x, y):
		#self.x_trans.setValue(x)
		#self.y_trans.setValue(y)
	
	#def setTranslateScale(self, xscale,yscale,zscale):
		#self.x_trans.setSingleStep(xscale)
		#self.y_trans.setSingleStep(yscale)
		#self.z_trans.setSingleStep(zscale)

	#def updateRotations(self,t3d):
		#rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		#convention = self.src.currentText()
		#if ( self.src_map[str(convention)] == EULER_SPIN ):
			#self.n3.setValue(rot[self.n3.getLabel()],True)
		
		#self.az.setValue(rot[self.az.getLabel()],True)
		#self.alt.setValue(rot[self.alt.getLabel()],True)
		#self.phi.setValue(rot[self.phi.getLabel()],True)
	
	#def sliderRotate(self):
		#self.target.loadRotation(self.getCurrentRotation())
	
	#def getCurrentRotation(self):
		#convention = self.src.currentText()
		#rot = {}
		#if ( self.current_src == EULER_SPIN ):
			#rot[self.az.getLabel()] = self.az.getValue()
			
			#n1 = self.alt.getValue()
			#n2 = self.phi.getValue()
			#n3 = self.n3.getValue()
			
			#norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			#n1 /= norm
			#n2 /= norm
			#n3 /= norm
			
			#rot[self.alt.getLabel()] = n1
			#rot[self.phi.getLabel()] = n2
			#rot[self.n3.getLabel()] = n3
			
		#else:
			#rot[self.az.getLabel()] = self.az.getValue()
			#rot[self.alt.getLabel()] = self.alt.getValue()
			#rot[self.phi.getLabel()] = self.phi.getValue()
		
		#return Transform3D(self.current_src, rot)
	
	#def set_src(self, val):
		#t3d = self.getCurrentRotation()
		
		#if (self.n3_showing) :
			#self.vbl.removeWidget(self.n3)
			#self.n3.deleteLater()
			#self.n3_showing = False
			#self.az.setRange(-360,360)
			#self.alt.setRange(-180,180)
			#self.phi.setRange(-360,660)
		
		#if ( self.src_map[str(val)] == EULER_SPIDER ):
			#self.az.setLabel('phi')
			#self.alt.setLabel('theta')
			#self.phi.setLabel('psi')
		#elif ( self.src_map[str(val)] == EULER_EMAN ):
			#self.az.setLabel('az')
			#self.alt.setLabel('alt')
			#self.phi.setLabel('phi')
		#elif ( self.src_map[str(val)] == EULER_IMAGIC ):
			#self.az.setLabel('alpha')
			#self.alt.setLabel('beta')
			#self.phi.setLabel('gamma')
		#elif ( self.src_map[str(val)] == EULER_XYZ ):
			#self.az.setLabel('xtilt')
			#self.alt.setLabel('ytilt')
			#self.phi.setLabel('ztilt')
		#elif ( self.src_map[str(val)] == EULER_MRC ):
			#self.az.setLabel('phi')
			#self.alt.setLabel('theta')
			#self.phi.setLabel('omega')
		#elif ( self.src_map[str(val)] == EULER_SPIN ):
			#self.az.setLabel('Omega')
			#self.alt.setRange(-1,1)
			#self.phi.setRange(-1,1)
			
			#self.alt.setLabel('n1')
			#self.phi.setLabel('n2')
			
			#self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			#self.n3.setRange(-1,1)
			#self.n3.setObjectName("n3")
			#self.vbl.addWidget(self.n3)
			#QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
			#self.n3_showing = True
		
		#self.current_src = self.src_map[str(val)]
		#self.updateRotations(t3d)
	
	#def load_src_options(self,widgit):
		#self.load_src()
		#for i in self.src_strings:
			#widgit.addItem(i)
	
	## read src as 'supported rotation conventions'
	#def load_src(self):
		## supported_rot_conventions
		#src_flags = []
		#src_flags.append(EULER_EMAN)
		#src_flags.append(EULER_SPIDER)
		#src_flags.append(EULER_IMAGIC)
		#src_flags.append(EULER_MRC)
		#src_flags.append(EULER_SPIN)
		#src_flags.append(EULER_XYZ)
		
		#self.src_strings = []
		#self.src_map = {}
		#for i in src_flags:
			#self.src_strings.append(str(i))
			#self.src_map[str(i)] = i

	#def setZoom(self,newscale):
		#self.zoom.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage3D()
 	if len(sys.argv)==1 : 
		pass
		e = EMData()
		e.set_size(64,64,64)
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

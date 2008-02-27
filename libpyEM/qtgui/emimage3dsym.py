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
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from time import *

from emimage3dobject import EMImage3DObject
from emimage3dobject import Camera

MAG_INCREMENT_FACTOR = 1.1

class EM3DSymViewer(EMImage3DObject):
	def __init__(self, parent=None):
		EMImage3DObject.__init__(self)
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMSymInspector(self)
	
	def getType(self):
		return "Symmetry Viewer"

	def init(self):
		self.data=None

		self.mmode=0
		self.cam = Camera()
		
		self.cube = False
		self.inspector=None
		
		self.sym_dl = 0
		self.spheredl = 0
		self.highresspheredl = 0

		self.glcontrast = 1.0
		self.glbrightness = 0.0
		
		self.rank = 1
		
		self.force_update = False
		self.force_force_update = False
		
		self.sym = None
		self.prop = None
		self.perturb = False
		self.nomirror = True
		
		self.radius = 50
		
		
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)
		
		
		#self.glbasicobjects = EMBasicOpenGLObjects()
		#print self.glbasicobjects.getSphereDL()
		
	def setRadius(self,radius):
		if ( radius > 0 ):
			self.radius = radius
			self.force_force_update = True
			self.parent.updateGL()
		else:
			print "Error, tried to set a zero or negative radius (",radius,")"
			exit(1)
		
	def genCurrentDisplayList(self):
		sym = self.inspector.getSym()
		prop = self.inspector.getProp()
		mirror = not self.inspector.getMirror()
		perturb = self.inspector.getPerturb()
		
		
		if self.sym == sym and self.prop == prop and self.nomirror == mirror and self.perturb == perturb and self.force_force_update == False: return
		else:
			self.sym = sym
			self.prop = prop
			self.nomirror = mirror
			self.perturb = perturb
		
		if self.spheredl == 0:
			self.spheredl=glGenLists(1)
				
			glNewList(self.spheredl,GL_COMPILE)
			gluSphere(self.gq,.5,4,2)
			glEndList()

		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
				
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()

		if ( self.sym_dl != 0 ): glDeleteLists( self.sym_dl, 1)
		
		self.sym_dl = glGenLists(1)
	
		if (self.sym_dl == 0):
			self.sym = None
			self.prop = None
			return #OpenGL is not initialized yet
		sym_object = parsesym(str(sym))
		if mirror == True : val = 0
		else: val = 1
		og = "eman" + ":prop=" + str(prop) + ":inc_mirror=" + str(val)
		if ( perturb == True ) : val = 1
		else: val = 0
		og += ":perturb="+ str(val)
		[og_name,og_args] = parsemodopt(og)
		eulers = sym_object.gen_orientations(og_name, og_args)

		glNewList(self.sym_dl,GL_COMPILE)
		for i in eulers:
			d = i.get_rotation()
			glPushMatrix()
			glRotate(d["az"],0,0,1)
			glRotate(d["alt"],1,0,0)
			glRotate(d["phi"],0,0,1)
			#print d["phi"],d["alt"],d["az"]
			glTranslate(0,0,self.radius)
			glCallList(self.spheredl)
			glPopMatrix()
		glEndList()
		
	def render(self):
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		glEnable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		self.cam.position()
		
		if ( self.sym_dl == 0 or self.force_update):
			self.genCurrentDisplayList()
			self.force_update = False
			if ( self.sym_dl == 0 ) : 
				print "error, you can't draw an empty list"
				return
		
		
		glColor(.9,.2,.8)
		# this is a nice light blue color (when lighting is on)
		# and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.2,.2,.8,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,50.0)
		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glPushMatrix()
		glCallList(self.sym_dl)
		glPopMatrix()
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		# FIXME the approach here is very inefficient
		glLoadIdentity()
		glTranslate(0,0,-2)
		glScale(2*self.radius,2*self.radius,1)
		glTranslate(-0.5,-0.5,0)
		self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)

	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMSymInspector(self)
		self.inspector.updateRotations(t3d)
	
	def getInspector(self):
		if not self.inspector : self.inspector=EMSymInspector(self)
		return self.inspector
		
	def regenDL(self, dummy=False):
		self.force_update = True
		self.parent.updateGL()

class EMSymViewerWidget(QtOpenGL.QGLWidget):
	
	allim=WeakKeyDictionary()
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMSymViewerWidget.allim[self]=0
		
		self.fov = 50 # field of view angle used by gluPerspective
		
		self.sliceviewer = EM3DSymViewer(self)
		self.sliceviewer.cam.setCamTrans("default_z",-100)

	def initializeGL(self):
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		
		GL.glClearColor(0,0,0,0)
		#GL.glClearAccum(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
	def paintGL(self):
		#glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		#self.sliceviewer.cam.setCamTrans("default_z",-100)
		#glTranslate(0,0,-100)
		
		glPushMatrix()
		self.sliceviewer.render()
		glPopMatrix()
		
		#glAccum(GL_ADD, self.sliceviewer.glbrightness)
		#glAccum(GL_ACCUM, self.sliceviewer.glcontrast)
		#glAccum(GL_RETURN, 1.0)
		
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,1,5000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.sliceviewer.resizeEvent()

	def showInspector(self,force=0):
		self.sliceviewer.showInspector(self,force)

	def closeEvent(self,event) :
		self.sliceviewer.closeEvent(event)
		
	def mousePressEvent(self, event):
		self.sliceviewer.mousePressEvent(event)
		self.emit(QtCore.SIGNAL("mousedown"), event)
		
	def mouseMoveEvent(self, event):
		self.sliceviewer.mouseMoveEvent(event)
		self.emit(QtCore.SIGNAL("mousedrag"), event)
	
	def mouseReleaseEvent(self, event):
		self.sliceviewer.mouseReleaseEvent(event)
		self.emit(QtCore.SIGNAL("mouseup"), event)
			
	def wheelEvent(self, event):
		self.sliceviewer.wheelEvent(event)

	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]
		

class EMSymInspector(QtGui.QWidget):
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
		
		
		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
	
		self.cubetog = QtGui.QPushButton("Cube")
		self.cubetog.setCheckable(1)
		self.vbl2.addWidget(self.cubetog)
		
		self.defaults = QtGui.QPushButton("Defaults")
		self.vbl2.addWidget(self.defaults)
		
		self.vbl.addWidget(self.getMainTab())
		
		self.n3_showing = False
		
		self.current_src = EULER_EMAN
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
		QtCore.QObject.connect(self.sym_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.symChanged)
		QtCore.QObject.connect(self.sym_text, QtCore.SIGNAL("editingFinished()"), target.regenDL)
		QtCore.QObject.connect(self.prop_text, QtCore.SIGNAL("editingFinished()"), target.regenDL)
		QtCore.QObject.connect(self.mirror_checkbox, QtCore.SIGNAL("stateChanged(int)"), target.regenDL)
		QtCore.QObject.connect(self.perturbtog, QtCore.SIGNAL("toggled(bool)"), self.perturbtoggled)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		#QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggleCube)
		QtCore.QObject.connect(self.defaults, QtCore.SIGNAL("clicked(bool)"), self.setDefaults)
		#QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		
	def perturbtoggled(self,val):
		self.target.regenDL()
	
	def setDefaults(self):
		self.x_trans.setValue(0.0)
		self.y_trans.setValue(0.0)
		self.z_trans.setValue(0.0)
		self.scale.setValue(1.0)
		self.glcontrast.setValue(1.0)
		self.glbrightness.setValue(0.0)
		
		self.az.setValue(0.0)
		self.alt.setValue(0.0)
		self.phi.setValue(0.0)

	def symChanged(self, sym):
		if sym == ' D ' or sym == ' C ' or sym == ' H ':
			self.sym_text.setEnabled(True)
		else:
			self.sym_text.setEnabled(False)
		
		self.target.regenDL()

	def getSym(self):
		sym = self.sym_map[str(self.sym_combo.currentText())]
		if sym in ['c','d','h']:
			sym = sym+self.sym_text.displayText()
		return sym
		
	def getProp(self):
		return float(self.prop_text.displayText())
	
	def getMirror(self):
		return self.mirror_checkbox.checkState() == Qt.Checked
	
	def getPerturb(self):
		return self.perturbtog.isChecked()

	def getMainTab(self):
	
		self.maintab = QtGui.QWidget()
		maintab = self.maintab
		maintab.vbl = QtGui.QVBoxLayout(self.maintab)
		maintab.vbl.setMargin(0)
		maintab.vbl.setSpacing(6)
		maintab.vbl.setObjectName("Main")
		
		
		self.scale = ValSlider(self,(0.01,30.0),"Zoom:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		maintab.vbl.addWidget(self.scale)
		
		self.hbl_sym = QtGui.QHBoxLayout()
		self.hbl_sym.setMargin(0)
		self.hbl_sym.setSpacing(6)
		self.hbl_sym.setObjectName("Sym")
		maintab.vbl.addLayout(self.hbl_sym)
		
		self.sym_combo = QtGui.QComboBox(maintab)
		self.symmetries = []
		self.symmetries.append(' Icosahedral ')
		self.symmetries.append(' Octahedral ')
		self.symmetries.append(' Tetrahedral ')
		self.symmetries.append(' D ')
		self.symmetries.append(' C ')
		self.symmetries.append(' H ')
		self.sym_map = {}
		self.sym_map[" Icosahedral "] = "icos"
		self.sym_map[" Octahedral "] = "oct"
		self.sym_map[" Tetrahedral "] = "tet"
		self.sym_map[" D "] = "d"
		self.sym_map[" C "] = "c"
		self.sym_map[" H "] = "h"
		for i in self.symmetries: self.sym_combo.addItem(i)
		self.hbl_sym.addWidget(self.sym_combo)
		
		self.sym_label = QtGui.QLabel()
		self.sym_label.setText('C/D/H sym')
		self.hbl_sym.addWidget(self.sym_label)
		
		
		self.hbl_sym2 = QtGui.QHBoxLayout()
		self.hbl_sym2.setMargin(0)
		self.hbl_sym2.setSpacing(6)
		self.hbl_sym2.setObjectName("Sym2")
		maintab.vbl.addLayout(self.hbl_sym2)
		
		self.mirror_checkbox = QtGui.QCheckBox("Mirror")
		self.hbl_sym2.addWidget(self.mirror_checkbox)
		
		self.perturbtog = QtGui.QPushButton("Perturb")
		self.perturbtog.setCheckable(1)
		self.hbl_sym2.addWidget(self.perturbtog)
		
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.sym_text = QtGui.QLineEdit(self)
		self.sym_text.setValidator(self.pos_int_validator)
		self.sym_text.setText("5")
		self.sym_text.setFixedWidth(50)
		self.hbl_sym.addWidget(self.sym_text)
		self.sym_text.setEnabled(False)
		
		self.prop_label = QtGui.QLabel()
		self.prop_label.setText('prop')
		self.hbl_sym.addWidget(self.prop_label)
		
		self.pos_double_validator = QtGui.QDoubleValidator(self)
		self.pos_double_validator.setBottom(0.05)
		self.prop_text = QtGui.QLineEdit(self)
		self.prop_text.setValidator(self.pos_double_validator)
		self.prop_text.setText("2.0")
		self.prop_text.setFixedWidth(50)
		self.hbl_sym.addWidget(self.prop_text)
		
		self.glcontrast = ValSlider(maintab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		maintab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(maintab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		maintab.vbl.addWidget(self.glbrightness)
		
		self.cbb = QtGui.QComboBox(maintab)
		self.vbl.addWidget(self.cbb)
		self.cbb.deleteLater()

		self.hbl_trans = QtGui.QHBoxLayout()
		self.hbl_trans.setMargin(0)
		self.hbl_trans.setSpacing(6)
		self.hbl_trans.setObjectName("Trans")
		maintab.vbl.addLayout(self.hbl_trans)
		
		self.x_label = QtGui.QLabel()
		self.x_label.setText('x')
		self.hbl_trans.addWidget(self.x_label)
		
		self.x_trans = QtGui.QDoubleSpinBox(maintab)
		self.x_trans.setMinimum(-10000)
		self.x_trans.setMaximum(10000)
		self.x_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.x_trans)
		
		self.y_label = QtGui.QLabel()
		self.y_label.setText('y')
		self.hbl_trans.addWidget(self.y_label)
		
		self.y_trans = QtGui.QDoubleSpinBox(maintab)
		self.y_trans.setMinimum(-10000)
		self.y_trans.setMaximum(10000)
		self.y_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.y_trans)
		
		self.z_label = QtGui.QLabel()
		self.z_label.setText('z')
		self.hbl_trans.addWidget(self.z_label)
		
		self.z_trans = QtGui.QDoubleSpinBox(maintab)
		self.z_trans.setMinimum(-10000)
		self.z_trans.setMaximum(10000)
		self.z_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.z_trans)
		
		self.hbl_src = QtGui.QHBoxLayout()
		self.hbl_src.setMargin(0)
		self.hbl_src.setSpacing(6)
		self.hbl_src.setObjectName("hbl")
		maintab.vbl.addLayout(self.hbl_src)
		
		self.label_src = QtGui.QLabel()
		self.label_src.setText('Rotation Convention')
		self.hbl_src.addWidget(self.label_src)
		
		self.src = QtGui.QComboBox(maintab)
		self.load_src_options(self.src)
		self.hbl_src.addWidget(self.src)
		
		# set default value -1 ensures that the val slider is updated the first time it is created
		self.az = ValSlider(maintab,(-360.0,360.0),"az",-1)
		self.az.setObjectName("az")
		maintab.vbl.addWidget(self.az)
		
		self.alt = ValSlider(maintab,(-180.0,180.0),"alt",-1)
		self.alt.setObjectName("alt")
		maintab.vbl.addWidget(self.alt)
		
		self.phi = ValSlider(maintab,(-360.0,360.0),"phi",-1)
		self.phi.setObjectName("phi")
		maintab.vbl.addWidget(self.phi)
		
		return maintab
	
	def setDefaults(self):
		self.x_trans.setValue(0.0)
		self.y_trans.setValue(0.0)
		self.z_trans.setValue(0.0)
		self.scale.setValue(1.0)
		self.glcontrast.setValue(1.0)
		self.glbrightness.setValue(0.0)
		
		self.az.setValue(0.0)
		self.alt.setValue(0.0)
		self.phi.setValue(0.0)

	def setXYTrans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
	
	def setTranslateScale(self, xscale,yscale,zscale):
		self.x_trans.setSingleStep(xscale)
		self.y_trans.setSingleStep(yscale)
		self.z_trans.setSingleStep(zscale)

	def updateRotations(self,t3d):
		rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		convention = self.src.currentText()
		if ( self.src_map[str(convention)] == EULER_SPIN ):
			self.n3.setValue(rot[self.n3.getLabel()],True)
		
		self.az.setValue(rot[self.az.getLabel()],True)
		self.alt.setValue(rot[self.alt.getLabel()],True)
		self.phi.setValue(rot[self.phi.getLabel()],True)
	
	def sliderRotate(self):
		self.target.loadRotation(self.getCurrentRotation())
	
	def getCurrentRotation(self):
		convention = self.src.currentText()
		rot = {}
		if ( self.current_src == EULER_SPIN ):
			rot[self.az.getLabel()] = self.az.getValue()
			
			n1 = self.alt.getValue()
			n2 = self.phi.getValue()
			n3 = self.n3.getValue()
			
			norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			n1 /= norm
			n2 /= norm
			n3 /= norm
			
			rot[self.alt.getLabel()] = n1
			rot[self.phi.getLabel()] = n2
			rot[self.n3.getLabel()] = n3
			
		else:
			rot[self.az.getLabel()] = self.az.getValue()
			rot[self.alt.getLabel()] = self.alt.getValue()
			rot[self.phi.getLabel()] = self.phi.getValue()
		
		return Transform3D(self.current_src, rot)
	
	def set_src(self, val):
		t3d = self.getCurrentRotation()
		
		if (self.n3_showing) :
			self.maintab.vbl.removeWidget(self.n3)
			self.n3.deleteLater()
			self.n3_showing = False
			self.az.setRange(-360,360)
			self.alt.setRange(-180,180)
			self.phi.setRange(-360,660)
		
		if ( self.src_map[str(val)] == EULER_SPIDER ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('psi')
		elif ( self.src_map[str(val)] == EULER_EMAN ):
			self.az.setLabel('az')
			self.alt.setLabel('alt')
			self.phi.setLabel('phi')
		elif ( self.src_map[str(val)] == EULER_IMAGIC ):
			self.az.setLabel('alpha')
			self.alt.setLabel('beta')
			self.phi.setLabel('gamma')
		elif ( self.src_map[str(val)] == EULER_XYZ ):
			self.az.setLabel('xtilt')
			self.alt.setLabel('ytilt')
			self.phi.setLabel('ztilt')
		elif ( self.src_map[str(val)] == EULER_MRC ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('omega')
		elif ( self.src_map[str(val)] == EULER_SPIN ):
			self.az.setLabel('Omega')
			self.alt.setRange(-1,1)
			self.phi.setRange(-1,1)
			
			self.alt.setLabel('n1')
			self.phi.setLabel('n2')
			
			self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			self.n3.setRange(-1,1)
			self.n3.setObjectName("n3")
			self.maintab.vbl.addWidget(self.n3)
			QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
			self.n3_showing = True
		
		self.current_src = self.src_map[str(val)]
		self.updateRotations(t3d)
	
	def load_src_options(self,widgit):
		self.load_src()
		for i in self.src_strings:
			widgit.addItem(i)
	
	# read src as 'supported rotation conventions'
	def load_src(self):
		# supported_rot_conventions
		src_flags = []
		src_flags.append(EULER_EMAN)
		src_flags.append(EULER_SPIDER)
		src_flags.append(EULER_IMAGIC)
		src_flags.append(EULER_MRC)
		src_flags.append(EULER_SPIN)
		src_flags.append(EULER_XYZ)
		
		self.src_strings = []
		self.src_map = {}
		for i in src_flags:
			self.src_strings.append(str(i))
			self.src_map[str(i)] = i
		
	
	def setColors(self,colors,current_color):
		for i in colors:
			self.cbb.addItem(i)

	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setScale(self,newscale):
		self.scale.setValue(newscale)
	
	
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMSymViewerWidget()
	window2=EMParentWin(window)
	window2.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())

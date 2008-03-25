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

# EMFloatingWidgets.py  Steve Ludtke  08/06/2006
# An experimental version of emimage.py that displays an image in a 3D context
# using texture mapping. Not a fully fleshed out class at this point

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider

from EMAN2 import *
from emimageutil import *
from emimage2d import *
from emimage3d import *
from emimagemx import *
from math import sqrt

from emimage import EMImage

from emglobjects import EMViewportDepthTools, Camera2, EMBasicOpenGLObjects, Camera
from emimage2dtex import *


height_plane = 500

class EMGLView3D:
	"""
	A view of an EMAN2 3D type, such as an isosurface or a 
	volume rendition, etc.
	"""
	def __init__(self, parent,image=None):
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.motiondull = 3.0
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		self.w = image.get_xsize()	# width of window
		self.h = image.get_ysize()	# height of window
		self.d = image.get_zsize()	# depth of the window
		self.sizescale = 1.0		# scale/zoom factor
		self.changefactor = 1.1		# used to zoom
		self.invchangefactor = 1.0/self.changefactor # used to invert zoom
		
		self.drawable = EMImage3DCore(image,self)		# the object that is drawable (has a draw function)
		self.drawable.cam.basicmapping = True
		self.drawable.cam.motiondull = 3.0
		#self.drawable.supressInspector = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.vdtools
		
		self.psets = []
		self.modelmatrices = []
		
	def width(self):
		try:
			return int(self.w)
		except:
			return 0
	
	def height(self):
		try:
			return int(self.h)
		except:
			return 0
	
	def depth(self):
		try:
			return int(self.d)
		except:
			return 0
	
	def setData(self,data):
		try: self.drawable.setData(data)
		except: pass
		
	def paintGL(self):
		self.psets = []
		self.planes = []
		self.modelmatrices = []
		self.cam.position()
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		
		glPushMatrix()
		glTranslatef(-self.width()/2.0,0,0)
		glRotatef(-90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.planes.append(('zy'))
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(self.width()/2.0,0,0)
		glRotatef(90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('yz'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
				
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,self.height()/2.0,0)
		glRotatef(-90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('xz'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,-self.height()/2.0,0)
		glRotatef(90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('zx'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,-self.depth()/2.0)
		glRotatef(180,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True):
				self.planes.append(('yx'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,self.depth()/2.0)
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('xy'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		self.drawable.render()
		glPopMatrix()
		
		if not lighting: glDisable(GL_LIGHTING)
		
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):	
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.drawable.inspector)
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)

	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)
	
	def update(self):
		pass
		self.parent.updateGL()
	
	def updateGL(self):
		self.parent.updateGL()
	
	def isinwin(self,x,y):
		val = False
		for i,p in enumerate(self.psets):
			if self.vdtools.isinwinpoints(x,y,p):
				val = True
				self.drawable.cam.plane = self.planes[i]
				self.vdtools.setModelMatrix(self.modelmatrices[i])
				break
		return val
	
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def leaveEvent(self):
		pass
	
	def toolTipEvent(self,event):
		pass
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
		
	def get_render_dims_at_depth(self, depth):
		return self.parent.get_render_dims_at_depth(depth)
	
class EMGLView2D:
	"""
	A view of a 2D drawable type, such as a single 2D image or a matrix of 2D images
	
	"""
	def __init__(self, parent,image=None):
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		
		if isinstance(image,list):
			if len(image) == 1:
				self.become2DImage(image[0])
			else:
				self.drawable = EMImageMXCore(image,self)
				self.w = image[0].get_xsize()
				self.h = image[0].get_ysize()
		elif isinstance(image,EMData):
			self.become2DImage(image)
		
		self.drawable.supressInspector = True
		self.initflag = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.sizescale = 1.0
		self.changefactor = 1.1
		self.invchangefactor = 1.0/self.changefactor
		
	def become2DImage(self,a):
		self.drawable = EMImage2DCore(a,self)
		#self.drawable.originshift = False
		self.w = a.get_xsize()
		self.h = a.get_ysize()
		
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def width(self):
		try:
			return int(self.sizescale*self.w)
		except:
			return 0
	
	def height(self):
		try:
			return int(self.sizescale*self.h)
		except:
			return 0
	
	def setData(self,data):
		self.drawable.setData(data)
		
	def initializeGL(self):
		self.drawable.initializeGL()
	
	def viewportHeight(self):
		return self.parent.height()	
	
	def viewportWidth(self):
		return self.parent.width()
	
	def testBoundaries(self):
		'''
		Called when the image is first drawn, this resets the dimensions of this object
		if it is larger than the current size of the viewport. It's somewhat of a hack,
		but it's early stages in the design
		'''
		h = self.vdtools.getMappedHeight()
		w = self.vdtools.getMappedWidth()
		
		if ( w > self.viewportWidth() ):
			self.w = self.viewportWidth()/self.sizescale
		if ( h > self.viewportHeight() ):
			self.h = self.viewportHeight()/self.sizescale
	
	def paintGL(self):
		self.cam.position()
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		if (self.initflag == True):
			self.testBoundaries()
			self.initflag = False

		
		#self.mediator.checkBoundaryIssues()
		if (self.updateFlag):
			self.drawable.resizeEvent(self.width(),self.height())
			self.updateFlag = False
		glPushMatrix()
		glTranslatef(-self.width()/2.0,-self.height()/2.0,0)
		try: self.drawable.render()
		except Exception, inst:
			print type(inst)     # the exception instance
			print inst.args      # arguments stored in .args
			print int
		glPopMatrix()
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		if self.drawFrame: self.vdtools.drawFrame()
		if not lighting: glDisable(GL_LIGHTING)
		
	def update(self):
		self.parent.updateGL()
	
	def updateGL(self):
		self.parent.updateGL()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):	
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.drawable.inspector)
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)
			#self.drawable.mousePressEvent(event)
		
		#self.updateGL()
	
	def scaleEvent(self,delta):
		if ( delta > 0 ):
			self.sizescale *= self.changefactor
		elif ( delta < 0 ):
			self.sizescale *= self.invchangefactor

		self.drawable.resizeEvent(self.width(),self.height())
	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.scaleEvent(event.delta())
			#print "updating",self.drawWidth(),self.drawHeight()
		else:
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)

		#self.updateGL()
	def emit(self, signal, event):
		self.parent.emit(signal,event)
	
	def leaveEvent(self):
		self.drawable.leaveEvent()
	
	def toolTipEvent(self,event):
		pass
	
	def isinwin(self,x,y):
		return self.vdtools.isinwin(x,y)

class EMGLViewQtWidget:
	def __init__(self, parent=None, qwidget=None, widget_parent=None):
		self.parent = parent
		self.qwidget = qwidget
		self.drawFrame = True
		self.mapcoords = True
		self.itex = 0
		self.genTexture = True
		self.click_debug = False
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		self.cam.motionRotate(0,0)
		self.borderwidth = 3.0
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.setQtWidget(qwidget)
		self.P_inv = None
		self.childreceiver = None
		self.widget_parent = widget_parent
		
		self.current = None
		self.previous = None
		
		self.e2children = []
		self.is_child = False
		
		self.vdtools = EMViewportDepthTools(self)
	
	def __del__(self):
		if (self.itex != 0 ):
			self.parent.deleteTexture(self.itex)
		
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def width(self):
		return self.qwidget.width()
	
	def height(self):
		return self.qwidget.height()
	
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def setQtWidget(self, widget, delete_current = False):
		if ( delete_current and self.qwidget != None ):
			self.qwidget.deleteLater()
		
		self.qwidget = widget
		
		if ( widget != None ):
			#self.qwidget.setVisible(True)
			self.qwidget.setEnabled(True)
			self.genTexture = True
			self.updateTexture()
			
	def updateTexture(self):
		if ( self.itex == 0 or self.genTexture == True ) : 
			if (self.itex != 0 ):
				#passpyth
				self.parent.deleteTexture(self.itex)
			self.genTexture = False
			##print "binding texture"
			#self.qwidget.setVisible(True)
			#self.qwidget.repaint()
			pixmap = QtGui.QPixmap.grabWidget(self.qwidget)
			#self.qwidget.setVisible(False)
			if (pixmap.isNull() == True ): print 'error, the pixmap was null'
			self.itex = self.parent.bindTexture(pixmap)
			if ( self.itex == 0 ): print 'Error - I could not generate the texture'
		
	def paintGL(self):
		#print "paintGL children"
		if (self.qwidget == None or self.itex == 0) :
			#print "no widget - paintGL children return" 
			return
		
		self.cam.position()
		
		# make sure the vdtools store the current matrices
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		
		glPushMatrix()
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
		glTexCoord2f(1.,0.)
		glVertex( self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
		glTexCoord2f(1.,1.)
		glVertex( self.qwidget.width()/2.0, self.qwidget.height()/2.0)
		glTexCoord2f(0.,1.)
		glVertex( -self.qwidget.width()/2.0,self.qwidget.height()/2.0)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		glPopMatrix()
	
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		if self.drawFrame:
			try: self.vdtools.drawFrame()
			except Exception, inst:
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args
				print int
		if (not lighting): glDisable(GL_LIGHTING)
		
		# now draw children if necessary - such as a qcombobox list view that has poppud up
		for i in self.e2children:
			glPushMatrix()
			try:
				i.paintGL()
			except Exception, inst:
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args
				print int
			glPopMatrix()
	
	def isinwin(self,x,y):
		for i in self.e2children:
			if i.isinwin(x,y):
				self.childreceiver = i
				return True
		
		return self.vdtools.isinwin(x,y)
	
	def eyeCoordsDif(self,x1,y1,x2,y2):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2)
			
	def mouseinwin(self,x,y,width,height):
		return self.vdtools.mouseinwin(x,y,width,height)

	def toolTipEvent(self,event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.toolTip(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: 
			QtGui.QToolTip.hideText()
			self.genTexture = True
			self.updateTexture()
			return
	
		p1 = QtCore.QPoint(event.x(),event.y())
		p2 = self.parent.mapToGlobal(p1)
		QtGui.QToolTip.showText(p2,cw.toolTip())
	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.wheelEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			qme=QtGui.QWheelEvent(lp,event.delta(),event.buttons(),event.modifiers(),event.orientation())
			QtCore.QCoreApplication.sendEvent(cw,qme)
			self.genTexture = True
			self.updateTexture()
	
	def mouseDoubleClickEvent(self, event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.mouseDoubleClickEvent(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: return
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		lp=cw.mapFromGlobal(gp)
		if (isinstance(cw,QtGui.QComboBox)):
			print "it's a combo"
		else:
			qme=QtGui.mouseDoubleClickEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			#self.qwidget.setVisible(True)
			QtCore.QCoreApplication.sendEvent(cw,qme)
			#self.qwidget.setVisible(False)
		self.genTexture = True
		self.updateTexture()
		
	def get_depth_for_height(self,height_plane):
		return self.parent.get_depth_for_height(height_plane)
	
	def mousePressEvent(self, event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mousePressEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
				
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			##print cw.objectName()
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				cw.showPopup()
				cw.hidePopup()
				widget = EMGLViewQtWidget(self.parent,None,cw);
				widget.setQtWidget(cw.view())
				widget.cam.loadIdentity()	
				widget.cam.setCamTrans("x",cw.geometry().x()-self.width()/2.0+cw.view().width()/2.0)
				widget.cam.setCamTrans("y",((self.height()/2.0-cw.geometry().y())-cw.view().height()/2.0))
				widget.cam.setCamTrans("z",0.1)
				widget.drawFrame = False
				self.e2children.append(widget)
				self.e2children[0].is_child = True
			else:
				qme=QtGui.QMouseEvent( event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child): QtCore.QCoreApplication.sendEvent(self.qwidget,qme)
				else: QtCore.QCoreApplication.sendEvent(cw,qme)
				
			self.genTexture = True
			self.updateTexture()
		
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mouseMoveEvent(event)
				self.childreceiver = None
				return
			else:
				l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
				cw=self.qwidget.childAt(l[0],l[1])
				self.current = cw
				if ( self.current != self.previous ):
					QtGui.QToolTip.hideText()
					if ( self.current != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Enter)
						QtCore.QCoreApplication.sendEvent(self.current,qme)
						
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
				
				self.previous = self.current
				if cw == None:
					QtGui.QToolTip.hideText()
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
						self.genTexture = True
						self.updateTexture()
					return
				gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				lp=cw.mapFromGlobal(gp)
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				QtCore.QCoreApplication.sendEvent(cw,qme)
			# FIXME
			# setting the genTexture flag true here causes the texture to be regenerated
			# when the mouse moves over it, which is inefficient.
			# The fix is to only set the genTexture flag when mouse movement
			# actually causes a change in the appearance of the widget (for instance, list boxes from comboboxes)
			self.genTexture = True
			self.updateTexture()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				#try:
				self.childreceiver.mouseReleaseEvent(event)
				self.childreceiver = None
				self.e2children.pop()
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				print "it's a combo"
				#cw.showPopup()
			else:
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child):
					##print self.qwidget
					##print self.qwidget.currentIndex().row()
					##print self.widget_parent
					##print self.qwidget.rect().left(),self.qwidget.rect().right(),self.qwidget.rect().top(),self.qwidget.rect().bottom()
					##print lp.x(),lp.y()
					self.widget_parent.setCurrentIndex(self.qwidget.currentIndex().row())
					#self.widget_parent.changeEvent(QtCore.QEvent())
					#self.widget_parent.highlighted(self.qwidget.currentIndex().row())
					#self.qwidget.commitData(self.qwidget.parent())
					##print self.qwidget.currentText()
					#self.widget_parent.setVisible(True)
					#self.widget_parent.setEnabled(True)
					#self.qwidget.setVisible(True)
					#QtCore.QCoreApplication.sendEvent(self.widget_parent,qme)
					#self.qwidget.setVisible(False)
					self.widget_parent.emit(QtCore.SIGNAL("activated(QString)"),self.widget_parent.itemText(self.qwidget.currentIndex().row()))
				else:
					#self.qwidget.setVisible(True)
					QtCore.QCoreApplication.sendEvent(cw,qme)
					#self.qwidget.setVisible(False)
			
			self.genTexture = True
			self.updateTexture()
		
	def leaveEvent(self):
		if (self.current != None) : 
			qme = QtCore.QEvent(QtCore.QEvent.Leave)
			QtCore.QCoreApplication.sendEvent(self.current,qme)
			self.current = None
			self.previouse = None
			self.genTexture = True
			self.updateTexture()
			
	def enterEvent():
		pass
	def timerEvent(self,event=None):
		pass
		#self.cam.motionRotate(.2,.2)


class EMFloatingWidgets(QtOpenGL.QGLWidget):
	def __init__(self,parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		# enable multisampling to combat aliasing
		fmt.setSampleBuffers(True)
		# stenciling is for object dependent shading
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.setMouseTracking(True)
		self.fov = 2*180*atan2(1,5)/pi
		
		self.floatwidget = EMFloatingWidgetsCore(self)
		
		self.cam = Camera()
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]

	
	def initializeGL(self):
		#print "initializeGL"
		glClearColor(0,0,0,0)
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
	
		# enable multisampling to combat aliasing
		if ( "GL_ARB_multisample" in glGetString(GL_EXTENSIONS) ): glEnable(GL_MULTISAMPLE)
		else: glDisable(GL_MULTISAMPLE)
		
	def paintGL(self):
		#print "paintGL"
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
			
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
		self.floatwidget.render()
		
	def resizeGL(self, width, height):
		#print "resizeGL"
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		#glFrustum(-1.*width/height,1.*width/height, -1.,1., 5.,15.)
		
		# fov angle is the given by
		#self.fov = 2*180*atan2(1,5)/pi
		# aspect ratio is given by
		self.aspect = float(self.width())/float(self.height())
		# this is the same as the glFrustum call above
		depth = self.get_depth_for_height(height_plane)
		gluPerspective(self.fov,self.aspect,depth-depth/4,depth+depth/4)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		try: self.floatwidget.resizeEvent(width,height)
		except: print "couldn't resize floatwidget"
		
	def mousePressEvent(self, event):
		self.floatwidget.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		self.floatwidget.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.floatwidget.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		self.floatwidget.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		self.floatwidget.wheelEvent(event)

	def toolTipEvent(self, event):
		self.floatwidget.wheelEvent(event)
		QtGui.QToolTip.hideText()
		

	def dragMoveEvent(self,event):
		print "received drag move event, but I don't do anything about it :("
		
	def event(self,event):
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		self.floatwidget.hoverEvent(event)

class EMFloatingWidgetsCore:
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		#print "init"
		self.parent = parent
	
		self.imtex=0
		self.current = None
		self.previous = None
	
		self.initFlag = True
		self.qwidgets = []
		
		#print "init done"
	
		
	def get_depth_for_height(self, height):
		try: 
			return self.parent.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			return 0

	def height(self):
		return self.parent.height()
	
	def width(self):
		return self.parent.width()

	def updateGL(self):
		pass
		#self.parent.updateGL()

	def addQtWidgetDrawer(self,widget):
		w = EMGLViewQtWidget(self)
		w.setQtWidget(widget)
		self.qwidgets.append(w)
		
		#print "initializeGL done"
	def render(self):
		
		if ( self.initFlag == True ):
			self.fd = QtGui.QFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("finished(int)"), self.finished)
			self.fd.show()
			self.fd.hide()
			self.qwidgets.append(EMGLViewQtWidget(self.parent))
			self.qwidgets[0].setQtWidget(self.fd)
			self.qwidgets[0].cam.setCamX(-100)
			self.initFlag = False
			

		
		for i in self.qwidgets:
			#print "getting opengl matrices"
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
			#print "paint child done"
		
		#print "paintGL done"
	def finished(self,val):
		if ( val == 1 ):
			for i in self.fd.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						w = EMGLView3D(self,a)
						self.qwidgets.append(w)
					else:
						w = EMGLView2D(self,a)
						self.qwidgets.append(w)
				else:
					w = EMGLView2D(self,a)
					self.qwidgets.append(w)
					
	def timer(self):
		pass
		#self.updateGL()
		
	def bindTexture(self,pixmap):
		return self.parent.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.parent.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.parent.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def resizeEvent(self, width, height):
		for i in self.qwidgets:
			i.set_update_P_inv()
	
	def mousePressEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
				self.updateGL()
				return
	
	def mouseMoveEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				self.current = i
				if (self.current != self.previous ):
					if ( self.previous != None ):
						self.previous.leaveEvent()
				i.mouseMoveEvent(event)
				self.previous = i
				self.updateGL()
				return
		
	def mouseReleaseEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
					
		
	def mouseDoubleClickEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
		
		
	def wheelEvent(self, event):
		for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.wheelEvent(event)
					self.updateGL()
					return

	def toolTipEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.toolTipEvent(event)
				self.updateGL()
				return
		
		QtGui.QToolTip.hideText()
		

	def dragMoveEvent(self,event):
		print "received drag move event"
		
	def event(self,event):
		#print "event"
		#QtGui.QToolTip.hideText()
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		#print "hoverEvent"
		if self.inspector :
			for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.hoverEvent(event)
					break
		self.updateGL()

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMFloatingWidgets()
	window2 = EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())

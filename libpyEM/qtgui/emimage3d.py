#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import Numeric
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary

class EMImage3D(QtOpenGL.QGLWidget):
	""" This class is not yet complete.
	A QT widget for rendering 3D EMData objects.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage2D.allim[self]=0
		
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

		
		self.data=None
		
		self.shapes={}
		self.shapechange=1
		
		self.inspector=None
		
		if image : 
			self.setData(image)
			self.show()
	
	def setData(self,data):
		"""Pass in a 3D EMData object"""
#		if not self.data and data: self.resize(data.get_xsize(),data.get_ysize())
		
		self.data=data
		if data==None:
			self.updateGL()
			return
		
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		
		self.minden=max(m0,mean-3.0*sigma)
		self.maxden=min(m1,mean+3.0*sigma)
		self.mindeng=max(m0,mean-5.0*sigma)
		self.maxdeng=min(m1,mean+5.0*sigma)
			
		self.showInspector()		# shows the correct inspector if already open
		self.updateGL()
		
	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		
		if self.shapechange:
			self.setupShapes()
			self.shapechange=0
		
		
		
		self.changec=self.data.get_attr("changecount")
				
	def resizeGL(self, width, height):
		side = min(width, height)
# 		GL.glViewport(0,0,self.width(),self.height())
# 	
# 		GL.glMatrixMode(GL.GL_PROJECTION)
# 		GL.glLoadIdentity()
# 		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
# 		GL.glMatrixMode(GL.GL_MODELVIEW)
# 		GL.glLoadIdentity()
		
	def setupShapes(self):
		# make our own cirle rather than use gluDisk or somesuch
		pass
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector3D(self)
		self.inspector.show()
	
	
	def mousePressEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
	
	def mouseMoveEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
# 		if self.rmousedrag and event.buttons()&Qt.RightButton:
# 			self.origin=(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
# 			self.rmousedrag=(event.x(),event.y())
# 			self.update()
		if event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
	
	def mouseReleaseEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
# 		if event.button()==Qt.RightButton:
# 			self.rmousedrag=None
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mouseup"), event)
				return


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
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vbl.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vbl.addWidget(self.maxs)
		
		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		self.vbl.addWidget(self.brts)
		
		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.vbl.addWidget(self.conts)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)
		QtCore.QObject.connect(self.invtog, QtCore.SIGNAL("toggled(bool)"), target.setInvert)
		QtCore.QObject.connect(self.ffttog, QtCore.SIGNAL("toggled(bool)"), target.setFFT)
		QtCore.QObject.connect(self.mmode, QtCore.SIGNAL("buttonClicked(int)"), target.setMMode)

	def newMin(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMin(val)

		self.updBC()
		self.busy=0
		
	def newMax(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMax(val)
		self.updBC()
		self.busy=0
	
	def newBrt(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0
		
	def newCont(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0

	def updBC(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		self.mins.setValue(x0)
		self.maxs.setValue(x1)
		self.target.setDenRange(x0,x1)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setLimits(self,lowlim,highlim,curmin,curmax):
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage3D()
 	if len(sys.argv)==1 : 
 		window.setData(test_image(size=(512,512)))

		# these lines are for testing shape rendering
# 		window.addShape("a",["rect",.2,.8,.2,20,20,80,80,2])
# 		window.addShape("b",["circle",.5,.8,.2,120,50,30.0,2])
# 		window.addShape("c",["line",.2,.8,.5,20,120,100,200,2])
# 		window.addShape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	else :
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
	

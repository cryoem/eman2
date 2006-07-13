from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import array

class EMImage(QtOpenGL.QGLWidget):
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.data=None
		self.scale=1.0
		self.minden=0
		self.maxden=1.0
		self.origin=(0,0)
		
		self.inspector=None
	
	def setData(self,data):
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
		
		if self.inspector: self.inspector.setLimits(self.minden,self.maxden,self.minden,self.maxden)
		
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setOrigin(self,x,y):
		self.origin=(x,y)
		self.updateGL()
		
	def setScale(self,val):
		self.scale=val
		self.updateGL()
		
	def setDenMin(self,val):
		self.minden=val
		self.updateGL()
		
	def setDenMax(self,val):
		self.maxden=val
		self.updateGL()

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
		GL.glLoadIdentity()
		GL.glTranslated(0.0, 0.0, -10.0)
		
		if self.data :
			a=self.data.render_amp8(self.origin[0],self.origin[1],self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
#			a=self.data.render_amp8(self.origin[0],self.origin[1],self.width(),self.height(),(self.width()-1)/4*4+4,1.0,1,254,self.minden,self.maxden,0)
#			GL.glPixelZoom(self.scale,self.scale)
			GL.glDrawPixels(self.width(),self.height(),GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
#			print self.width(),self.height(),len(a)
			hist=array.array("I")
			hist.fromstring(a[-1024:])
			if self.inspector : self.inspector.setHist(hist)
	
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport((width - side) / 2, (height - side) / 2, side, side)
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GL.glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
		GL.glMatrixMode(GL.GL_MODELVIEW)
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			if not self.inspector : self.inspector=EMImageInspector(self)
			if self.inspector: self.inspector.setLimits(self.minden,self.maxden,self.minden,self.maxden)
			self.inspector.show()
	
	def mouseMoveEvent(self, event):
		pass
	
	def mouseReleaseEvent(self, event):
		pass

class Histogram(QtGui.QWidget):
	def __init__(self,parent):
		QtGui.QWidget.__init__(self,parent)
		self.brush=QtGui.QBrush(Qt.black)
		
		self.histdata=None
		self.setMinimumSize(QtCore.QSize(258,128))
	
	def setData(self,data):
		self.histdata=data
		self.norm=max(self.histdata[1:-1])
		self.update()
	
	def paintEvent (self, event):
		if not self.histdata : return
		p=QtGui.QPainter()
		p.begin(self)
		p.setPen(Qt.black)
		for i,j in enumerate(self.histdata):
			p.drawLine(i,127,i,127-j*126/self.norm)
		p.end()

class EMImageInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vboxlayout = QtGui.QVBoxLayout(self)
		self.vboxlayout.setMargin(0)
		self.vboxlayout.setSpacing(6)
		self.vboxlayout.setObjectName("vboxlayout")
		
		self.hist = Histogram(self)
		self.hist.setObjectName("hist")
		self.vboxlayout.addWidget(self.hist)
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.vboxlayout.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vboxlayout.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vboxlayout.addWidget(self.maxs)
		
		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		self.vboxlayout.addWidget(self.brts)
		
		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.vboxlayout.addWidget(self.conts)
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)

	def newMin(self,val):
		self.target.setDenMin(val)
		self.minval=val
		
	def newMax(self,val):
		self.target.setDenMax(val)
		self.maxval=val
	
	def newBrt(self,val):
		pass
		
	def newCont(self,val):
		pass

	def setHist(self,hist):
		self.hist.setData(hist)

	def setLimits(self,lowlim,highlim,curmin,curmax):
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage()
	window.data=test_image()
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
	

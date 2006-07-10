from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys

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
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setDenOrigin(self,x,y):
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
			a=self.data.render_amp8(self.origin[0],self.origin[1],self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,1,254,self.minden,self.maxden,0)
#			print "%d %d   %d"%(self.width(),self.height(),len(a))
			GL.glDrawPixels(self.width(),self.height(),GL.GL_LUMINANCE,GL.GL_BYTE,a)
	
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport((width - side) / 2, (height - side) / 2, side, side)
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GL.glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
		GL.glMatrixMode(GL.GL_MODELVIEW)
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			if self.inspector : self.inspector.show()
			else:
				self.inspector=EMImageInspector(self)
				self.inspector.show()
	
	def mouseMoveEvent(self, event):
		pass
	
	def mouseReleaseEvent(self, event):
		pass
	
class EMImageInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		
		self.vboxlayout = QtGui.QVBoxLayout(self)
		self.vboxlayout.setMargin(0)
		self.vboxlayout.setSpacing(6)
		self.vboxlayout.setObjectName("vboxlayout")
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.vboxlayout.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vboxlayout.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vboxlayout.addWidget(self.maxs)
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), target.setDenMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), target.setDenMax)


	def setlims(target):
		data=target.data
		if data==None : return
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		
		self.mins.setRange(max(m0,mean-3.0*sigma),min(m1,mean+3.0*sigma))
		

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
	

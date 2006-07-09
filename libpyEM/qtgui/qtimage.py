from PyQt4 import QtCore, QtGui, QtOpenGL
from OpenGL import GL
from math import *
from EMAN2 import *
from valslider import ValSlider
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
	
	def setData(data):
		self.data=data
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		
		if data==None:
			self.updateGL()
			return
		
		self.minden=max(m0,mean-3.0*sigma)
		self.maxden=min(m1,mean+3.0*sigma)
		self.updateGL()
		
	def setrange(x0,x1):
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setorigin(x,y):
		self.origin=(x,y)
		self.updateGL()
		
	def initializeGL(self):
		GL.glClearColor(0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
		GL.glLoadIdentity()
		GL.glTranslated(0.0, 0.0, -10.0)
		
		if self.data :
			a=self.data.render_amp8(origin[0],origin[1],self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,1,254,self.minden,self.maxden,0)
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
		pass
	
	def mouseMoveEvent(self, event):
		pass
	
	def mouseReleaseEvent(self, event):
		pass

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = QtGLImage()
	window.data=test_image()
	window.show()
	
	w2=QtGui.QWidget()
	w2.resize(256,128)
	
	w3=ValSlider(w2)
	w3.resize(256,24)
	w2.show()
	
	sys.exit(app.exec_())
	

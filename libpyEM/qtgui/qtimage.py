from PyQt4 import QtCore, QtGui, QtOpenGL
from OpenGL import GL
from math import *
from EMAN2 import *
import sys

class QtGLImage(QtOpenGL.QGLWidget):
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		self.data=None
	
	def initializeGL(self):
		GL.glClearColor(.5,0.,.5,0.)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
		GL.glLoadIdentity()
		GL.glTranslated(0.0, 0.0, -10.0)
		
		if self.data :
			a=self.data.render_amp8(0,0,self.width(),self.height(),(self.width()-1)/4*4+4,1.0,1,254,self.data.get_attr("minimum"),self.data.get_attr("maximum"),0)
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
	

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = QtGLImage()
	window.data=test_image()
	window.show()
	sys.exit(app.exec_())
	
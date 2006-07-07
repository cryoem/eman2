from PyQt4 import QtCore, QtGui, QtOpenGL
from OpenGL import GL
from math import *
import sys

class QtGLImage(QtOpenGL.QGLWidget):
	def __init__(self, parent=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
	
	def initializeGL(self):
		GL.glClearColor(.5,0.,.5,0.)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
		GL.glLoadIdentity()
		GL.glTranslated(0.0, 0.0, -10.0)
	
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
    window.show()
    sys.exit(app.exec_())
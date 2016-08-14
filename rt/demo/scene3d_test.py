#!/usr/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from emscene3d import EMScene3D, EMInspector3D, EMInspectorControlShape
from emshapeitem3d import *

from EMAN2 import *

class GLdemo(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.widget = EMScene3D()

		self.cube1 = EMCube(50.0)
		self.widget.insertNewNode("Cube1", self.cube1)    # Something to Render something..... an EMItem3D

		self.cube2 = EMCube(50.0)
		self.widget.insertNewNode("Cube2", self.cube2)
		
		self.sphere = EMSphere(50.0)
		self.widget.insertNewNode("Sphere", self.sphere)
		
		self.cylinder = EMCylinder(50, 300)
		self.widget.insertNewNode("Cylinder", self.cylinder)
		
		self.cone = EMCone(50, 300)
		self.widget.insertNewNode("Cone", self.cone)
		
		self.line = EMLine(10,10,10,200,200,200, 16)
		self.widget.insertNewNode("Line", self.line)
		
		self.text1 = EM3DText('3D text', 75)
		self.widget.insertNewNode("Text1", self.text1)
		
		
		#Show inspector
		#self.widget.showInspector()
		
		# QT stuff to display the widget
		vbox = QtGui.QVBoxLayout()
		vbox.addWidget(self.widget)
		self.setLayout(vbox)
		self.setGeometry(300, 300, 600, 600)
		self.setWindowTitle('BCM EM Viewer')
	
	def show_inspector(self):
		self.inspector.show()
		
if __name__ == "__main__":
	import sys
	app = QtGui.QApplication(sys.argv)
	window = GLdemo()
	window.show()
	app.exec_()
	

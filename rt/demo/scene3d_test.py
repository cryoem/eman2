#!/usr/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from emscene3d import EMScene3D, EMInspector3D, EMInspectorControlShape
from emshapeitem3d import *

from EMAN2 import *

class GLdemo(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.widget = EMScene3D()
		#self.widget.camera.useprespective(50, 0.5)
#		self.cube1 = glCube(50.0)
		self.cube1 = EMCube(50.0)
		self.widget.addChild(self.cube1)    # Something to Render something..... (this could just as well be one of Ross's SGnodes)
		#self.widget.activatenode(cube1)
#		self.cube2 = glCube(50.0)
		self.cube2 = EMCube(50.0)
		self.widget.addChild(self.cube2)
		#self.widget.activatenode(cube2)

		self.inspector = EMInspector3D(self.widget)
		self.widget.set_inspector(self.inspector)
		
		rootnode = self.inspector.add_tree_node("root node", self.widget)
		self.inspector.add_tree_node("cube1", self.cube1, rootnode)
		self.inspector.add_tree_node("cube2", self.cube2, rootnode)
		
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
	window.show_inspector()
	app.exec_()
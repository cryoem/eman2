#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import sys

class Shell(QtGui.QTextEdit):
	
	def __init__(self):
		QtGui.QTextEdit.__init__(self)
	
	def keyPressEvent(self,event):
#		print self.textCursor().position()
		print self.document().findBlock(4).text()
#		print self.toPlainText()
		if self.find("\0012") : 
			cur=self.textCursor()
			cur.movePosition(QtGui.QTextCursor.End)
			self.setTextCursor(cur)
		QtGui.QTextEdit.keyPressEvent(self,event)

if __name__ == '__main__':
	app = QtGui.QApplication([])
	window = Shell()
	window.show()
	
	sys.exit(app.exec_())

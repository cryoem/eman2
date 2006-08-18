#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import sys
import imp
from subprocess import Popen,PIPE
import traceback
from EMAN2 import *
from emimagemx import *
from emimage2d import *

class Shell(QtGui.QTextEdit):
	
	def __init__(self):
		QtGui.QTextEdit.__init__(self)
		self.setPlainText("Welcome to the EMAN2 interactive Python environment\nversion0.01 alpha\n\n>>> ")
		
		self.local={}
#		self.connect(self,QtCore.SIGNAL("CursorPositionChanged"), self.newCursor)
#		self.python=Popen("python",stdin=PIPE,stdout=PIPE,stderr=PIPE,universal_newlines=True)
	
#	def newCursor(self):
#		print self.textCursor().position()
	
	def keyPressEvent(self,event):
		lpos=self.document().end().previous().position()
		cur=self.textCursor()
		if cur.position()<lpos+4: 
			cur.movePosition(QtGui.QTextCursor.End)
			self.setTextCursor(cur)
		if cur.position()==lpos+4 and event.key()==Qt.Key_Backspace : return
		
		if event.key()==Qt.Key_Return :
			com=str(self.document().end().previous().text()[4:])
#			print com
			try: 
				ret=eval(com,globals(),self.local)
			except SyntaxError:
				pass
			except:
				self.append(traceback.format_exc()+">>> ")
				return
			else:
				self.append(repr(ret)+"\n>>> ")
				return
			
			try:
				exec com in globals(),self.local
			except:
				self.append(traceback.format_exc()+">>> ")
			else:
				self.append(">>> ")
			return
			
			
		QtGui.QTextEdit.keyPressEvent(self,event)

if __name__ == '__main__':
	app = QtGui.QApplication([])
	window = Shell()
	window.show()
	
	sys.exit(app.exec_())

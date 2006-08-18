#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import sys
import imp
from subprocess import Popen,PIPE
import traceback
from EMAN2 import *
from emimage2d import *
from emimagemx import *

class Shell(QtGui.QTextEdit):
	
	def __init__(self):
		QtGui.QTextEdit.__init__(self)
		self.setFontFamily("Courier")
		self.setFontPointSize(12)
		self.setPlainText("Welcome to the EMAN2 interactive Python environment\nversion0.01 alpha\n\n>>> ")
		
		self.com=""
		self.local={}
#		self.connect(self,QtCore.SIGNAL("CursorPositionChanged"), self.newCursor)
#		self.python=Popen("python",stdin=PIPE,stdout=PIPE,stderr=PIPE,universal_newlines=True)
	
#	def newCursor(self):
#		print self.textCursor().position()
	
	def keyPressEvent(self,event):
		"This is where most of the magic that allows the user to only type on the current line is"
		lpos=self.document().end().previous().position()
		cur=self.textCursor()
		if cur.position()<lpos+4: 
			cur.movePosition(QtGui.QTextCursor.End)
			self.setTextCursor(cur)
		if cur.position()==lpos+4 and event.key()==Qt.Key_Backspace : return
		
		if event.key()==Qt.Key_Return :
			self.com+=str(self.document().end().previous().text()[4:])
# 			print self.com
			if self.com.count("\n")>0 and self.com[-1]!="\n" : 
				self.com+="\n"
				self.append("--> ")
				return
			
			try: 
				ret=eval(self.com,globals(),self.local)
			except SyntaxError:
				pass
			except:
				self.append(traceback.format_exc()+">>> ")
				self.com=""
				return
			else:
				self.append(repr(ret)+"\n>>> ")
				self.com=""
				return

			try:
				exec self.com in globals(),self.local
			except SyntaxError,s:
				if s[0][:14]=="unexpected EOF" :
					self.com+="\n"
					self.append("--> ")
					return
				self.append(traceback.format_exc()+">>> ")
				self.com=""
			except:
				self.append(traceback.format_exc()+">>> ")
				self.com=""
			else:
				self.com=""
				self.append(">>> ")
			return
			
			
		QtGui.QTextEdit.keyPressEvent(self,event)
	
	def insertFromMimeData(self,source):
		"This deals with 'paste' events"	

if __name__ == '__main__':
	app = QtGui.QApplication([])
	window = Shell()
	window.show()
	
	sys.exit(app.exec_())

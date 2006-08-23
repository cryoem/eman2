#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import sys
import imp
from subprocess import Popen,PIPE
import traceback
from EMAN2 import *
from emimage import *

class Shell(QtGui.QTextEdit):
	
	def __init__(self):
		QtGui.QTextEdit.__init__(self)
		self.setFontFamily("Courier")
		self.setFontPointSize(12)
		self.setPlainText("Welcome to the EMAN2 interactive Python environment\nversion0.01 alpha\n\n>>> ")
		self.history=[]
		self.histc=0
		
		self.com=""
		self.local={}
#		self.connect(self,QtCore.SIGNAL("CursorPositionChanged"), self.newCursor)
#		self.python=Popen("python",stdin=PIPE,stdout=PIPE,stderr=PIPE,universal_newlines=True)
	
#	def newCursor(self):
#		print self.textCursor().position()
	
	def keyPressEvent(self,event):
		"This is where most of the magic that allows the user to only type on the current line is"
		if event.key()==Qt.Key_Up :
			cur=self.textCursor()
			self.histc-=1
			cur.movePosition(QtGui.QTextCursor.End)
			cur.select(QtGui.QTextCursor.BlockUnderCursor)
			self.setTextCursor(cur)
			cur.deleteChar()
			try:
				self.append(">>> "+self.history[self.histc])
			except:
				self.append(">>> "+self.history[0])
				self.histc+=1
			return
		if event.key()==Qt.Key_Down:
			cur=self.textCursor()
			self.histc+=1
			if self.histc>0: self.histc=0
			cur.movePosition(QtGui.QTextCursor.End)
			cur.select(QtGui.QTextCursor.BlockUnderCursor)
			self.setTextCursor(cur)
			cur.deleteChar()
			if self.histc==0 : self.append(">>> ")
			else:
				try:
					self.append(">>> "+self.history[self.histc])
				except:
					self.append(">>> ")
					self.histc=0
			return
		lpos=self.document().end().previous().position()
		cur=self.textCursor()
		if cur.position()<lpos+4: 
			cur.movePosition(QtGui.QTextCursor.End)
			self.setTextCursor(cur)
		if cur.position()==lpos+4 and event.key()==Qt.Key_Backspace : return
		
		# This is where the code is actually executed
		if event.key()==Qt.Key_Return :
			self.history.append(str(self.document().end().previous().text()[4:]))
			self.histc=0
			self.com+=str(self.document().end().previous().text()[4:])
# 			print self.com
			# for multi-line entries, a blank line is required before we execute
			if self.com.count("\n")>0 and self.com[-1]!="\n" : 
				self.com+="\n"
				self.append("--> ")
				return
			
			# try executing it as an expression
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

			# expression didn't work, now try it as a statement
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
		lns=str(source.text().replace(">>> ","")).split("\n")
		for i in lns: 
			self.append("P>> "+i)
			self.history.append(i)
			self.histc=0
		self.com+="\n".join(lns)
		

if __name__ == '__main__':
	app = QtGui.QApplication([])
	window = Shell()
	window.show()
	
	sys.exit(app.exec_())

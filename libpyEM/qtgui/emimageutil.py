#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import Qt
from math import *
import numpy

class EventRerouter:
	def __init__(self,target=None):
		self.target = target
		
	def set_target(self,target):
		self.target = target

	def mousePressEvent(self, event):
		self.target.mousePressEvent(event)
			
	def wheelEvent(self,event):
		self.target.wheelEvent(event)
	
	def mouseMoveEvent(self,event):
		self.target.mouseMoveEvent(event)

	def mouseReleaseEvent(self,event):
		self.target.mouseReleaseEvent(event)
		
	def mouseDoubleClickEvent(self,event):
		self.target.mouseDoubleClickEvent(event)
		
	def keyPressEvent(self,event):
		self.target.keyPressEvent(event)
		
	def dropEvent(self,event):
		self.target.dropEvent(event)
		
	def closeEvent(self,event) :
		self.target.closeEvent(event)
		
	def dragEnterEvent(self,event):
		self.target.dragEnterEvent(event)
		
	def keyPressEvent(self,event):
		self.target.keyPressEvent(event)

	def get_core_object(self):
		return self.target
	
	def get_target(self):
		return self.target # use this one instead of the above

class EMParentWin(QtGui.QWidget):
	"""Used to give the opengl widgets a parent, necessary for OSX Leopard"""
	def __init__(self,child=None):
		QtGui.QWidget.__init__(self,None)
		
		self.child = child
		self.resize(child.width()+20,child.height()+20)
		self.setMaximumSize(8000,8000)

		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(4)
		self.hbl.setSpacing(6)
		self.hbl.addWidget(self.child)
		self.setLayout(self.hbl)
		
#		self.setWindowTitle(self.tr("Whatever"))

	def closeEvent(self, e):
		try:
			self.child.closeEvent(e)
			#self.child.inspector.close()
		except: pass
		QtGui.QWidget.closeEvent(self,e)
	
	def get_qt_widget(self):
		return self.child
	
class ImgHistogram(QtGui.QWidget):
	""" A small fixed-size histogram widget"""
	def __init__(self,parent):
		QtGui.QWidget.__init__(self,parent)
		self.brush=QtGui.QBrush(Qt.black)
		
		self.font=QtGui.QFont("Helvetica", 12);
		self.probe=None
		self.histdata=None
		self.setMinimumSize(QtCore.QSize(258,128))
	
	def set_data(self,data,minden,maxden):
		self.histdata=data
#		self.norm=max(self.histdata)
		self.norm=0
		self.minden=minden
		self.maxden=maxden
		for i in self.histdata: self.norm+=float(i)*i
		self.norm-=max(self.histdata)**2
		self.norm=sqrt(self.norm/255)*3.0
		self.total=sum(self.histdata)
		if self.norm==0 : self.norm=1.0
		if self.total==0 : self.histdata=None
		self.update()
	
	def paintEvent (self, event):
		if self.histdata==None : return
		p=QtGui.QPainter()
		p.begin(self)
		p.setBackground(QtGui.QColor(16,16,16))
		p.eraseRect(0,0,self.width(),self.height())
		p.setPen(Qt.darkGray)
		for i,j in enumerate(self.histdata):
			p.drawLine(i,127,i,127-j*126/self.norm)
		
		# If the user has dragged, we need to show a value
		if self.probe :
			p.setPen(Qt.blue)
			p.drawLine(self.probe[0]+1,0,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setPen(Qt.red)
			p.drawLine(self.probe[0]+1,127,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setFont(self.font)
			p.drawText(200,20,"x=%d"%(self.probe[0]))
			p.drawText(200,36,"%1.2f"%(self.probe[0]/255.0*(self.maxden-self.minden)+self.minden))
			p.drawText(200,52,"y=%d"%(self.probe[1]))
			p.drawText(200,68,"%1.2f%%"%(100.0*float(self.probe[1])/self.total))
		
		p.setPen(Qt.black)
		p.drawRect(0,0,257,128)
		p.end()

	def mousePressEvent(self, event):
		if event.button()==Qt.LeftButton:
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.update()
			
	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton:
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.update()
	
	def mouseReleaseEvent(self, event):
		self.probe=None

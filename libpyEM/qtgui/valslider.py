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

import sys
from PyQt4 import QtCore, QtGui

def clamp(x0,val,x1):
	return max(min(val,x1),x0)

class ValSlider(QtGui.QWidget):
	"""The valslider class represents a connected text widget and horizontal slider.
	setValue(float) - to programatically change the value
	emit valueChanged(float)
	"""
	def __init__(self, parent, range=None, label=None, value=0,labelwidth=30):
		if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		if range : self.range=list(range)
		else : self.range=[0,1.0]
		self.value=value
		self.ignore=0
		self.intonly=0
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if label:
			self.label = QtGui.QLabel(self)
			self.setLabel(label)
			
			sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(0),QtGui.QSizePolicy.Policy(0))
#			sizePolicy.setHorizontalStretch(1)
#			sizePolicy.setVerticalStretch(0)
#			sizePolicy.setHeightForWidth(self.text.sizePolicy().hasHeightForWidth())
			self.label.setSizePolicy(sizePolicy)
#			self.label.setAlignment(QtCore.Qt.AlignRight+QtCore.Qt.AlignVCenter)
			self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
			self.label.setMinimumSize(QtCore.QSize(labelwidth,20))
			self.label.setObjectName("label")
			self.hboxlayout.addWidget(self.label)
		
		
		self.text = QtGui.QLineEdit(self)
		
		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(0))
		sizePolicy.setHorizontalStretch(1)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.text.sizePolicy().hasHeightForWidth())
		self.text.setSizePolicy(sizePolicy)
		self.text.setMinimumSize(QtCore.QSize(80,0))
		self.text.setObjectName("text")
		self.hboxlayout.addWidget(self.text)
		
		self.slider = QtGui.QSlider(self)
		
		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(0))
		sizePolicy.setHorizontalStretch(7)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.slider.sizePolicy().hasHeightForWidth())
		self.slider.setSizePolicy(sizePolicy)
		self.slider.setMinimumSize(QtCore.QSize(100,0))
		self.slider.setMaximum(4095)
		self.slider.setSingleStep(16)
		self.slider.setPageStep(256)
		self.slider.setOrientation(QtCore.Qt.Horizontal)
		self.slider.setObjectName("slider")
		self.hboxlayout.addWidget(self.slider)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		QtCore.QObject.connect(self.slider, QtCore.SIGNAL("valueChanged(int)"), self.sliderChange)
		QtCore.QObject.connect(self.slider, QtCore.SIGNAL("sliderReleased()"), self.sliderReleased)
		QtCore.QObject.connect(self.slider, QtCore.SIGNAL("sliderPressed()"), self.sliderPressed)

	def setRange(self,minv,maxv):
		if maxv<=minv : maxv=minv+.001
		self.range=[float(minv),float(maxv)]
		self.updates()
		#self.slider.setRange(*self.range)
		
	def getRange(self): return self.range

	def setValue(self,val,quiet=0):
		if val <= self.range[0]:
			self.range[0] = val
			#self.updates()
		if val >= self.range[1]:
			self.range[1] = val
			#self.updates()
		
		if self.intonly : 
			if self.value==int(val+.5) : return
			self.value=int(val+.5)
		else: 
			if self.value==val : return
			self.value=val

		self.updateboth()
		if not quiet : self.emit(QtCore.SIGNAL("valueChanged"),self.value)
	
	def getValue(self):
		return self.value
	
	def setIntonly(self,flag):
		self.intonly=flag
		self.updateboth()
		
	def textChange(self):
		if self.ignore : return
		x=self.text.text()
		if len(x)==0 : return
		if x[0]=='<' : 
			try: self.range[1]=float(x[1:])
			except: pass
			self.updateboth()
		elif x[0]=='>' : 
			try: self.range[0]=float(x[1:])
			except: pass
			self.updateboth()
		else:
			try:
				if (self.intonly) :
					if self.value==int(float(x)+.5) : return
					self.value=int(float(x)+.5)
				else : 
					self.value=float(x)
				self.updates()
				self.emit(QtCore.SIGNAL("valueChanged"),self.value)
				self.emit(QtCore.SIGNAL("textChanged"),self.value)
			except:
				self.updateboth()
				
	def sliderChange(self,x):
		if self.ignore : return
		ov=self.value
		self.value=(self.slider.value()/4095.0)*(self.range[1]-self.range[0])+self.range[0]
		if self.intonly : 
			self.value=int(self.value+.5)
			if self.value==ov : return
		self.updatet()
		self.emit(QtCore.SIGNAL("valueChanged"),self.value)
	
	def sliderReleased(self):
		self.emit(QtCore.SIGNAL("sliderReleased"),self.value)
		
	def sliderPressed(self):
		self.emit(QtCore.SIGNAL("sliderPressed"),self.value)
	
	def setLabel(self,label):
		self.label.setText(label)
	
	def getLabel(self):
		return str(self.label.text())
		
	def updates(self):
		self.ignore=1
		self.slider.setValue(clamp(0,(self.value-self.range[0])/(self.range[1]-self.range[0])*4095.0,4095.0))
		self.ignore=0

	def updatet(self):
		self.ignore=1
		self.text.setText(str(self.value)[:self.text.width()/10-1])
		self.ignore=0
		
	def updateboth(self):
		self.updates()
		self.updatet()

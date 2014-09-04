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

import sys, math, weakref
from OpenGL.GL import *
from OpenGL import GLU
from PyQt4 import QtCore, QtGui, QtOpenGL 
from PyQt4.QtCore import Qt

leftarrow = [
    '12 10 2 1',
    'b c #0000ff',
    'c c None',
    'ccccccccbbbb',
    'ccccccbbbbbc',
    'ccccbbbbbccc',
    'ccbbbbbccccc',
    'bbbbbccccccc',
    'bbbbbccccccc',
    'ccbbbbbccccc',
    'ccccbbbbbccc',
    'ccccccbbbbbc',
    'ccccccccbbbb'
]    

rightarrow = [
    '12 10 2 1',
    'b c #0000ff',
    'c c None',
    'bbbccccccccc',
    'cbbbbbcccccc',
    'cccbbbbbcccc',
    'cccccbbbbbcc',
    'cccccccbbbbb',
    'cccccccbbbbb',
    'cccccbbbbbcc',
    'cccbbbbbcccc',
    'cbbbbbcccccc',
    'bbbbcccccccc'
]

def clamp(x0,val,x1):
	return int(max(min(val,x1),x0))

class ValSlider(QtGui.QWidget):
	"""The valslider class represents a connected text widget and horizontal slider.
	showenable - if -1, no enable box shown, if 0, shown unchecked, if 1 shown and checked
	setValue(float) - to programatically change the value
	emit valueChanged(float)
	"""
	def __init__(self, parent=None, rng=None, label=None, value=0,labelwidth=30,showenable=-1,rounding=3):
		#if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		#print label, "allocated"
		
		if rng : self.rng=list(rng)
		else : self.rng=[0,1.0]
		if value==None : value=0.0
		self.value=value
		self.oldvalue=value-1.0
		self.ignore=0
		self.intonly=0
		self.rounding = rounding
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.setEnabled)
		
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
		
		self.updateboth()
		if showenable>=0 : self.setEnabled(showenable)
		
	#def __del__(self):
		#print self.getLabel(), " freed"
#		QtGui.QWidget.__del__(self)

	def setEnabled(self,ena):
		self.slider.setEnabled(ena)
		self.text.setEnabled(ena)
		try : return self.enablebox.setChecked(ena)
		except: pass
		self.emit(QtCore.SIGNAL("enableChanged"),ena) 

	def getEnabled(self):
		try : return self.enablebox.isChecked()
		except: return True

	def setRange(self,minv,maxv):
		if maxv<=minv : maxv=minv+.001
		if maxv==minv : minv=math.fabs(minv)*0.999
		self.rng=[float(minv),float(maxv)]
		self.updates()
		#self.slider.setRange(*self.rng)
		
	def getRange(self): return self.rng

	def setValue(self,val,quiet=0):
		if val==None : return

		if val <= self.rng[0]:
			self.rng[0] = val
			#self.updates()
		if val >= self.rng[1]:
			self.rng[1] = val
			#self.updates()
		
		if self.intonly : 
			if self.value==int(val+.5) : return
			self.value=int(val+.5)
		else: 
			if self.value==val : return
			self.value=val

#		if self.intonly : print self.value,val
		self.updateboth()
		if not quiet and self.value!=self.oldvalue: self.emit(QtCore.SIGNAL("valueChanged"),self.value)
		self.oldvalue=self.value
	
	def getValue(self):
		return self.value
	
	def setIntonly(self,flag):
		self.intonly=flag
		if flag : 
			self.rounding=0
			self.text.setMinimumSize(QtCore.QSize(50,0))
		else:
			self.rounding=3
			self.text.setMinimumSize(QtCore.QSize(80,0))
			
		self.updateboth()
		
	def textChange(self):
		if self.ignore : return
		x=self.text.text()
		if len(x)==0 : return
		if x[0]=='<' : 
			try: self.rng[1]=float(x[1:])
			except: pass
			self.updateboth()
		elif x[0]=='>' : 
			try: self.rng[0]=float(x[1:])
			except: pass
			self.updateboth()
		else:
			try:
				if (self.intonly) :
					if self.value==int(float(x)+.5) : return
					self.value=int(float(x)+.5)
				else : 
					self.value=float(x)
#				print "new text ",self.value

				# If the user enters a value out of range, adjust the slider range
				if self.value>self.rng[1] : self.rng[1]=self.value
				if self.value<self.rng[0] : self.rng[0]=self.value
				
				self.updates()
				if self.value!=self.oldvalue:
					self.emit(QtCore.SIGNAL("valueChanged"),self.value)
					self.emit(QtCore.SIGNAL("textChanged"),self.value)
				self.oldvalue=self.value
			except:
				self.updateboth()
				
	def sliderChange(self,x):
		if self.ignore : return
		ov=self.value
		self.value=(self.slider.value()/4095.0)*(self.rng[1]-self.rng[0])+self.rng[0]
		if self.intonly : 
			self.value=int(self.value+.5)
			if self.value==ov : return
		self.updatet()
		if self.value!=self.oldvalue: self.emit(QtCore.SIGNAL("valueChanged"),self.value)
		self.oldvalue=self.value
	
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
#		print "updates: ",self.value,self.rng
		self.slider.setValue(clamp(0,(self.value-self.rng[0])/float(self.rng[1]-self.rng[0])*4095.0,4095.0))
		self.ignore=0

	def updatet(self):
		self.ignore=1
		#self.text.setText(str(self.value)[:self.text.width()/10-1]) # Stupido!!! This will turn a sci number to a real
		if self.intonly : self.text.setText(str(int(self.value)))
		else : self.text.setText(str(round(self.value,self.rounding)))
		self.ignore=0
		
	def updateboth(self):
		self.updates()
		self.updatet()

class ValBox(QtGui.QWidget):
	"""A ValSlider without the slider part. Everything is the same except that the slider doesn't exist,
	so for virtually all purposes it could be used as a drop-in replacement.
	"""
	def __init__(self, parent=None, rng=None, label=None, value=0,labelwidth=30,showenable=-1):
		#if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		if rng : self.rng=list(rng)
		else : self.rng=[0,1.0]
		if value==None : value=0
		self.value=value
		self.ignore=0
		self.intonly=0
		self.digits=5		# digits of precision to display
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.setEnabled)
			
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
		self.text.setMinimumSize(QtCore.QSize(60,0))
		self.text.setObjectName("text")
		self.hboxlayout.addWidget(self.text)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		
		if showenable>=0 : self.setEnabled(showenable)
		self.updateboth()

	def setEnabled(self,ena):
		self.text.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enableChanged"),ena) 

	def getEnabled(self):
		try : return self.enablebox.isChecked()
		except: return True
		
	def setRange(self,minv,maxv):
		if maxv<=minv : maxv=minv+.001
		self.rng=[float(minv),float(maxv)]
		self.updates()
		#self.slider.setRange(*self.rng)
		
	def getRange(self): return self.rng

	def setDigits(self,digits): 
		self.digits=int(digits)

	def setValue(self,val,quiet=0):
		if val==None: 
			self.text.setText("")
			return
		if val <= self.rng[0]:
			self.rng[0] = val
			#self.updates()
		if val >= self.rng[1]:
			self.rng[1] = val
			#self.updates()
		
		if self.intonly : 
			if self.value==int(val+.5) : return
			self.value=int(val+.5)
		else: 
			if self.value==val : return
			self.value=val

#		if self.intonly : print self.value,val
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
			try: self.rng[1]=float(x[1:])
			except: pass
			self.updateboth()
		elif x[0]=='>' : 
			try: self.rng[0]=float(x[1:])
			except: pass
			self.updateboth()
		else:
			try:
				if (self.intonly) :
					if self.value==int(float(x)+.5) : return
					self.value=int(float(x)+.5)
				else : 
					self.value=float(x)
#				print "new text ",self.value
				self.emit(QtCore.SIGNAL("valueChanged"),self.value)
				self.emit(QtCore.SIGNAL("textChanged"),self.value)
			except:
				self.updateboth()
				
	def setLabel(self,label):
		self.label.setText(label)
	
	def getLabel(self):
		return str(self.label.text())
		
	def updatet(self):
		self.ignore=1
		fmt="%%1.%dg"%self.digits
		self.text.setText((fmt%self.value))
		self.ignore=0
		
	def updateboth(self):
		self.updatet()

class StringBox(QtGui.QWidget):
	"""A ValBox but it takes arbitrary text. Basically maintains the label/enable functionality for a QLineEdit widget
	"""
	def __init__(self, parent=None, label=None, value="",labelwidth=30,showenable=-1):
		#if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		if value==None : value=""
		self.ignore=0
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.setEnabled)
			
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
		self.text.setText(value)
		
		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(0))
		sizePolicy.setHorizontalStretch(1)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.text.sizePolicy().hasHeightForWidth())
		self.text.setSizePolicy(sizePolicy)
		self.text.setMinimumSize(QtCore.QSize(80,0))
		self.text.setObjectName("text")
		self.hboxlayout.addWidget(self.text)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		
		if showenable>=0 : self.setEnabled(showenable)

	def setEnabled(self,ena):
		self.text.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enableChanged"),ena) 

	def getEnabled(self):
		try : return self.enablebox.isChecked()
		except: return True

	def setValue(self,val,quiet=0):
		if self.getValue()==val : return
		self.text.setText(val)
		if not quiet : self.emit(QtCore.SIGNAL("valueChanged"),val)
	
	def getValue(self):
		return str(self.text.text())
		
	def textChange(self):
		if self.ignore : return
		self.emit(QtCore.SIGNAL("valueChanged"),self.getValue())
		self.emit(QtCore.SIGNAL("textChanged"),self.getValue())
				
	def setLabel(self,label):
		self.label.setText(label)
	
	def getLabel(self):
		return str(self.label.text())
		

class CheckBox(QtGui.QWidget):
	"""A QCheckBox with a label
	"""
	def __init__(self, parent=None, label=None, value="",labelwidth=30,showenable=-1):
		#if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		if value==None : value=False
		if value!=False and value!=True : value=True
		
		self.ignore=0
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.setEnabled)
			
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
		
		
		self.check = QtGui.QCheckBox(self)
#		print "'%s'"%str(value)
		self.check.setChecked(value)
		self.hboxlayout.addWidget(self.check)
		
		QtCore.QObject.connect(self.check, QtCore.SIGNAL("stateChanged(int)"), self.boolChanged)
		
		if showenable>=0 : self.setEnabled(showenable)

	def setEnabled(self,ena):
		self.check.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enableChanged"),ena) 
		
	def getEnabled(self):
		try : return self.enablebox.isChecked()
		except: return True
		
	def setValue(self,val,quiet=0):
		try:
			if int(val)!=0: val=True
			else : val=False
		except:
			if isinstance(val,str) : 
				if val.lower() in ("true","yes","t","y") : val=True
				else : val=False
				
		if self.getValue()==val : return
		self.check.setChecked(val)
		if not quiet : self.emit(QtCore.SIGNAL("valueChanged"),val)
	
	def getValue(self):
		return bool(self.check.isChecked())
	
		
	def boolChanged(self,newv):
		if self.ignore : return
		self.emit(QtCore.SIGNAL("valueChanged"),bool(self.check.isChecked()))
				
	def setLabel(self,label):
		self.label.setText(label)
	
	def getLabel(self):
		return str(self.label.text())


class RangeSlider(QtGui.QWidget):
	"""This is an int slider with two values in a fixed range (v0,v1) in a fixed range (min,max). Each value
	can be set individually or the pair can be moved up and down together. The values are displayed at
	the top and bottom of the vertical slider.
	"""
	def __init__(self, parent=None, rng=(0,100), value=(25,75)):
		#if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		self.rng=tuple(rng)
		self.value=tuple(value)
		self.mdownloc=None
		if len(rng)!=2 or len(value)!=2 : raise Exception,"RangeSlider needs a valid range and value)"

		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(0),QtGui.QSizePolicy.Policy(7))
		sizePolicy.setHorizontalStretch(0)
		sizePolicy.setVerticalStretch(7)
		sizePolicy.setHeightForWidth(False)
		self.setSizePolicy(sizePolicy)
#		self.setMinimumSize(QtCore.QSize(12,80))

		#self.vboxlayout = QtGui.QVBoxLayout(self)
		#self.vboxlayout.setMargin(0)
		#self.vboxlayout.setSpacing(6)
		#self.vboxlayout.setObjectName("vboxlayout")
		
		## Label on top
		#self.toptxt = QtGui.QLabel(str(value[1]))
		#self.vboxlayout.addWidget(self.toptxt)
		
		## canvas in the middle
		#self.slider = QtGui.QWidget(self)	
		#sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(0),QtGui.QSizePolicy.Policy(7))
		#sizePolicy.setHorizontalStretch(0)
		#sizePolicy.setVerticalStretch(7)
		#sizePolicy.setHeightForWidth(False)
		#self.slider.setSizePolicy(sizePolicy)
		#self.slider.setMinimumSize(QtCore.QSize(15,100))
		#self.vboxlayout.addWidget(self.slider)
		
		## label on the bottom
		#self.bottxt = QtGui.QLabel(str(value[0]))
		#self.vboxlayout.addWidget(self.bottxt)
		
		#self.text = QtGui.QLineEdit(self)
		
		
		
		#QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		#QtCore.QObject.connect(self.slider, QtCore.SIGNAL("valueChanged(int)"), self.sliderChange)
		#QtCore.QObject.connect(self.slider, QtCore.SIGNAL("sliderReleased()"), self.sliderReleased)
		#QtCore.QObject.connect(self.slider, QtCore.SIGNAL("sliderPressed()"), self.sliderPressed)
		
		#self.updateboth()

	def sizeHint(self): return QtCore.QSize(15,100)

	def paintEvent(self,event):
		"""Redraws the widget"""
		p=QtGui.QPainter(self)
		p.setPen(Qt.gray)
		p.drawRect(1,1,self.size().width()-2,self.size().height()-2)
		p.setPen(Qt.black)
		p.drawRect(0,0,self.size().width()-2,self.size().height()-2)
		p.setPen(Qt.blue)
		p.drawLine(3,self.vtoy(self.value[0]),self.size().width()-4,self.vtoy(self.value[0]))
		p.drawLine(3,self.vtoy(self.value[1]),self.size().width()-4,self.vtoy(self.value[1]))
		p.drawLine(self.size().width()/2,self.vtoy(self.value[0]),self.size().width()/2,self.vtoy(self.value[1]))
		
	def mousePressEvent(self,event):
		y=event.y()
		v0y=self.vtoy(self.value[0])
		v1y=self.vtoy(self.value[1])
		
		print y,v0y,v1y
		
		# outside the current range, no effect
		if y>v0y+3 : return
		if y<v1y-3 : return
		
		# middle region, drag both v0 and v1
		if (v0y-v1y<6 and y<v0y and y>v1y) or (y<v0y-3 and y>v1y+3) :
			self.mdownloc=(3,event.x(),event.y(),self.value[0],self.value[1])
		
		# drag v0
		elif (y>=v0y-3) :
			self.mdownloc=(1,event.x(),event.y(),self.value[0],self.value[1])
		
		# drag v1
		else :
			self.mdownloc=(2,event.x(),event.y(),self.value[0],self.value[1])
			
	def mouseMoveEvent(self,event):
		y=event.y()
		if self.mdownloc!=None :
			v0=self.mdownloc[3]
			v1=self.mdownloc[4]
			
			if self.mdownloc[0]&1 :
				v0=self.ytov(y-self.mdownloc[2]+self.vtoy(self.mdownloc[3]))

			if self.mdownloc[0]&2 :
				v1=self.ytov(y-self.mdownloc[2]+self.vtoy(self.mdownloc[4]))

			if v0>=self.rng[0] and v1<=self.rng[1] : self.setValue(v0,v1)

#			print self.value
			
	def mouseReleaseEvent(self,event):
		self.mdownloc=None
		
	def vtoy(self,value):
		"returns the y screen coordinate corresponding to value"
		h=self.size().height()
		return h-int((h-4)*(value-self.rng[0])/(self.rng[1]-self.rng[0])+3)
		
	def ytov(self,y):
		"returns the value corresponding to the y screen coordinate"
		h=self.size().height()
		return int((h-y+3)*(self.rng[1]-self.rng[0])/(h-4)+self.rng[0])


	def setRange(self,minv,maxv,quiet=False):
		if maxv<=minv : maxv=minv+.001
		self.rng=(float(minv),float(maxv))
		self.setValue(self.value[0],self.value[1],quiet)

	def getRange(self): return self.rng

	def setValue(self,v0,v1,quiet=False):
		v0=clamp(self.rng[0],v0,self.rng[1])
		v1=clamp(self.rng[0],v1,self.rng[1])
		if v1<=v0 : v1=v0+1
		self.value=(v0,v1)
		self.update()
		if not quiet : self.emit(QtCore.SIGNAL("valueChanged"), self.value)

	def getValue(self):
		return self.value
		
class EMSpinWidget(QtGui.QWidget):
	"""
	Widget for a unbounded spin box using left and right arrow keys. When the value is changed
	valueChanged(int) is emited
	@param value is the starting value
	@param coeff is controls the exponential growth rate when the arrow is held down
	@param maxarrowwidth is the size of the arrow buttons
	"""
	def __init__(self, value, coeff, rounding=2, maxarrowwidth=20, postivemode=False, wheelstep=1):
		QtGui.QWidget.__init__(self)
		self.value = value
		self.coeff = coeff
		self.powercoeff = 0.0
		self.maxarrowwidth = maxarrowwidth
		self.postivemode = postivemode
		self.wheelstep=wheelstep
		self.rounding = rounding
		
		shbox = QtGui.QHBoxLayout()
		self.lbutton = QtGui.QPushButton("",self)
		self.lbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(leftarrow)))
		self.lbutton.setAutoRepeat(True)
		self.lbutton.setAutoRepeatDelay(200)
		shbox.addWidget(self.lbutton)
		self.numbox = QtGui.QLineEdit(str(self.value), self)
		shbox.addWidget(self.numbox)
		self.lbutton.setMaximumHeight(self.numbox.height())
		self.lbutton.setMaximumWidth(self.maxarrowwidth)
		self.rbutton = QtGui.QPushButton("",self)
		self.rbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(rightarrow)))
		self.rbutton.setAutoRepeat(True)
		self.rbutton.setAutoRepeatDelay(200)
		self.rbutton.setMaximumHeight(self.numbox.height())
		self.rbutton.setMaximumWidth(self.maxarrowwidth)
		shbox.addWidget(self.rbutton)
		self.setLayout(shbox)
		
		QtCore.QObject.connect(self.numbox,QtCore.SIGNAL("editingFinished()"),self._on_editfinish)
		QtCore.QObject.connect(self.lbutton,QtCore.SIGNAL("clicked()"),self._on_clickleft)
		QtCore.QObject.connect(self.rbutton,QtCore.SIGNAL("clicked()"),self._on_clickright)
		QtCore.QObject.connect(self.lbutton,QtCore.SIGNAL("released()"),self._on_unclickleft)
		QtCore.QObject.connect(self.rbutton,QtCore.SIGNAL("released()"),self._on_unclickright)
	
	def setValue(self, value, quiet=1):
		self.value = value
		self.numbox.setText(str(round(self.value, self.rounding)))
		if not quiet: self.emit(QtCore.SIGNAL("valueChanged(int)"),self.value)
	
	def getValue(self):
		return self.value
	
	def setRounding(self, rounding):
		self.rounding = rounding
		
	def setIncrement(self, incr):
		self.coeff = incr
		
	def wheelEvent(self, event):
		if event.delta() > 0:
			self.setValue(self.value+self.wheelstep)
		else:
			value = self.value-self.wheelstep
			if self.postivemode:
				if value > 0:
					self.setValue(value)
			else:
				self.setValue(value)
		self.emit(QtCore.SIGNAL("valueChanged(int)"),self.value)
		
	def _on_clickleft(self):
		self.value = self.value - self.coeff*(2**self.powercoeff)
		self.numbox.setText(str(round(self.value, self.rounding)))
		self.powercoeff += 0.1
		self.emit(QtCore.SIGNAL("valueChanged(int)"),self.value)
		
	def _on_clickright(self):
		self.value = self.value + self.coeff*(2**self.powercoeff)
		self.numbox.setText(str(round(self.value,self.rounding)))
		self.powercoeff += 0.1
		self.emit(QtCore.SIGNAL("valueChanged(int)"),self.value)
	
	def _on_unclickleft(self):
		if not self.lbutton.isDown():
			self.powercoeff = 0.0
			
	def _on_unclickright(self):
		if not self.rbutton.isDown():
			self.powercoeff = 0.0
			
	def _on_editfinish(self):
		try:
			oldvalue = self.value
			self.value = float(self.numbox.text())
			if self.value != oldvalue: self.emit(QtCore.SIGNAL("valueChanged(int)"),float(self.numbox.text()))
			
		except ValueError:
			self.numbox.setText(str(self.value))
			print "ERROR!!! You must enter a number"
		
class EMQTColorWidget(QtGui.QWidget):
	"""
	A widget displaying a color box that is used to control colors
	multiple boxes can be implemented
	"""
	def __init__(self, parent=None, red=255, green=255, blue=255, width=30, height=30):
		QtGui.QWidget.__init__(self, parent)
		self.width = width
		self.height = height
		self.color = QtGui.QColor(red,green,blue)
		self.setMinimumHeight(self.width+2)
		self.setMinimumWidth(self.height+2)
		self.setAcceptDrops(True)
	
	def setColor(self, color):
		self.color = color
		self.update()
		
	def getColor(self):
		return self.color
		
	def paintEvent(self, e):
		qp = QtGui.QPainter()
		qp.begin(self)
		self._draw_widget(qp)
		qp.end()
		
	def _draw_widget(self, qp):
		pen = QtGui.QPen(QtGui.QColor(0, 0, 0))
		pen.setWidth(2)
		qp.setPen(pen)
		qp.setBrush(self.color)
		qp.drawRect(0, 0, self.width, self.height)
	
	def dragEnterEvent(self, e):
		e.accept()
		
	def dropEvent(self, e):
		self.color = QtGui.QColor(e.mimeData().colorData())
		self.update()
		self.emit(QtCore.SIGNAL("newcolor(QColor)"), self.color)

	def mouseMoveEvent(self, e):

		if e.buttons() != QtCore.Qt.LeftButton:
			return
		self._dragdrop(e)

		
	def _dragdrop(self, e):
		mimeData = QtCore.QMimeData()
		mimeData.setColorData(self.color)
		drag = QtGui.QDrag(self)
		drag.setMimeData(mimeData)
		drag.setHotSpot(e.pos() - self.rect().topLeft())
		dropAction = drag.start(QtCore.Qt.MoveAction)
		
	def mousePressEvent(self, event):
		if event.buttons() != QtCore.Qt.RightButton:
			self.inicolor = self.color
			self.colrodialog = EMQtColorDialog(self.color)
			QtCore.QObject.connect(self.colrodialog,QtCore.SIGNAL("currentColorChanged(const QColor &)"),self._on_colorchange)
			QtCore.QObject.connect(self.colrodialog,QtCore.SIGNAL("colorSelected(const QColor &)"),self._on_colorselect)
			QtCore.QObject.connect(self.colrodialog,QtCore.SIGNAL("canceled()"),self._on_cancel)
			QtCore.QObject.connect(self.colrodialog,QtCore.SIGNAL("newconnection()"), self._on_additonal_connect)
		else:
			self._dragdrop(event)
			
	def _on_colorchange(self, color):
		if color.isValid():
			self.color = color
			self.update()
			self.emit(QtCore.SIGNAL("newcolor(QColor)"), self.color)
			
	def _on_colorselect(self, color):
		if color.isValid():
			self.color = color
			self.update()
			self.emit(QtCore.SIGNAL("newcolor(QColor)"), self.color)
			
	def _on_cancel(self):
		self.color = self.inicolor
		self.update()
		self.emit(QtCore.SIGNAL("newcolor(QColor)"), self.color)
		
	def _on_additonal_connect(self):
		QtCore.QObject.disconnect(self.colrodialog,QtCore.SIGNAL("currentColorChanged(const QColor &)"),self._on_colorchange)
		QtCore.QObject.disconnect(self.colrodialog,QtCore.SIGNAL("colorSelected(const QColor &)"),self._on_colorselect)
		QtCore.QObject.disconnect(self.colrodialog,QtCore.SIGNAL("canceled()"),self._on_cancel)
		QtCore.QObject.disconnect(self.colrodialog,QtCore.SIGNAL("newconnection()"), self._on_additonal_connect)

def singleton(cls):
	instances = {}
	def getinstance(inicolor):
		if cls not in instances:
			instances[cls] = cls(inicolor)
		else:
			instances[cls].emit(QtCore.SIGNAL("newconnection()"))
			instances[cls].setCurrentColor(inicolor)
			if instances[cls].hidden:
				instances[cls].show()
		return instances[cls]
	return getinstance
    
@singleton
class EMQtColorDialog(QtGui.QColorDialog):
	"""
	The Is to create a non-modal color dialog. Only one color dialog is allowed at once, so I use the singltion pattern
	"""
	def __init__(self, inicolor):
		QtGui.QColorDialog.__init__(self, inicolor)
		self.hidden = False
		self.setWindowFlags(QtCore.Qt.Window)
		self.show()
	
	def hideEvent(self, e):
		QtGui.QColorDialog.hideEvent(self, e)
		self.emit(QtCore.SIGNAL("canceled()"))
		self.hidden = True

class EMLightControls(QtOpenGL.QGLWidget):
	"""
	Widget to set the postion of a light in 3D. When the light position is moved  it emits
	its position in spherical corridinants (theta and phi) as the signal: lightPositionMoved
	Its position can be set via: setAngularPosition
	@param light, the glLight the this widget uses
	"""
	def __init__(self, light, parent=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
		self.light = light
		self.theta = 0
		self.phi = 0
		self.lightposition = [0.0, 0.0, 1.0, 0.0]
		self.ambient = [1.0,1.0,1.0,1.0]
		self.diffuse = [1.0,1.0,1.0,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
	
	def initializeGL(self):
		glClearColor(0.0, 0.0, 0.0, 0.0)		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_LIGHTING)
		glEnable(self.light)
		self.setLightColor()
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		# Draw the flashlight
		self.drawFlashLight()
		# Draw the sphere
		glMaterialfv(GL_FRONT, GL_AMBIENT, [0.3,0.3,0.3,1.0])
		glMaterialfv(GL_FRONT, GL_DIFFUSE, [1.0,1.0,1.0,1.0])
		glMaterialfv(GL_FRONT, GL_SPECULAR, [1.0,1.0,1.0,1.0])
		glLightfv(self.light, GL_POSITION, self.lightposition)
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		quad = GLU.gluNewQuadric()
		GLU.gluQuadricDrawStyle(quad, GLU.GLU_FILL)
		GLU.gluQuadricNormals(quad, GLU.GLU_SMOOTH)
		GLU.gluSphere(quad, 3.0, 30, 30)
		glFlush()
		
		
	def drawFlashLight(self):
		glPushMatrix()
		glRotate(-self.phi, 1.0, 0.0, 0.0)
		glRotate(-self.theta, 0.0, 1.0, 0.0)
		glTranslate(0,0,4)
		glMaterialfv(GL_FRONT, GL_AMBIENT, [0.0,0.0,0.0,1.0])
		glMaterialfv(GL_FRONT, GL_DIFFUSE, [0.0,0.0,0.0,1.0])
		glMaterialfv(GL_FRONT, GL_SPECULAR, [0.0,0.0,0.0,1.0])
		glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,1.0,0.0,1.0])
		flashlight = GLU.gluNewQuadric()
		GLU.gluQuadricDrawStyle(flashlight, GLU.GLU_FILL)
		GLU.gluCylinder(flashlight, 0.3, 0.3, 1.5, 30,30)
		GLU.gluCylinder(flashlight, 0.5, 0.3, 0.4, 30,30)
		glPushMatrix()
		glTranslate(0,0,1.5)
		GLU.gluDisk(flashlight, 0.0, 0.3, 30,30)
		glPopMatrix()
		glMaterialfv(GL_FRONT, GL_EMISSION, self.diffuse)
		GLU.gluDisk(flashlight, 0.0, 0.5, 30,30)
		glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,0.0,0.0,1.0])
		glPopMatrix()
		
	def resizeGL(self, width, height):
		glViewport(0,0,width,height)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		GLU.gluPerspective(60.0, (float(width)/float(height)), 1.0, 100.0)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glTranslate(0,0,-10.0)
	
	def mousePressEvent(self, event):
		self.init_x = event.x()
		self.init_y = event.y()
		
	def mouseMoveEvent(self, event):
		""" When the mouse moves, move the light"""
		self.theta -= event.x() - self.init_x
		self.phi -= event.y() - self.init_y
		self.theta = self.theta % 360
		self.phi = self.phi % 360

		self.init_x = event.x()
		self.init_y = event.y()
		self.update()
		self.emit(QtCore.SIGNAL("lightPositionMoved"), [self.theta, self.phi])
		
	def setPosition(self, x, y, z, w, quiet=False):
		self.lightposition = [x, y, z, w]
		if not quiet: self.update()
		
	def getPosition(self):
		return self.lightposition
		
	def getAngularPosition(self):
		return [self.theta, self.phi]
	
	def setAngularPosition(self, theta, phi):
		""" The position is in spherical coords"""
		self.theta = theta
		self.phi = phi
		z = math.sin(math.radians(self.theta + 90))*math.cos(math.radians(self.phi))
		y = math.sin(math.radians(self.theta + 90))*math.sin(math.radians(self.phi))
		x = math.cos(math.radians(self.theta + 90))
		self.setPosition(x, y, z, 0.0, quiet=True)
		self.update()
	
	def setLightColor(self):
		glLightfv(self.light, GL_AMBIENT, self.ambient)
		glLightfv(self.light, GL_DIFFUSE, self.diffuse)
		glLightfv(self.light, GL_SPECULAR, self.specular)
		
	def setAmbient(self, r, g, b):
		self.ambient = [r, g, b, 1.0]
		self.setLightColor()
		self.update()
	
	def setDiffuse(self, r, g, b):
		self.diffuse = [r, g, b, 1.0]
		self.setLightColor()
		self.update()
		
	def setSpecular(self, r, g, b):
		self.specular = [r, g, b, 1.0]
		self.setLightColor()
		self.update()
		
class CameraControls(QtOpenGL.QGLWidget):
	"""
	Widget to set the camera position. When clipping planes are moved a farMoved(float) or nearMoved(float) signal is emmited.
	This widget is an observer of the Camera object (used in the Scenegraph). To update call updateWidget
	@param scenegraph, the scenegraph this widget oberrves
	"""
	def __init__(self, parent=None, scenegraph=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
		self.scenegraph = weakref.ref(scenegraph)
		self.texture = glGenTextures(1)
		
	def __del__(self):
		glDeleteTextures(self.texture)
		
	def initializeGL(self):
		glClearColor(0.0, 0.0, 0.0, 0.0)		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)

	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		self._drawZslice()
		self._drawViewingVolume()
		glFlush()
		
	def resizeGL(self, width, height):
		self.width = width
		self.height = height
		#self.setMaximumWidth(width)
		#self.setMaximumHeight(height)
		glViewport(0,0,width,height)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(-width/2, width/2, -height/2, height/2, -100.0, 100.0)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
	def mousePressEvent(self, event):
		self.init_x = event.x()
		
	def mouseMoveEvent(self, event):
		""" Move the clipping planes"""
		self.movement = float(event.x() - self.init_x)*self.scale
		if math.fabs(event.x()-(self.near_clipping + self.width/2)) > math.fabs(event.x()-(self.far_clipping + self.width/2)):
			self.emit(QtCore.SIGNAL("farMoved(float)"), self.movement)
		else:
			self.emit(QtCore.SIGNAL("nearMoved(float)"), self.movement)
		self.init_x = event.x()
		
	def _drawViewingVolume(self):
		# These current scheme draws what looks like a frustrum, which is a bit of a lie, b/c usually a orthographic projection is used. But is it hard to depict such
		#situation using a camera analogy
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glColor3f(0.7, 0.7, 0.0)
		sixtydegrees = math.sin(math.radians(self.scenegraph().camera.getFovy()))
		size = min(self.width, self.height)
		self.scale = float(self.scenegraph().camera.getWidth())/float(size)*self.scenegraph().camera.getViewPortWidthScaling()
		origin = 0.0
		#print self.scenegraph().camera.getClipNear()
		self.near_clipping = origin + (self.scenegraph().camera.getClipNear() + self.scenegraph().camera.getZclip())/self.scale
		self.far_clipping = origin + (self.scenegraph().camera.getClipFar() + self.scenegraph().camera.getZclip())/self.scale
		glBegin(GL_LINES)
		glVertex(self.near_clipping, -self.height/2.2, 0)
		glVertex(self.near_clipping, self.height/2.2, 0)
		glVertex(self.far_clipping, -self.height/2.2, 0)
		glVertex(self.far_clipping, self.height/2.2, 0)
		glEnd()
		
	def _drawZslice(self):
		""" Draw the Z slice in a box as the texture"""
		#glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.pixels = self.scenegraph().pixels
		glBindTexture(GL_TEXTURE_2D, self.texture)
		
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, self.pixels[2], self.pixels[3], 0, GL_RGBA, GL_UNSIGNED_BYTE, self.pixels[4])
		
		
		glEnable(GL_TEXTURE_2D)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBindTexture(GL_TEXTURE_2D, self.texture)
		aspectratio = float(self.scenegraph().camera.getHeight())/float(self.scenegraph().camera.getWidth())

		glBegin(GL_QUADS)
		glTexCoord2f(0.0,0.0)
		size = min(self.width, self.height)
		glVertex(-size/2,-aspectratio*size/2,-1)
		glTexCoord2f(1.0,0.0)
		glVertex(size/2,-aspectratio*size/2,-1)
		glTexCoord2f(1.0,1.0)
		glVertex(size/2,aspectratio*size/2,-1)
		glTexCoord2f(0.0,1.0)
		glVertex(-size/2,aspectratio*size/2,-1)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		
	def updateWidget(self):
		self.update()
		
class EMANToolButton(QtGui.QToolButton):
	"""
	This widget is a toolbutton that is toogleable. Only one button can be down at any given point.
	The only weakness is only one of these button groups can be used at any given time. It is possible 
	to fix this, but when I need such functionality
	"""
	toolpanellist = []
	def __init__(self):
		QtGui.QToolButton.__init__(self)
		self.setMinimumWidth(30)
		self.setMinimumHeight(30)
		self.weakrefself = weakref.ref(self)
		EMANToolButton.toolpanellist.append(self.weakrefself)
		
	def __del__(self):
		EMANToolButton.toolpanellist.remove(self.weakrefself)
	
	def setSelfAsUnique(self):
		for tool in EMANToolButton.toolpanellist:
			if tool() != self:
				tool().setDown(False)
	
	def setDown(self, state, quiet=True):
		QtGui.QToolButton.setDown(self, state)
		if state: self.setSelfAsUnique()
		if not quiet: self.emit(QtCore.SIGNAL("clicked(int)"), self.isDown())
		
	def mousePressEvent(self, event):
		# Toggle the button on and off
		if not self.isDown():
			self.setSelfAsUnique()
			self.setDown(True, quiet=False)
		
	def mouseReleaseEvent(self, event):
		pass

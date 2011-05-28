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
from PyQt4.QtCore import Qt

def clamp(x0,val,x1):
	return int(max(min(val,x1),x0))

class ValSlider(QtGui.QWidget):
	"""The valslider class represents a connected text widget and horizontal slider.
	showenable - if -1, no enable box shown, if 0, shown unchecked, if 1 shown and checked
	setValue(float) - to programatically change the value
	emit valueChanged(float)
	"""
	def __init__(self, parent=None, rng=None, label=None, value=0,labelwidth=30,showenable=-1):
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
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.enabletog)
		
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
		if showenable>=0 : self.enabletog(showenable)
		
	#def __del__(self):
		#print self.getLabel(), " freed"
#		QtGui.QWidget.__del__(self)

	def enabletog(self,ena):
		self.slider.setEnabled(ena)
		self.text.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enablechange"),ena) 

	def setRange(self,minv,maxv):
		if maxv<=minv : maxv=minv+.001
		self.rng=[float(minv),float(maxv)]
		self.updates()
		#self.slider.setRange(*self.rng)
		
	def getRange(self): return self.rng

	def setValue(self,val,quiet=0):
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
		self.text.setText(str(self.value)[:self.text.width()/10-1])
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
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.enabletog)
			self.enabletog(showenable)
			
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
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		
		self.updateboth()

	def enabletog(self,ena):
		self.text.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enablechange"),ena) 
		
	def setRange(self,minv,maxv):
		if maxv<=minv : maxv=minv+.001
		self.rng=[float(minv),float(maxv)]
		self.updates()
		#self.slider.setRange(*self.rng)
		
	def getRange(self): return self.rng

	def setValue(self,val,quiet=0):
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
		self.text.setText(str(self.value)[:self.text.width()/10-1])
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
		self.value=value
		self.ignore=0
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		if showenable>=0 :
			self.enablebox=QtGui.QCheckBox(self)
			self.enablebox.setChecked(showenable)
			self.hboxlayout.addWidget(self.enablebox)
			QtCore.QObject.connect(self.enablebox, QtCore.SIGNAL("toggled(bool)"), self.enabletog)
			self.enabletog(showenable)
			
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
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		
		self.updateboth()

	def enabletog(self,ena):
		self.text.setEnabled(ena)
		self.emit(QtCore.SIGNAL("enablechange"),ena) 
		
	def setValue(self,val,quiet=0):
		if self.value==val : return
		self.value=val
		self.updateboth()
		if not quiet : self.emit(QtCore.SIGNAL("valueChanged"),self.value)
	
	def getValue(self):
		return self.value
	
		
	def textChange(self):
		if self.ignore : return
		self.value=self.text.text()
		self.emit(QtCore.SIGNAL("valueChanged"),self.value)
		self.emit(QtCore.SIGNAL("textChanged"),self.value)
				
	def setLabel(self,label):
		self.label.setText(label)
	
	def getLabel(self):
		return str(self.label.text())
		
	def updatet(self):
		self.ignore=1
#		if self.validate!=None : self.validate()
		self.ignore=0
		
	def updateboth(self):
		self.updatet()


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
	
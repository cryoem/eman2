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
from EMAN2 import *
from valslider import ValSlider
from emanimationutil import Animator
import weakref
import copy
import sys
import math

class EMTransformPanel:
	def __init__(self,target,parent):
		self.target = weakref.ref(target)
		self.parent = weakref.ref(parent)
		
		self.label_src = QtGui.QLabel(parent)
		self.label_src.setText('Rotation Convention')
		
		self.src = QtGui.QComboBox(parent)
		self.load_src_options(self.src)
		
		self.x_label = QtGui.QLabel()
		self.x_label.setText('x')
		
		self.x_trans = QtGui.QDoubleSpinBox(parent)
		self.x_trans.setMinimum(-10000)
		self.x_trans.setMaximum(10000)
		self.x_trans.setValue(0.0)
	
		self.y_label = QtGui.QLabel()
		self.y_label.setText('y')
		
		self.y_trans = QtGui.QDoubleSpinBox(parent)
		self.y_trans.setMinimum(-10000)
		self.y_trans.setMaximum(10000)
		self.y_trans.setValue(0.0)
		
		self.z_label = QtGui.QLabel()
		self.z_label.setText('z')
		
		self.z_trans = QtGui.QDoubleSpinBox(parent)
		self.z_trans.setMinimum(-10000)
		self.z_trans.setMaximum(10000)
		self.z_trans.setValue(0.0)
		
		self.az = ValSlider(parent,(-360.0,360.0),"az",-1)
		self.az.setObjectName("az")
		self.az.setValue(0.0)
		
		self.alt = ValSlider(parent,(-180.0,180.0),"alt",-1)
		self.alt.setObjectName("alt")
		self.alt.setValue(0.0)
		
		self.phi = ValSlider(parent,(-360.0,360.0),"phi",-1)
		self.phi.setObjectName("phi")
		self.phi.setValue(0.0)
		
		self.scale = ValSlider(parent,(0.01,30.0),"Zoom:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		
		self.n3_showing = False
		
		self.current_src = "eman"
		
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), self.target().set_scale)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), self.target().set_cam_x)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), self.target().set_cam_y)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), self.target().set_cam_z)
		
		
	def set_defaults(self):
		self.x_trans.setValue(0.0)
		self.y_trans.setValue(0.0)
		self.z_trans.setValue(0.0)
		self.scale.setValue(1.0)
		self.az.setValue(0.0)
		self.alt.setValue(0.0)
		self.phi.setValue(0.0)
		
	def slider_rotate(self):
		self.target().load_rotation(self.get_current_rotation())
		
	def get_current_rotation(self):
		convention = self.src.currentText()
		rot = {}
		if ( self.current_src == "spin" ):
			rot[self.az.getLabel()] = self.az.getValue()
			
			n1 = self.alt.getValue()
			n2 = self.phi.getValue()
			n3 = self.n3.getValue()
			
			norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			n1 /= norm
			n2 /= norm
			n3 /= norm
			
			rot[self.alt.getLabel()] = n1
			rot[self.phi.getLabel()] = n2
			rot[self.n3.getLabel()] = n3
			
		else:
			rot[self.az.getLabel()] = self.az.getValue()
			rot[self.alt.getLabel()] = self.alt.getValue()
			rot[self.phi.getLabel()] = self.phi.getValue()
		
		rot["type"] = self.current_src
		
		return Transform(rot)
	
	def addWidgets(self,target):
		
		target.addWidget(self.scale)
		self.hbl_trans = QtGui.QHBoxLayout()
		self.hbl_trans.setMargin(0)
		self.hbl_trans.setSpacing(6)
		self.hbl_trans.setObjectName("Trans")
		self.hbl_trans.addWidget(self.x_label)
		self.hbl_trans.addWidget(self.x_trans)
		self.hbl_trans.addWidget(self.y_label)
		self.hbl_trans.addWidget(self.y_trans)
		self.hbl_trans.addWidget(self.z_label)
		self.hbl_trans.addWidget(self.z_trans)
		
		target.addLayout(self.hbl_trans)
		
		self.hbl_src = QtGui.QHBoxLayout()
		self.hbl_src.setMargin(0)
		self.hbl_src.setSpacing(6)
		self.hbl_src.setObjectName("hbl")
		self.hbl_src.addWidget(self.label_src)
		self.hbl_src.addWidget(self.src)
		
		
		target.addLayout(self.hbl_src)
		target.addWidget(self.az)
		target.addWidget(self.alt)
		target.addWidget(self.phi)
	
	def set_src(self, val):
		t3d = self.get_current_rotation()
		
		if (self.n3_showing) :
			self.parent().get_transform_layout().removeWidget(self.n3)
			self.n3.deleteLater()
			self.n3_showing = False
			self.az.setRange(-360,360)
			self.alt.setRange(-180,180)
			self.phi.setRange(-360,660)
		
		if ( self.src_map[str(val)] == "spider" ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('psi')
		elif ( self.src_map[str(val)] == "eman" ):
			self.az.setLabel('az')
			self.alt.setLabel('alt')
			self.phi.setLabel('phi')
		elif ( self.src_map[str(val)] == "imagic"):
			self.az.setLabel('alpha')
			self.alt.setLabel('beta')
			self.phi.setLabel('gamma')
		elif ( self.src_map[str(val)] == "xyz"):
			self.az.setLabel('xtilt')
			self.alt.setLabel('ytilt')
			self.phi.setLabel('ztilt')
		elif ( self.src_map[str(val)] == "mrc" ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('omega')
		elif ( self.src_map[str(val)] == "spin" ):
			self.az.setLabel('omega')
			self.alt.setRange(-1,1)
			self.phi.setRange(-1,1)
			
			self.alt.setLabel('n1')
			self.phi.setLabel('n2')
			
			self.n3 = ValSlider(self.parent(),(-360.0,360.0),"n3",-1)
			self.n3.setRange(-1,1)
			self.n3.setObjectName("n3")
			self.parent().get_transform_layout().addWidget(self.n3)
			QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
			self.n3_showing = True
		
		self.current_src = self.src_map[str(val)]
		self.update_rotations(t3d)
	
	def load_src_options(self,widgit):
		self.load_src()
		for i in self.src_strings:
			widgit.addItem(i)
			
	def load_src(self):
		# supported_rot_conventions
		src_flags = []
		src_flags.append("eman")
		src_flags.append("spider")
		src_flags.append("imagic")
		src_flags.append("mrc")
		src_flags.append("spin")
		src_flags.append("xyz")
		
		self.src_strings = []
		self.src_map = {}
		for i in src_flags:
			self.src_strings.append(str(i))
			self.src_map[str(i)] = i
			
	def update_rotations(self,t3d):
		rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		convention = self.src.currentText()
		if ( self.src_map[str(convention)] == "spin" ):
			self.n3.setValue(rot[self.n3.getLabel()],True)
		
		self.az.setValue(rot[self.az.getLabel()],True)
		self.alt.setValue(rot[self.alt.getLabel()],True)
		self.phi.setValue(rot[self.phi.getLabel()],True)
		
	def set_scale(self,newscale):
		self.scale.setValue(newscale)
		
	def set_xy_trans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
		
	def set_xyz_trans(self, x, y,z):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
		self.z_trans.setValue(z)
	
import weakref

class EMParentWin(QtGui.QWidget,Animator):
	"""
	This class adds a status bar with a size grip to QGLWidgets on Mac OS X, 
	to provide a visual cue that the window can be resized. This is accomplished
	by placing a QGLWidget and a QStatusBar onto a QWidget. The EMGLWidget 
	class hides the fact that this workaround is required.
	"""
	
	def __init__(self,child,enable_timer=False):
		"""
		@param child: the central GL display widget
		@type child: EMGLWidget
		@param enable_timer: not used... historical purposes???
		"""
		#TODO: figure out why the enable_timer parameter isn't being used
		QtGui.QWidget.__init__(self,None)
		Animator.__init__(self)

		self.child = weakref.ref(child) # Either EMGLWidget.parent_widget or EMParentWin.child must be a weakref for garbage collection purposes.
		self.resize(child.width(),child.height())
		self.setMaximumSize(8000,8000)

		self.hbl = QtGui.QVBoxLayout()
		
		self.hbl.setSpacing(0)
		self.hbl.addWidget(self.child(),100)
		if get_platform() == "Darwin": # because OpenGL widgets in Qt don't leave room in the bottom right hand corner for the resize tool
			self.status = QtGui.QStatusBar()
			self.status.setSizeGripEnabled(True)
			self.hbl.addWidget(self.status,0)
			self.margin = 0
		else:
			self.margin = 5
		self.hbl.setMargin(self.margin)
		self.setLayout(self.hbl)
		
		
	def __del__(self):
		pass
	
	def get_margin(self):
		return 50
	

	def closeEvent(self, e):
		try:
			self.child().closeEvent(e)
			#self.child.inspector.close()
		except: pass
		QtGui.QWidget.closeEvent(self,e)
		
	def resizeEvent(self,event):
		self.child().resizeEvent(event)
	
	def get_qt_widget(self):
		return self.child()
	
	def update(self):
		self.child().updateGL()
		
	def updateGL(self):
		self.child().updateGL()
	
	def keyPressEvent(self,event):
		self.child().keyPressEvent(event)

	#def mousePressEvent(self, event):
		#print "mouse",event
		#self.child().mousePressEvent(event)


	#def width(self):
		##print "asked for width!"
		#return self.child.width()
	
	#def height(self):
		##print "asked for height!"
		#return self.child.height()
		
	def initGL(self):
		self.child().glInit()
	
class ImgHistogram(QtGui.QWidget):
	""" A small fixed-size histogram widget"""
	def __init__(self,parent):
		QtGui.QWidget.__init__(self,parent)
		
		self.brush=QtGui.QBrush(Qt.black)
		self.font=QtGui.QFont("Helvetica", 12);
		self.probe=None
		self.histdata=None
		self.setMinimumSize(QtCore.QSize(258,128))
		self.volume = False
		
		self.upposition = 0
		self.upcolor = 0
		self.uppresent = 0
		self.guiparent = parent
		
		self.mpressed = None

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
	
	def setProbe(self, value):
		self.volume = False
		
		if (self.maxden-self.minden) != 0:
			x = int(255.0*(value - self.minden)/(self.maxden-self.minden))
			if x > 255: x = 255
			if x < 0: x = 0
			self.threshold = value
			self.probe=(x,self.histdata[x])
		else:
			self.threshold = value
		
		value=7
		self.update()
	def setColor(self,value):

		self.upcolor[self.uppresent]=value
		self.update()

	def setDynamicProbe(self, position, color, present, value):
		self.volume=True
		
		self.guiparent.level.setText("Level: %1.3f"%value)

		self.upposition = position
		self.upcolor = color
		self.uppresent = present
		
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
		
		if self.histdata==None : return
		p=QtGui.QPainter()
		p.begin(self)
		p.setBackground(QtGui.QColor(16,16,16))
		p.eraseRect(0,0,self.width(),self.height())
		p.setPen(Qt.darkGray)
		for i,j in enumerate(self.histdata):
			p.drawLine(i,127,i,127-j*126/self.norm)

		# If the user has dragged, we need to show a value

		if (self.probe) and (self.volume == False):
			p.setPen(Qt.blue)
			p.drawLine(self.probe[0]+1,0,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setPen(Qt.red)
			p.drawLine(self.probe[0]+1,127,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setFont(self.font)
			#p.drawText(200,20,"x=%d"%(self.probe[0]))
			#p.drawText(200,36,"%1.2f"%self.threshold)
			#p.drawText(200,52,"y=%d"%(self.probe[1]))
			#p.drawText(200,68,"%1.2f%%"%(100.0*float(self.probe[1])/self.total))
			p.drawText(150,20,"thres=%1.2f"%self.threshold)
			p.drawText(150,36,"p(t)=%1.2f%%"%(100.0*float(self.probe[1])/self.total))
			
		if  (self.upposition !=0)and (self.volume):
			#self.uppresent=0 
			#self.upposition=[[185,121],[185,18],[255,18]] 
			
			p.setPen(QtGui.QColor(0,0,0))
			p.setBrush(QtGui.QBrush(Qt.SolidPattern))
			p.setBrush(QtGui.QColor(self.upcolor[0][0],self.upcolor[0][1],self.upcolor[0][2]))
			p.drawRect(self.upposition[0][0]-2.5,self.upposition[0][1]-2.5,5,5)
			p.setBrush(QtGui.QBrush(Qt.NoBrush))
			
			p.setPen(QtGui.QColor(0,0,0))
			p.setBrush(QtGui.QBrush(Qt.SolidPattern))
			p.setBrush(QtGui.QColor(self.upcolor[1][0],self.upcolor[1][1],self.upcolor[1][2]))
			p.drawRect(self.upposition[1][0]-2.5,self.upposition[1][1]-2.5,5,5)
			p.setBrush(QtGui.QBrush(Qt.NoBrush))
			
			p.setPen(QtGui.QColor(0,0,0))
			p.setBrush(QtGui.QBrush(Qt.SolidPattern))
			p.setBrush(QtGui.QColor(self.upcolor[2][0],self.upcolor[2][1],self.upcolor[2][2]))
			p.drawRect(self.upposition[2][0]-2.5,self.upposition[2][1]-2.5,5,5)
			
			p.setBrush(QtGui.QBrush(Qt.NoBrush))
			
			p.setPen(Qt.yellow)
			presentx=self.upposition[self.uppresent][0]
			presenty=self.upposition[self.uppresent][1]
			self.uppositioned=sorted(self.upposition)
			
			p.drawLine(self.uppositioned[0][0],self.uppositioned[0][1],self.uppositioned[1][0],self.uppositioned[1][1])
			p.drawLine(self.uppositioned[1][0],self.uppositioned[1][1],self.uppositioned[2][0],self.uppositioned[2][1])
			
			
			self.upcolorpreserve = copy.deepcopy(self.upcolor) 
			for i in range(3):
				if (self.uppositioned[i][0]==presentx) and (self.uppositioned[i][1]==presenty):
					self.uppresent = i
			
			
			for i in range(3):
				if self.upposition[0]==self.uppositioned[i]:
					self.upcolor[0]=self.upcolorpreserve[i]
					
			for i in range(3):
				if self.upposition[1]==self.uppositioned[i]:
					self.upcolor[1]=self.upcolorpreserve[i]
					
			for i in range(3):
				if self.upposition[2]==self.uppositioned[i]:
					self.upcolor[2]=self.upcolorpreserve[i]		
				
			self.upposition.sort()
			 
			presentthreshold = (self.upposition[self.uppresent][0]-3)/252.0*(self.maxden-self.minden)+self.minden
			self.guiparent.level.setText("Level: %1.3f"%presentthreshold)
			#self.guiparent.cappingcolor.setRgb(self.upcolor[self.uppresent])
			
			#p.setPen(QtGui.QColor(self.upcolor[1][0],self.upcolor[1][1],self.upcolor[1][2]))
			#p.setPen(QtGui.QColor(self.upcolor[2][0],self.upcolor[2][1],self.upcolor[2][2]))

		p.setPen(Qt.black)
		p.drawRect(0,0,258,128)
		p.end()
	
	def mousePressEvent(self, event):

		if (event.button()==Qt.LeftButton):
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.threshold = (self.probe[0]/255.0*(self.maxden-self.minden)+self.minden)
			if (self.volume == False):
				self.emit(QtCore.SIGNAL("thresholdChanged(float)"), self.threshold)
			self.update()
		
		if (event.button()==Qt.LeftButton) and (self.volume):
			if ((event.x()-self.upposition[0][0])**2+(event.y()-self.upposition[0][1])**2)<49:
				self.uppresent =0
				self.mpressed = True
			elif ((event.x()-self.upposition[1][0])**2+(event.y()-self.upposition[1][1])**2)<49:
				self.uppresent =1
				self.mpressed = True
			elif ((event.x()-self.upposition[2][0])**2+(event.y()-self.upposition[2][1])**2)<49:
				self.uppresent =2
				self.mpressed = True
			else: 
				self.mpressed = False
				
			if (self.mpressed) :
				self.upposition[self.uppresent][0]=event.x()
				self.upposition[self.uppresent][1]=event.y()
				if self.upposition[self.uppresent][0] > 255: self.upposition[self.uppresent][0] =255 
				if self.upposition[self.uppresent][0] < 3: self.upposition[self.uppresent][0] =3
				if self.upposition[self.uppresent][1] > 125: self.upposition[self.uppresent][1] =125 
				if self.upposition[self.uppresent][1] < 3: self.upposition[self.uppresent][1] =3
				
			self.update()
			
	def mouseMoveEvent(self, event):
		if (event.buttons()&Qt.LeftButton) and (self.volume==False):
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.threshold = (self.probe[0]/255.0*(self.maxden-self.minden)+self.minden)
			self.emit(QtCore.SIGNAL("thresholdChanged(float)"), self.threshold)
			self.update()
		
		if (Qt.LeftButton) and (self.volume):
			if self.mpressed:
				self.upposition[self.uppresent][0]=event.x()
				self.upposition[self.uppresent][1]=event.y()
				if self.upposition[self.uppresent][0] > 255: self.upposition[self.uppresent][0] =255 
				if self.upposition[self.uppresent][0] < 3: self.upposition[self.uppresent][0] =3
				if self.upposition[self.uppresent][1] > 125:self.upposition[self.uppresent][1] =125 
				if self.upposition[self.uppresent][1] < 3: self.upposition[self.uppresent][1] =3
				
				self.update()
				
	def mouseReleaseEvent(self, event):
		#self.probe=None
		self.mpressed = False
		pass
		
class EMMetaDataTable(object):
	"""This is basically a factory class that will return an instance of QtWidget
	"""
	def __new__(cls,parent,metadata):
		'''
		metadata should be a dict
		'''
		if not isinstance(metadata,dict): raise
		
		left = [str(k) for k in metadata.keys()]
		right = [str(v) for v in metadata.values()]
		
		from emform import EMParamTable, ParamDef,EMFormWidget
		
		params = []
		a = EMParamTable(name="Metadata",desc_short="",desc_long="Meta data associated with this image")
		pleft = ParamDef(name="key",vartype="stringlist",desc_short="Key",desc_long="The key of the metadata value object",property=None,choices=left)
		pright = ParamDef(name="value",vartype="stringlist",desc_short="Value",desc_long="The value of the metadata object as a string",property=None,choices=right)

		
		a.append(pleft)
		a.append(pright)
		params.append(a)
		
		form = EMFormWidget(parent,params,disable_ok_cancel=True)
		return form
	
	
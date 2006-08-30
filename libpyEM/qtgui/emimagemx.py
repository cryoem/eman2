#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import Numeric
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary

class EMImageMX(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImageMX.allim[self]=0
		
		self.data=None
		self.datasize=(1,1)
		self.scale=1.0
		self.minden=0
		self.maxden=1.0
		self.invert=0
		self.fft=None
		self.mindeng=0
		self.maxdeng=1.0
		self.origin=(0,0)
		self.nperrow=8
		self.nshow=-1
		self.mousedrag=None
		self.nimg=0
		self.changec={}
		
		self.inspector=None
		if data: 
			self.setData(data)
			self.show()
			
			
	def setData(self,data):
		if not self.data and data:
			try:
				if len(data)<self.nperrow : w=len(data)*(data[0].get_xsize()+2)
				else : w=self.nperrow*(data[0].get_xsize()+2)
				if w>0 : self.resize(w,512)
			except: pass
		self.data=data
		if data==None or len(data)==0:
			self.updateGL()
			return
		
		self.nimg=len(data)
		
		self.minden=data[0].get_attr("mean")
		self.maxden=self.minden
		self.mindeng=self.minden
		self.maxdeng=self.minden
		
		for i in data:
			if i.get_zsize()!=1 :
				self.data=None
				self.updateGL()
				return
			mean=i.get_attr("mean")
			sigma=i.get_attr("sigma")
			m0=i.get_attr("minimum")
			m1=i.get_attr("maximum")
		
			self.minden=min(self.minden,max(m0,mean-3.0*sigma))
			self.maxden=max(self.maxden,min(m1,mean+3.0*sigma))
			self.mindeng=min(self.mindeng,max(m0,mean-5.0*sigma))
			self.maxdeng=max(self.maxdeng,min(m1,mean+5.0*sigma))
		
		self.showInspector()		# shows the correct inspector if already open
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setOrigin(self,x,y):
		"""Set the display origin within the image"""
		self.origin=(x,y)
		self.updateGL()
		
	def setScale(self,newscale):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		
#		yo=self.height()-self.origin[1]-1
		yo=self.origin[1]
#		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+yo)-self.height()/2)
		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(yo-self.height()/2)+self.height()/2)
#		print self.origin,newscale/self.scale,yo,self.height()/2+yo
		
		self.scale=newscale
		self.updateGL()
		
	def setDenMin(self,val):
		self.minden=val
		self.updateGL()
		
	def setDenMax(self,val):
		self.maxden=val
		self.updateGL()

	def setNPerRow(self,val):
		self.nperrow=val
		self.updateGL()
		
	def setNShow(self,val):
		self.nshow=val
		self.updateGL()

	def setInvert(self,val):
		if val: self.invert=1
		else : self.invert=0
		self.updateGL()

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		for i in self.data:
			self.changec[i]=i.get_attr("changecount")
		
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)
		
		GL.glPixelZoom(1.0,-1.0)
		n=len(self.data)
		x,y=-self.origin[0],self.height()-self.origin[1]-1
		hist=Numeric.zeros(256)
		for i in range(n):
			w=int(min(self.data[i].get_xsize()*self.scale,self.width()))
			h=int(min(self.data[i].get_ysize()*self.scale,self.height()))
			if x>0 and x<self.width() and y>0 and y<self.height() :
				a=self.data[i].render_amp8(0,0,w,h,(w-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,2)
				GL.glRasterPos(x,y)
				GL.glDrawPixels(w,h,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
				hist2=Numeric.fromstring(a[-1024:],'i')
				hist+=hist2
			elif x+w>0 and y-h<self.height() and x<self.width() and y>0:
				if x<0 : 
					x0=-x/self.scale
					x1=w+x
				else : 
					x0=0
					x1=w
				if y>self.height()-1 : y1=h-y+self.height()-1
				else : y1=h
				x0,x1,y1=int(x0),int(x1),int(y1)
				a=self.data[i].render_amp8(x0,0,x1,y1,(x1-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,2)
				if x<0 : xx=0
				else : xx=x
				if y>=self.height() : yy=self.height()-1
				else : yy=y
				GL.glRasterPos(xx,yy)
				GL.glDrawPixels(x1,y1,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
				hist2=Numeric.fromstring(a[-1024:],'i')
				hist+=hist2
				
			
			if (i+1)%self.nperrow==0 : 
				y-=h+2.0
				x=-self.origin[0]
			else: x+=w+2.0
		if self.inspector : self.inspector.setHist(hist,self.minden,self.maxden)
	
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		if not self.inspector : self.inspector=EMImageMxInspector2D(self)
		self.inspector.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
		self.inspector.show()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.RightButton:
			self.mousedrag=(event.x(),event.y())
	
	def mouseMoveEvent(self, event):
		if self.mousedrag:
			self.origin=(self.origin[0]+self.mousedrag[0]-event.x(),self.origin[1]-self.mousedrag[1]+event.y())
			self.mousedrag=(event.x(),event.y())
			self.update()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.RightButton:
			self.mousedrag=None

class EMImageMxInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vboxlayout")
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.vbl.addWidget(self.hist)

		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)
		
		self.mmeas = QtGui.QPushButton("App")
		self.mmeas.setCheckable(1)
		self.hbl2.addWidget(self.mmeas)

		self.mdel = QtGui.QPushButton("Del")
		self.mdel.setCheckable(1)
		self.hbl2.addWidget(self.mdel)

		self.mmove = QtGui.QPushButton("Move")
		self.mmove.setCheckable(1)
		self.hbl2.addWidget(self.mmove)

		self.bg=QtGui.QButtonGroup()
		self.bg.setExclusive(1)
		self.bg.addButton(self.mmeas)
		self.bg.addButton(self.mdel)
		self.bg.addButton(self.mmove)

		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl)
		
		self.lbl = QtGui.QLabel("#/row:")
		self.lbl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
		self.hbl.addWidget(self.lbl)
		
		self.nrow = QtGui.QSpinBox(self)
		self.nrow.setObjectName("nrow")
		self.nrow.setRange(1,50)
		self.nrow.setValue(8)
		self.hbl.addWidget(self.nrow)
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vbl.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vbl.addWidget(self.maxs)
		
		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		self.vbl.addWidget(self.brts)
		
		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.vbl.addWidget(self.conts)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		QtCore.QObject.connect(self.nrow, QtCore.SIGNAL("valueChanged(int)"), target.setNPerRow)
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)


	def newMin(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMin(val)

		self.updBC()
		self.busy=0
		
	def newMax(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMax(val)
		self.updBC()
		self.busy=0
	
	def newBrt(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0
		
	def newCont(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0

	def updBC(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		self.mins.setValue(x0)
		self.maxs.setValue(x1)
		self.target.setDenRange(x0,x1)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setLimits(self,lowlim,highlim,curmin,curmax):
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage()
	if len(sys.argv)==1 : window.setData([test_image(),test_image(1),test_image(2),test_image(3)])
	else :
		a=EMData.read_images(sys.argv[1])
		window.setData(a)
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
	

#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import array

class EMImage(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.data=None
		self.datasize=(1,1)
		self.scale=1.0
		self.minden=0
		self.maxden=1.0
		self.mindeng=0
		self.maxdeng=1.0
		self.origin=(0,0)
		self.nperrow=6
		self.nshow=0
		self.mousedrag=None
		
		self.inspector=None
		self.inspectorl=None
		self.inspector3=None
	
	def setData(self,data):
		"""You may pass a single 2D image, a list of 2D images or a single 3D image"""
		self.data=data
		if data==None:
			self.updateGL()
			return
		
		# If we have a list of 2D images
		if isinstance(data,list) :
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
		# If we have a single 2D image
		elif data.get_zsize()==1:
			mean=data.get_attr("mean")
			sigma=data.get_attr("sigma")
			m0=data.get_attr("minimum")
			m1=data.get_attr("maximum")
			
			self.minden=max(m0,mean-3.0*sigma)
			self.maxden=min(m1,mean+3.0*sigma)
			self.mindeng=max(m0,mean-5.0*sigma)
			self.maxdeng=min(m1,mean+5.0*sigma)

			self.datasize=(data.get_xsize(),data.get_ysize())
		# if we have a single 3D image
		elif data.get_zsize()>1 :
			pass
		# Someone passed something wierd
		else :
			self.data=None
			self.updateGL()
			return
		
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
		if isinstance(self.data,list) :
			yo=self.height()-self.origin[1]-1
			self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+yo)-self.height()/2)
		else : self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+self.origin[1])-self.height()/2)
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

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		
		
		if isinstance(self.data,list) :
			if self.nshow==-1 :
				GL.glPixelZoom(1.0,-1.0)
				n=len(self.data)
				x,y=-self.origin[0],self.height()-self.origin[1]-1
				for i in range(n):
					w=int(min(self.data[i].get_xsize()*self.scale,self.width()))
					h=int(min(self.data[i].get_ysize()*self.scale,self.height()))
#					print i,x,y,w,h
					if x>0 and x<self.width() and y>0 and y<self.height() :
						a=self.data[i].render_amp8(0,0,w,h,(w-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
						GL.glRasterPos(x,y)
						GL.glDrawPixels(w,h,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
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
						a=self.data[i].render_amp8(x0,0,x1,y1,(x1-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
#						a=self.data[i].render_amp8(int(-x),int(y-self.height()),min(w,int(x+w)),int(h-y+self.height()),min(w-1,int(x+w-1))/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
						if x<0 : xx=0
						else : xx=x
						if y>=self.height() : yy=self.height()-1
						else : yy=y
						GL.glRasterPos(xx,yy)
						GL.glDrawPixels(x1,y1,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
						
					
					if (i+1)%self.nperrow==0 : 
						y-=h+2.0
						x=-self.origin[0]
					else: x+=w+2.0
			else:
				a=self.data[self.nshow].render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
				GL.glRasterPos(0,self.height()-1)
				GL.glPixelZoom(1.0,-1.0)
				GL.glDrawPixels(self.width(),self.height(),GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
				hist=array.array("I")
				hist.fromstring(a[-1024:])
				if self.inspectorl : self.inspectorl.setHist(hist,self.minden,self.maxden)
		else :
			a=self.data.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
#			a=self.data.render_amp8(self.origin[0],self.origin[1],self.width(),self.height(),(self.width()-1)/4*4+4,1.0,1,254,self.minden,self.maxden,0)
			GL.glRasterPos(0,self.height()-1)
			GL.glPixelZoom(1.0,-1.0)
			GL.glDrawPixels(self.width(),self.height(),GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
#			print "X ",GL.glGetFloatv(GL.GL_CURRENT_RASTER_POSITION),GL.glGetBooleanv(GL.GL_CURRENT_RASTER_POSITION_VALID)
#			print self.width(),self.height(),len(a)
			hist=array.array("I")
			hist.fromstring(a[-1024:])
			if self.inspector : self.inspector.setHist(hist,self.minden,self.maxden)
	
	def resizeGL(self, width, height):
		side = min(width, height)
#		GL.glViewport((width - side) / 2, (height - side) / 2, side, side)
		GL.glViewport(0,0,self.width(),self.height())
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
#		GL.glMatrixMode(GL.GL_PROJECTION)
#		GL.glLoadIdentity()
#		GL.glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
#		GL.glMatrixMode(GL.GL_MODELVIEW)
	
	def showInspector(self,force=0):
		if not force and self.inspector==None and self.inspectorl==None and self.inspector3==None : return
		if isinstance(self.data,list) or isinstance(self.data,tuple) :
			if self.inspector : self.inspector.hide()
			if self.inspector3 : self.inspector3.hide()
			if not self.inspectorl : self.inspectorl=EMImageMxInspector2D(self)
			self.inspectorl.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
			self.inspectorl.show()
		elif self.data.get_zsize()==1 :
			if self.inspectorl : self.inspectorl.hide()
			if self.inspector3 : self.inspector3.hide()
			if not self.inspector : self.inspector=EMImageInspector2D(self)
			self.inspector.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
			self.inspector.show()
		else:
			pass	# 3d not done yet
	
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

class ImgHistogram(QtGui.QWidget):
	def __init__(self,parent):
		QtGui.QWidget.__init__(self,parent)
		self.brush=QtGui.QBrush(Qt.black)
		
		self.font=QtGui.QFont("Helvetica", 10);
		self.probe=None
		self.histdata=None
		self.setMinimumSize(QtCore.QSize(258,128))
	
	def setData(self,data,minden,maxden):
		self.histdata=data
#		self.norm=max(self.histdata)
		self.norm=0
		self.minden=minden
		self.maxden=maxden
		for i in self.histdata: self.norm+=i*i
		self.norm-=max(self.histdata)**2
		self.norm=sqrt(self.norm/255)*3.0
		self.total=sum(self.histdata)
		if self.norm==0 : self.norm=1.0
		if self.total==0 : self.histdata=None
		self.update()
	
	def paintEvent (self, event):
		if not self.histdata : return
		p=QtGui.QPainter()
		p.begin(self)
		p.setPen(Qt.darkGray)
		for i,j in enumerate(self.histdata):
			p.drawLine(i,127,i,127-j*126/self.norm)
		
		# If the user has dragged, we need to show a value
		if self.probe :
			p.setPen(Qt.blue)
			p.drawLine(self.probe[0]+1,0,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setPen(Qt.darkRed)
			p.drawLine(self.probe[0]+1,127,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setFont(self.font)
			p.drawText(200,20,"x=%d"%(self.probe[0]))
			p.drawText(200,34,"%1.2f"%(self.probe[0]/255.0*(self.maxden-self.minden)+self.minden))
			p.drawText(200,48,"y=%d"%(self.probe[1]))
			p.drawText(200,62,"%1.2f%%"%(100.0*float(self.probe[1])/self.total))
		
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
		self.nrow.setValue(6)
		self.hbl.addWidget(self.nrow)
		
		self.lbl = QtGui.QLabel("N:")
		self.lbl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
		self.hbl.addWidget(self.lbl)
		
		self.imgn = QtGui.QSpinBox(self)
		self.imgn.setObjectName("imgn")
		self.imgn.setRange(-1,50)
		self.imgn.setValue(0)
		self.imgn.setSpecialValueText("All")
		self.hbl.addWidget(self.imgn)
		
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
		QtCore.QObject.connect(self.imgn, QtCore.SIGNAL("valueChanged(int)"), target.setNShow)
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
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value+self.brts.value))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value-self.brts.value))
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


class EMImageInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist)
		
		self.vbl2 = QtGui.QVBoxLayout(self)
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
		
		self.invtog = QtGui.QPushButton("Invert")
		self.invtog.setCheckable(1)
		self.vbl2.addWidget(self.invtog)
		
		self.ffttog = QtGui.QPushButton("FFT")
		self.ffttog.setCheckable(1)
		self.vbl2.addWidget(self.ffttog)

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
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value+self.brts.value))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value-self.brts.value))
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
	if len(sys.argv)==1 : window.setData(test_image(size=(512,512)))
	else :
		a=EMData.read_images(sys.argv[1])
		if len(a)==1 : window.setData(a[0])
		else : window.setData(a)
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
	

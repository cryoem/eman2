#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu)
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

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer

class EMImageMX(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display stacks of 2D images
	in 'matrix' form on the display. The middle mouse button will bring up a
	control-panel. The QT event loop must be running for this object to function
	properly.
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
		self.gamma=1.0
		self.origin=(0,0)
		self.nperrow=8
		self.nshow=-1
		self.mousedrag=None
		self.nimg=0
		self.changec={}
		self.mmode=0
		self.targetorigin=None
		
		self.coords=[]
		self.valstodisp=["Img #"]
		self.setAcceptDrops(True)
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		
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
		self.timer.start(25)
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setOrigin(self,x,y):
		"""Set the display origin within the image"""
		self.origin=(x,y)
		self.targetorigin=None
		self.updateGL()
		
	def setScale(self,newscale):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		
		if self.targetorigin : 
			self.origin=self.targetorigin
			self.targetorigin=None
			
		if self.data and len(self.data)>0 and (self.data[0].get_ysize()*newscale>self.height() or self.data[0].get_xsize()*newscale>self.width()):
			newscale=min(float(self.height())/self.data[0].get_ysize(),float(self.width())/self.data[0].get_xsize())
			
#		yo=self.height()-self.origin[1]-1
		yo=self.origin[1]
#		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+yo)-self.height()/2)
#		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(yo-self.height()/2)+self.height()/2)
		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+self.origin[1])-self.height()/2)
#		print self.origin,newscale/self.scale,yo,self.height()/2+yo
		
		self.scale=newscale
		self.updateGL()
		
	def setDenMin(self,val):
		self.minden=val
		self.updateGL()
		
	def setDenMax(self,val):
		self.maxden=val
		self.updateGL()

	def setGamma(self,val):
		self.gamma=val
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
	
	def timeout(self):
		"""Called a few times eeach second when idle for things like automatic scrolling"""
		if self.targetorigin :
			vec=(self.targetorigin[0]-self.origin[0],self.targetorigin[1]-self.origin[1])
			h=hypot(vec[0],vec[1])
			if h<25.0 :
				self.origin=self.targetorigin
				self.targetorigin=None
			vec=(vec[0]/h,vec[1]/h)
			self.origin=(self.origin[0]+vec[0]*10.0,self.origin[1]+vec[1]*10.0)
			self.updateGL()
		
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		for i in self.data:
			self.changec[i]=i.get_attr("changecount")
		
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)
		
#		GL.glPixelZoom(1.0,-1.0)
		n=len(self.data)
		x,y=-self.origin[0],-self.origin[1]
		hist=numpy.zeros(256)
		if len(self.coords)>n : self.coords=self.coords[:n]
		for i in range(n):
			w=int(min(self.data[i].get_xsize()*self.scale,self.width()))
			h=int(min(self.data[i].get_ysize()*self.scale,self.height()))
			if x<self.width() and y<self.height():
				if x>=0 and y>=0:
					a=self.data[i].render_amp8(0,0,w,h,(w-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,6)
					GL.glRasterPos(x,y)
					GL.glDrawPixels(w,h,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
					hist2=numpy.fromstring(a[-1024:],'i')
					hist+=hist2
					
					ty=y
					for v in self.valstodisp:
						if v=="Img #" : self.renderText(x,ty,"%d"%i)
						else : 
							av=self.data[i].get_attr(v)
							if isinstance(av,float) : avs="%1.4g"%av
							else: avs=str(av)
							try: self.renderText(x,ty,str(avs))
							except: self.renderText(x,ty,"------")
						ty+=16
				elif x+w>0 and y+h>0:
					tx=int(max(x,0))
					ty=int(max(y,0))
					tw=int(w-tx+x)
					th=int(h-ty+y)
					a=self.data[i].render_amp8(int(-min(x/self.scale,0)),int(-min(y/self.scale,0)),tw,th,(tw-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,6)
					GL.glRasterPos(tx,ty)
					GL.glDrawPixels(tw,th,GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
					hist2=numpy.fromstring(a[-1024:],'i')
					hist+=hist2
				
			
			try: self.coords[i]=(x+self.origin[0],y+self.origin[1],self.data[i].get_xsize()*self.scale,self.data[i].get_ysize()*self.scale)
			except: self.coords.append((x+self.origin[0],y+self.origin[1],self.data[i].get_xsize()*self.scale,self.data[i].get_ysize()*self.scale))
			
			if (i+1)%self.nperrow==0 : 
				y+=h+2.0
				x=-self.origin[0]
			else: x+=w+2.0
		if self.inspector : self.inspector.setHist(hist,self.minden,self.maxden)
	
	def renderText(self,x,y,s):
		GL.glRasterPos(x+2,y+2)
		for c in s:
			GLUT.glutBitmapCharacter(GLUT.GLUT_BITMAP_9_BY_15,ord(c))

	
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
		#print width/(self.data[0].get_xsize()*self.scale)
		self.nperrow=int(width/(self.data[0].get_xsize()*self.scale))
		if self.nperrow<1 : self.nperrow=1
		#except: pass
		
		if self.data and len(self.data)>0 and (self.data[0].get_ysize()*self.scale>self.height() or self.data[0].get_xsize()*self.scale>self.width()):
			self.scale=min(float(self.height())/self.data[0].get_ysize(),float(self.width())/self.data[0].get_xsize())
		
	def scrollTo(self,n):
		"""Moves image 'n' to the center of the display"""
#		print self.origin,self.coords[0],self.coords[1]
#		try: self.origin=(self.coords[n][0]-self.width()/2,self.coords[n][1]+self.height()/2)
#		try: self.origin=(self.coords[8][0]-self.width()/2-self.origin[0],self.coords[8][1]+self.height()/2-self.origin[1])
		try: self.targetorigin=(self.coords[n][0]-self.width()/2+self.data[0].get_xsize()*self.scale/2,self.coords[n][1]-self.height()/2+self.data[0].get_ysize()*self.scale/2)
		except: return
#		print n,self.origin
#		self.updateGL()
		
	def setValDisp(self,v2d):
		"""Pass in a list of strings describing image attributes to overlay on the image, in order of display"""
		v2d.reverse()
		self.valstodisp=v2d
		self.updateGL()
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		if not self.inspector : self.inspector=EMImageMxInspector2D(self)
		self.inspector.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
		self.inspector.show()

	def scrtoimg(self,vec):
		"""Converts screen location (ie - mouse event) to pixel coordinates within a single
		image from the matrix. Returns (image number,x,y) or None if the location is not within any
		of the contained images. """
		absloc=((vec[0]+self.origin[0]),(self.height()-(vec[1]-self.origin[1])))
		for i,c in enumerate(self.coords):
			if absloc[0]>c[0] and absloc[1]>c[1] and absloc[0]<c[0]+c[2] and absloc[1]<c[1]+c[3] :
				return (i,(absloc[0]-c[0])/self.scale,(absloc[1]-c[1])/self.scale)
		return None
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
		
	def dragEnterEvent(self,event):
#		f=event.mimeData().formats()
#		for i in f:
#			print str(i)
		
		if event.source()==self:
			event.setDropAction(Qt.MoveAction)
			event.accept()
		elif event.provides("application/x-eman"):
			event.setDropAction(Qt.CopyAction)
			event.accept()

	
	def dropEvent(self,event):
		lc=self.scrtoimg((event.pos().x(),event.pos().y()))
		if event.source()==self:
#			print lc
			n=int(event.mimeData().text())
			if not lc : lc=[len(self.data)]
			if n>lc[0] : 
				self.data.insert(lc[0],self.data[n])
				del self.data[n+1]
			else : 
				self.data.insert(lc[0]+1,self.data[n])
				del self.data[n]
			event.setDropAction(Qt.MoveAction)
			event.accept()
		elif event.provides("application/x-eman"):
			x=loads(event.mimeData().data("application/x-eman"))
			if not lc : self.data.append(x)
			else : self.data.insert(lc[0],x)
			self.setData(self.data)
			event.acceptProposedAction()


	def mousePressEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
#		print lc
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.RightButton:
			self.mousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
			if self.mmode=="drag" and lc:
				xs=self.data[lc[0]].get_xsize()
				ys=self.data[lc[0]].get_ysize()
				drag = QtGui.QDrag(self)
				mimeData = QtCore.QMimeData()
				mimeData.setData("application/x-eman", dumps(self.data[lc[0]]))
				mimeData.setText( str(lc[0])+"\n")
				di=QImage(self.data[lc[0]].render_amp8(0,0,xs,ys,xs*4,1.0,0,255,self.minden,self.maxden,1.0,14),xs,ys,QImage.Format_RGB32)
				mimeData.setImageData(QtCore.QVariant(di))
				drag.setMimeData(mimeData)

# This (mini image drag) looks cool, but seems to cause crashing sometimes in the pixmap creation process  :^(
				#di=QImage(self.data[lc[0]].render_amp8(0,0,xs,ys,xs*4,1.0,0,255,self.minden,self.maxden,14),xs,ys,QImage.Format_RGB32)
				#if xs>64 : pm=QtGui.QPixmap.fromImage(di).scaledToWidth(64)
				#else: pm=QtGui.QPixmap.fromImage(di)
				#drag.setPixmap(pm)
				#drag.setHotSpot(QtCore.QPoint(12,12))
				
				dropAction = drag.start()
#				print dropAction
			
			elif self.mmode=="del" and lc:
				del self.data[lc[0]]
				self.setData(self.data)

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
		
		self.vals = QtGui.QMenu()
		self.valsbut = QtGui.QPushButton("Values")
		self.valsbut.setMenu(self.vals)
		
		try:
			self.vals.clear()
			vn=self.target.data[0].get_attr_dict().keys()
			vn.sort()
			for i in vn:
				action=self.vals.addAction(i)
				action.setCheckable(1)
				action.setChecked(0)
		except:
			pass
		
		action=self.vals.addAction("Img #")
		action.setCheckable(1)
		action.setChecked(1)
		
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
		
		#self.mmeas = QtGui.QPushButton("Meas")
		#self.mmeas.setCheckable(1)
		#self.hbl2.addWidget(self.mmeas)

		self.mdel = QtGui.QPushButton("Del")
		self.mdel.setCheckable(1)
		self.hbl2.addWidget(self.mdel)

		self.mdrag = QtGui.QPushButton("Drag")
		self.mdrag.setCheckable(1)
		self.hbl2.addWidget(self.mdrag)

		self.bg=QtGui.QButtonGroup()
		self.bg.setExclusive(1)
#		self.bg.addButton(self.mmeas)
		self.bg.addButton(self.mdel)
		self.bg.addButton(self.mdrag)

		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl)
		
		self.hbl.addWidget(self.valsbut)
		
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
		
		self.gammas = ValSlider(self,(.5,2.0),"Gam:")
		self.gammas.setObjectName("gamma")
		self.gammas.setValue(1.0)
		self.vbl.addWidget(self.gammas)

		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		QtCore.QObject.connect(self.vals, QtCore.SIGNAL("triggered(QAction*)"), self.newValDisp)
		QtCore.QObject.connect(self.nrow, QtCore.SIGNAL("valueChanged(int)"), target.setNPerRow)
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)
		QtCore.QObject.connect(self.gammas, QtCore.SIGNAL("valueChanged"), self.newGamma)
		
		#QtCore.QObject.connect(self.mmeas, QtCore.SIGNAL("clicked(bool)"), self.setMeasMode)
		QtCore.QObject.connect(self.mdel, QtCore.SIGNAL("clicked(bool)"), self.setDelMode)
		QtCore.QObject.connect(self.mdrag, QtCore.SIGNAL("clicked(bool)"), self.setDragMode)

	def newValDisp(self):
		v2d=[str(i.text()) for i in self.vals.actions() if i.isChecked()]
		self.target.setValDisp(v2d)

	def setMeasMode(self,i):
		self.target.mmode="meas"
	
	def setDelMode(self,i):
		self.target.mmode="del"
	
	def setDragMode(self,i):
		self.target.mmode="drag"

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
	
	def newGamma(self,val):
		if self.busy : return
		self.busy=1
		self.target.setGamma(val)
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
	window = EMImageMX()
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
	

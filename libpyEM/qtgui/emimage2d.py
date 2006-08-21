#!/bin/env python

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import Numeric
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary

class EMImage2D(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage2D.allim[self]=0
		
# 		try: 
# 			if EMImage2D.gq : pass
# 		except:
# 			EMImage2D.gq=GLU.gluNewQuadric()
# 			GLU.gluQuadricDrawStyle(EMImage2D.gq,GLU.GLU_FILL)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_SMOOTH)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_NONE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_OUTSIDE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_INSIDE)
# 			GLU.gluQuadricTexture(EMImage2D.gq,GL.GL_FALSE)

		
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
		self.rmousedrag=None
		self.mmode=0
		
		self.shapes={}
		self.shapechange=1
		
		self.inspector=None
		
		if image : 
			self.setData(image)
			self.show()
	
	def setData(self,data):
		"""You may pass a single 2D image, a list of 2D images or a single 3D image"""
		if not self.data and data: self.resize(data.get_xsize(),data.get_ysize())
		
		self.data=data
		if data==None:
			self.updateGL()
			return
		
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		
		self.minden=max(m0,mean-3.0*sigma)
		self.maxden=min(m1,mean+3.0*sigma)
		self.mindeng=max(m0,mean-5.0*sigma)
		self.maxdeng=min(m1,mean+5.0*sigma)

		self.datasize=(data.get_xsize(),data.get_ysize())
	
		if self.fft : 
			self.setFFT(1)
		
		self.showInspector()		# shows the correct inspector if already open
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		self.updateGL()
	
	def setDenMin(self,val):
		self.minden=val
		self.updateGL()
		
	def setDenMax(self,val):
		self.maxden=val
		self.updateGL()

	def setOrigin(self,x,y):
		"""Set the display origin within the image"""
		self.origin=(x,y)
		self.updateGL()
		
	def setScale(self,newscale):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+self.origin[1])-self.height()/2)
		self.scale=newscale
		self.updateGL()
		
	def setInvert(self,val):
		if val: self.invert=1
		else : self.invert=0
		self.updateGL()
		
	def setFFT(self,val):
		if val:
			try: 
				self.fft=self.data.do_fft()
				self.fft.set_value_at(0,0,0,0)
				self.fft.set_value_at(1,0,0,0)
				self.fft.process_inplace("eman1.xform.fourierorigin",{})
				self.fft=self.fft.get_fft_amplitude()
			
				mean=self.fft.get_attr("mean")
				sigma=self.fft.get_attr("sigma")
				m0=self.fft.get_attr("minimum")
				m1=self.fft.get_attr("maximum")
			
				self.fminden=0
				self.fmaxden=min(m1,mean+5.0*sigma)
				self.fmindeng=0
				self.fmaxdeng=min(m1,mean+8.0*sigma)
				
				self.ominden=self.minden
				self.omaxden=self.maxden
				
				self.showInspector()
			except: 
				self.fft=None
		else: 
			self.fft=None
			self.minden=self.ominden
			self.maxden=self.omaxden
			self.showInspector()
		self.updateGL()

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
	
	def paintGL(self):
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
#		GL.glLoadIdentity()
#		GL.glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		
		if self.shapechange:
			self.setupShapes()
			self.shapechange=0
		
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)
		
		if self.fft : a=self.fft.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,2)
		else : a=self.data.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,2)
		GL.glRasterPos(0,self.height()-1)
		GL.glPixelZoom(1.0,-1.0)
		GL.glDrawPixels(self.width(),self.height(),GL.GL_LUMINANCE,GL.GL_UNSIGNED_BYTE,a)
		hist=Numeric.fromstring(a[-1024:],'i')
		if self.inspector : 
			if self.invert: self.inspector.setHist(hist,self.maxden,self.minden) 
			else: self.inspector.setHist(hist,self.minden,self.maxden)
	
		GL.glPushMatrix()
		GL.glTranslate(-self.origin[0],-self.origin[1],0)
		GL.glScale(self.scale,self.scale,self.scale)
		GL.glCallList(1)
		GL.glPopMatrix()
		self.changec=self.data.get_attr("changecount")
				
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	def setupShapes(self):
		# make our own cirle rather than use gluDisk or somesuch
		GL.glNewList(2,GL.GL_COMPILE)
		GL.glBegin(GL.GL_LINE_LOOP)
		d2r=pi/180.0
		for i in range(90): GL.glVertex(sin(i*d2r*4.0),cos(i*d2r*4.0))
		GL.glEnd()
		GL.glEndList()
		
		GL.glNewList(1,GL.GL_COMPILE)
		GL.glPushMatrix()
		for k,s in self.shapes.items():
			GL.glColor(s[1],s[2],s[3])
			
			if s[0]=="rect":
				GL.glLineWidth(s[8])
				GL.glBegin(GL.GL_LINE_LOOP)
				GL.glVertex(s[4],s[5])
				GL.glVertex(s[6],s[5])
				GL.glVertex(s[6],s[7])
				GL.glVertex(s[4],s[7])
				GL.glEnd()
			elif s[0]=="line":
				GL.glPushMatrix()
				GL.glColor(s[1],s[2],s[3])
				GL.glLineWidth(s[8])
				GL.glBegin(GL.GL_LINES)
				GL.glVertex(s[4],s[5])
				GL.glVertex(s[6],s[7])
				GL.glEnd()
				GL.glPopMatrix()
			elif s[0]=="label":
				if s[8]<0 :
					GL.glPushMatrix()
					GL.glColor(0,0,0)
					GL.glTranslate(s[4],s[5],0)
					GL.glScale(s[7]/100.0/self.scale,s[7]/100.0/self.scale,s[7]/100.0/self.scale)
					GL.glLineWidth(-s[8]+2)
					for i in s[6]:
						GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_ROMAN,ord(i))
					GL.glPopMatrix()
				
				GL.glPushMatrix()
				GL.glColor(s[1],s[2],s[3])
				GL.glTranslate(s[4],s[5],0)
#				GL.glScale(s[7]/100.0,s[7]/100.0,s[7]/100.0)
				GL.glScale(s[7]/100.0/self.scale,s[7]/100.0/self.scale,s[7]/100.0/self.scale)
				GL.glLineWidth(fabs(s[8]))
				for i in s[6]:
					GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_ROMAN,ord(i))
				GL.glPopMatrix()
			elif s[0]=="circle":
				GL.glPushMatrix()
				GL.glColor(s[1],s[2],s[3])
				GL.glLineWidth(s[7])
				GL.glTranslate(s[4],s[5],0)
				GL.glScale(s[6],s[6],s[6])
				GL.glCallList(2)
				GL.glPopMatrix()

			
		GL.glPopMatrix()
		GL.glEndList()
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector2D(self)
		if self.fft : self.inspector.setLimits(self.fmindeng,self.fmaxdeng,self.fminden,self.fmaxden)
		else : self.inspector.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
		self.inspector.show()
	
	def setMMode(self,m):
		self.mmode=m
	
	def addShape(self,k,s):
		"""Add a 'shape' object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed.
		
		0         1  2  3  4  5     6     7     8
		"rect"    R  G  B  x0 y0    x1    y1    linew
		"line"    R  G  B  x0 y0    x1    y1    linew
		"label"   R  G  B  x0 y0    text  size	linew
		"circle"  R  G  B  x0 y0    r     linew
		"""
		self.shapes[k]=s
		self.shapechange=1
		self.updateGL()
	
	def addShapes(self,d):
		self.shapes.update(d)
		self.shapechange=1
		self.updateGL()
	
	def delShapes(self,k):
		try:
			for i in k:
				del a[k]
		except: del a[k]
		self.shapechange=1
		self.updateGL()
	
	def scrtoimg(self,vec):
		return ((vec[0]+self.origin[0])/self.scale,(self.height()-(vec[1]-self.origin[1]))/self.scale)
	
	def mousePressEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.RightButton:
			self.rmousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
			elif self.mmode==1 :
				try: 
					del self.shapes["MEASL"]
				except: pass
				self.addShape("MEAS",("line",.5,.1,.5,lc[0],lc[1],lc[0]+1,lc[1],2))
	
	def mouseMoveEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
		if self.rmousedrag and event.buttons()&Qt.RightButton:
			self.origin=(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
			self.rmousedrag=(event.x(),event.y())
			self.update()
		elif event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
			elif self.mmode==1 :
				self.addShape("MEAS",("line",.5,.1,.5,self.shapes["MEAS"][4],self.shapes["MEAS"][5],lc[0],lc[1],2))
				dx=lc[0]-self.shapes["MEAS"][4]
				dy=lc[1]-self.shapes["MEAS"][5]
				self.addShape("MEASL",("label",.5,.1,.5,lc[0]+2,lc[1]+2,"%d,%d - %d,%d\n%1.1f,%1.1f (%1.2f)"%(self.shapes["MEAS"][4],self.shapes["MEAS"][5],lc[0],lc[1],dx,dy,hypot(dx,dy)),10,-1))
	
	def mouseReleaseEvent(self, event):
		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.RightButton:
			self.rmousedrag=None
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mouseup"), event)
				return
			elif self.mmode==1 :
				self.addShape("MEAS",("line",.5,.1,.5,self.shapes["MEAS"][4],self.shapes["MEAS"][5],lc[0],lc[1],2))


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
		
		self.vbl2 = QtGui.QVBoxLayout()
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

		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)
		
		self.mapp = QtGui.QPushButton("App")
		self.mapp.setCheckable(1)
		self.hbl2.addWidget(self.mapp)

		self.mmeas = QtGui.QPushButton("Meas")
		self.mmeas.setCheckable(1)
		self.hbl2.addWidget(self.mmeas)

		self.mmode=QtGui.QButtonGroup()
		self.mmode.setExclusive(1)
		self.mmode.addButton(self.mapp,0)
		self.mmode.addButton(self.mmeas,1)
		
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
		QtCore.QObject.connect(self.invtog, QtCore.SIGNAL("toggled(bool)"), target.setInvert)
		QtCore.QObject.connect(self.ffttog, QtCore.SIGNAL("toggled(bool)"), target.setFFT)
		QtCore.QObject.connect(self.mmode, QtCore.SIGNAL("buttonClicked(int)"), target.setMMode)

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
	window = EMImage2D()
	if len(sys.argv)==1 : 
		window.setData(test_image(size=(512,512)))

		# these lines are for testing shape rendering
# 		window.addShape("a",["rect",.2,.8,.2,20,20,80,80,2])
# 		window.addShape("b",["circle",.5,.8,.2,120,50,30.0,2])
# 		window.addShape("c",["line",.2,.8,.5,20,120,100,200,2])
# 		window.addShape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	else :
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
	

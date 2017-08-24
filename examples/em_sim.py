#!/usr/bin/env python
# Muyuan Chen 2017-08
import numpy as np

from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtGui, QtCore, QtOpenGL
from PyQt4.QtCore import Qt
from emimage2d import EMImage2DWidget
from EMAN2 import *
from emapplication import EMApp


class Microscope(QtOpenGL.QGLWidget):
	
	def __init__(self, parent=None, imgwindow=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
		self.win_size=[500,1000]
		self.setMinimumSize(self.win_size[0], self.win_size[1])
		self.imgwindow=imgwindow
		
		#### set up the first two lens
		self.source=s=.9
		d0=.15
		d1=.25
		f0=.085
		#f0=(d0*d1-d0*f1)/(d1-f1+d0)
		f1=(d0*d1-f0*d1-f0*d0)/(d0-f0)
		
		### y_pos, focal_pos
		self.lens=[
			[s-d0, f0],
			[s-d0-d1, f1],
			
			[0.1, -1], ### f==-1 : specimen stage
			
			[-0.03, .07],
			[-0.3, .07],
			[-0.56, .07],
			[-.9, -2] ### f==-2 : detector
			]
		
		self.tmp_pts=[]
		self.drag_lens=-1
		self.hold_shift=False
		self.beam_dist=0
		self.defocus=0
		self.beam_lastxy=[0,0]
		
		if imgwindow:
			imgwindow.show()
			img=self.draw_wave()
			self.imgwindow.set_data(from_numpy(img))

	def initializeGL(self):
		glClearColor(1.0, 1.0, 1.0, 1.0)
		glClearDepth(1.0)
		glLineStipple(1, 0xAAAA);

		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(40.0, 1.0, 1.0, 30.0)
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable( GL_LINE_SMOOTH );
		glEnable( GL_POLYGON_SMOOTH );
		glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
		glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

	def resizeGL(self, w, h):
	
		glViewport(0, 0, w, h)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(40.0, 1.0, 1.0, 30.0)
		self.win_size=[w,h]

	def paintGL(self):

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glLoadIdentity()
		
		
		
		for l in self.lens:
			self.draw_lens(l[0], l[1]);
			
		
		#glColor3f( 1., .5, .5 );
		#glLineWidth(2.)
		#glBegin(GL_LINE_STRIP)
		#for p in self.tmp_pts:	
			#glVertex(p[0],p[1])
		#glEnd()
		
		self.draw_path()
		glFlush()
		
	def draw_wave(self):
		sz=256
		raw=np.zeros(sz)
		wid=2
		w2=4
		raw[sz/2-wid-w2:sz/2+wid-w2]=1
		raw[sz/2-wid+w2:sz/2+wid+w2]=1

		wavelen=2.
		
		
		mult=400
#		pln=np.round(.13*mult)
#		lens=np.round(np.array([[.27,.07], [.26, .09], [.34, .07]])*mult)
#		print pln, lens
		l=np.array(self.lens)[3:]
		pln=np.round((self.lens[2][0]-l[0,0])*mult)
		lens=np.round(mult*np.array([[l[i,0]-l[i+1,0], l[i,1]] for i in range(len(l)-1)]))
#		print pln,lens
		
		imgs=[]

		ix=np.arange(sz)
		iz=np.arange(1,pln)[:,None,None]

		dst=np.sqrt((ix-ix[:,None])**2 +iz**2)
		cpx=raw[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)*(1/dst**2)
		img=np.sum(cpx, axis=1)
		imgs.append(img)
		for il,ln in enumerate(lens):
			f=ln[1]
			ps=((ix-sz/2)**2)/(f*2)*(2*np.pi/wavelen)
			proj=imgs[-1][-1]
			proj_ps=proj*np.exp(-1j*(np.pi-ps))

			img0=[]
			zmax=ln[0]
			iz=np.arange(1,zmax)[:,None,None]

			dst=np.sqrt((ix-ix[:,None])**2 +iz**2)
			cpx1=proj_ps[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)*(1/dst**2)

			img=np.sum(cpx1, axis=1)

			imgs.append(img)
		
		final=np.vstack(imgs)
		final/=np.sum(abs(final),axis=1)[:,None]
		return abs(final)[::-1,:].copy()







	
	#### draw vertex and keep track of the distance
	def draw_vertex(self, x, y):
		glVertex(x, y, 0)
		lx,ly=self.beam_lastxy
		if self.beam_dist<0: ### start tracking
			self.beam_dist=0
		else: ### add distance
			self.beam_dist+=np.sqrt((lx-x)**2+(ly-y)**2)
		self.beam_lastxy=[x,y]
	
	def draw_path(self):
		
		src=[self.source] ## input position
		theta=[np.tan(np.pi/3.)] ## input angle
		diverge=[False for s in src] ### whether the beam at previous layer diverges
		
		bign=1e10 ### a big number. so we do not reach inf when the beam is paralleled
		
		
		beamstop=False ### stop the beam at some point..
		scatter=False ### scattered beam from last stage
		ymax=self.lens[-1][0]
		
		self.beam_dist=0
		dist=[0] ### track beam distance
		for il,l in enumerate(self.lens):
			news=[]
			newt=[]
			newv=[]
			
			for si,s in enumerate(src):
				
				if s==None: ### no beam
					news.append(None)
					newt.append(0)
					newv.append(False)
					continue
				
				self.beam_dist=-1
				if l[1]==-2: ### detector stage always stops the beam
					beamstop=True
					maxwidth=np.inf ### max width of the column
				else:
					beamstop=False
					maxwidth=.6 ### max width of the column
					
				tt=theta[si]
				dvg=diverge[si]
				
			
				#### start from the focal point at each iteration
				d0=s-l[0]
			
				#### deal with parallel beams
				if l[1]<0 or abs(1./l[1]-1./d0)<1./bign:
					d1=bign
				else:
					d1=1./(1./l[1]-1./d0)
				
					
				if abs(d0)>bign:
					d0=np.sign(d0)*bign
					
				w=tt*d0 ### width of beam when it hit the lens
				tt=w/d1 ### output angle
				
				
				#### deal with beam scattering
				if not scatter:
					s0=s
					x0=0
				else:
					## beam cannot be scattered at the first stage.
					x0=s-self.lens[il-1][0]
					s0=self.lens[il-1][0]
					#print si, x0, s0, d0
				
				ym=s-d0/abs(w)
				xm=1
				if ym<ymax:
					xm=xm/ym*ymax
					ym=ymax

				glColor3f( .5, 1., .2 + (si>0)*.7 );
				glBegin(GL_LINES)
				if d1>0: ### output beam converge to next focal point
					
					if (scatter and abs(x0)<bign) or (d0>0 and dvg==False): 
						#### lens below last focal point, draw input beam
						
						if abs(w)>maxwidth:
							#### input beam miss this lens. stop the beam
							
							self.draw_vertex(-x0, s0)
							self.draw_vertex(-xm, ym)
							
							self.draw_vertex(x0, s0)
							self.draw_vertex(xm, ym)
							beamstop=True
						
						else:
							#### input beam hit the lens
							self.draw_vertex(-x0, s0)
							self.draw_vertex(-w, l[0])
							
							self.draw_vertex(x0, s0)
							self.draw_vertex(w, l[0])
					
					if beamstop==False:
					
						#### draw output beam
						self.draw_vertex(-w, l[0])
						
						if il<len(self.lens)-1 and l[0]-d1<=self.lens[il+1][0]:
							#### output beam stop at the next lens
							
							l1=self.lens[il+1][0]
							
							d2=l[0]-l1
							w1=w-tt*d2
							
							self.draw_vertex(-w1, l1)
							self.draw_vertex(w1, l1)
						else:
							#### output beam extend to next focal point
							self.draw_vertex(0, l[0]-d1)
							self.draw_vertex(0, l[0]-d1)
							
						self.draw_vertex(w, l[0])
					dvg=False
				
				else:
					#### output beam diverging
					dvg=True
					if il<len(self.lens)-1:
						l1=self.lens[il+1][0]
						mult=(l1-l[0])/d1-1
						if w*(2+mult)>maxwidth:
							#mult=abs(1./w)
							mult=(ymax-l[0])/d1-1
							#print mult
							beamstop=True
					else:
						mult=abs(1./w)
					
					for sign in [-1,1]:

						self.draw_vertex(sign*x0, s0)
						if abs(w)>maxwidth: ### input beam miss next lens
							self.draw_vertex(sign*xm, ym)
						
						else: ### output beam hit next lens
							self.draw_vertex(sign*w, l[0])
							
							self.draw_vertex(sign*w, l[0])
							self.draw_vertex(sign*w*(2+mult), l[0]+d1*(1+mult))
					
				glEnd()
				s=l[0]-d1
				dist[si]+=self.beam_dist
				if l[1]==-1: ### specimen stage, generate scattered and unscatterd beams
					
					news.append(s)
					newt.append(tt)
					newv.append(False)
					
					news.append(l[0]+abs(w))
					newt.append(np.tan(np.pi/4.))
					newv.append(False)
					dist.append(dist[si])
					
					news.append(l[0]-abs(w))
					newt.append(np.tan(np.pi/4.))
					newv.append(False)
					dist.append(dist[si])
					scatter=True
				
				else:
					if not beamstop:
						news.append(s)
						newt.append(tt)
						newv.append(dvg)
					else:
						news.append(None)
						newt.append(0)
						newv.append(False)
			src=news
			theta=newt
			diverge=newv
			if l[1]>=0:
				scatter=False
		#print ", ".join(["{:.2f}".format(d) for d in dist]), abs(dist[0]*2-dist[1]-dist[2])
		#self.defocus=dist[0]*2-dist[1]-dist[2]
		
	
	def draw_lens(self, y=0, focal=.1, scale=1.2):
		
		pts=[]
		
		if focal>=0: ### real lens
			arc=(1./focal)/200.
			t=np.arange(-np.pi*arc, np.pi*arc+1e-5, np.pi*arc/20.)+np.pi/2.
			p0=np.vstack([np.cos(t), np.sin(t)]).T
			p0[:,1]-=p0[0,1]
			l=p0[0,0]-p0[-1,0]
			n=len(p0)
			p1=np.vstack([np.cos(t[::-1]), -np.sin(t[::-1])]).T
			p1[:,1]-=p1[0,1]
			pts=np.vstack([p0,p1])
			pts=pts/l*scale
			pts+=[0,y]
			#print pts
			glColor3f( .5, .5, 1. );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts.tolist())
			glDrawArrays(GL_LINE_STRIP, 0, len(pts))
			
			
			glPushAttrib(GL_ENABLE_BIT)
			glEnable(GL_LINE_STIPPLE)
			glBegin(GL_LINES)
			glVertex(pts[0][0], pts[0][1], 0)
			glVertex(pts[n][0], pts[n][1], 0)
			glEnd()
			glPopAttrib()
		
		elif focal==-1: ### specimen stage
			dx=scale/2.
			dy=.02
			pts=np.array([[dx,dy],[-dx,dy],[-dx,-dy], [dx,-dy], [dx,dy]])
			pts+=[0,y]
			glColor3f( 1., .8, 1. )
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts.tolist())
			glDrawArrays(GL_POLYGON, 0, len(pts))
			glColor3f( 1., .5, 1. )
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts.tolist())
			glDrawArrays(GL_LINE_STRIP, 0, len(pts))
			
		elif focal==-2: ### detector
			dx=scale/2.
			dy=.02
			pts=np.array([[dx,dy],[-dx,dy],[-dx,0], [dx,0], [dx,dy]])
			pts+=[0,y]
			glColor3f( .6, .6, .6 )
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts.tolist())
			glDrawArrays(GL_POLYGON, 0, len(pts))
			glColor3f( .1, .1, .1 )
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts.tolist())
			glDrawArrays(GL_LINE_STRIP, 0, len(pts))
			
			glColor3f( .5, .5, .5);
			glBegin(GL_LINES)
			glVertex(-2, y, 0)
			glVertex(2, y, 0)
			glEnd()
			
			
		return pts
		
	def scr_to_img(self, pp):
		winsz=np.array(self.win_size, dtype=float)
		p0=np.array([pp.x()-winsz[0]/2., winsz[1]/2.-pp.y()])
		p=p0/winsz*2.
		return p
		
	def mousePressEvent(self, QMouseEvent):
		p=self.scr_to_img(QMouseEvent.pos())
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			self.hold_shift=True
		else:
			self.hold_shift=False
			
		#print p
		for i,l in enumerate(self.lens):
			#print abs(p[1]-l[0])
			if abs(p[1]-l[0])<0.05:
				self.drag_lens=i
				
				break
			
		if self.lens[self.drag_lens][1]<=0:
			self.hold_shift=False
		#self.tmp_pts.append(p.tolist())
		#self.updateGL();


	def mouseMoveEvent(self, QMouseEvent):
		if self.drag_lens<0:
			return
		else:
			p=self.scr_to_img(QMouseEvent.pos())
			if self.hold_shift:
				dy=abs(p[1]-self.lens[self.drag_lens][0])
				dy=min(.18, dy)
				dy=0.01/max(.01, dy)
				self.lens[self.drag_lens][1]=dy
			else:
				err=.06
				ypos=p[1]
				if self.drag_lens>0:
					ypos=min(ypos, self.lens[self.drag_lens-1][0]-err)
				else:
					ypos=min(ypos, self.source-err)
				if self.drag_lens<len(self.lens)-1:
					ypos=max(ypos, self.lens[self.drag_lens+1][0]+err)

				self.lens[self.drag_lens][0]=ypos
				
			
			#### constrain parallel beam after first two lens
			d0=self.source-self.lens[0][0]
			d1=self.lens[0][0]-self.lens[1][0]
			
			f0=self.lens[0][1]
			f1=self.lens[1][1]
			if self.drag_lens==0:
				f=(d0*d1-f0*d1-f0*d0)/(d0-f0)
				self.lens[1][1]=max(0.02, f)
			elif self.drag_lens==1:
				f=(d0*d1-d0*f1)/(d1-f1+d0)
				self.lens[0][1]=max(0.02, f)
				
			self.updateGL()
		
	
	
	def mouseReleaseEvent(self, QMouseEvent):
		if self.drag_lens>=0:
			l=self.lens[self.drag_lens]
			print "lens {:d}: py={:.2f}, f={:.2f}".format(
				self.drag_lens, l[0], l[1])
			self.drag_lens=-1
			
			if self.imgwindow:
				img=self.draw_wave()
				self.imgwindow.set_data(from_numpy(img))



class MainWindow(QtGui.QMainWindow):
	
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		
		self.imgview = EMImage2DWidget()
		img=np.random.rand(512,512)
		self.imgview.set_data(from_numpy(img))
		
		widget = Microscope(self, self.imgview)    
		self.setCentralWidget(widget)

		
		

def main():
	app = EMApp()
#	app=QtGui.QApplication([""])
	
	window = MainWindow()
	window.show()
	
	app.execute()
	
	
	
if __name__ == '__main__':
	main()
	
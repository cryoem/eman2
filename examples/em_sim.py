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
from emshape import EMShape
from emplot2d import EMPlot2DWidget



class Microscope(QtOpenGL.QGLWidget):
	
	def __init__(self, parent=None, imgwindow=None, pltwindow=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
		self.win_size=[500,1000]
		self.setMinimumSize(self.win_size[0], self.win_size[1])
		self.imgwindow=imgwindow
		self.pltwindow=pltwindow
		
		#### set up the first two lens
		self.source=s=.9
		d0=.15
		d1=.25
		f0=.085
		#f0=(d0*d1-d0*f1)/(d1-f1+d0)
		f1=(d0*d1-f0*d1-f0*d0)/(d0-f0)
		
		### y_pos, focal_pos
		self.lens=[
			### condenser lens 1
			[s-d0, f0],
			
			### condenser lens 2
			[s-d0-d1, f1],
			
			### f==-1 : specimen stage
			[0.102, -1], 
			
			### objective lens
			[-0.03, .061],
			
			### objective aperture
			### for aperture, -4<f<-3. aperture size = f+4 (from 0 to 1)
			[-0.14, -3.75],  
			
			### intermediate lens
			[-0.238, .056],
			#[-0.2, .2],

			
			### intermediate aperture
			[-0.354, -3.01],  
			
			### projector lens
			[-0.462, .06],
			
			### f==-2 : detector
			[-.714, -2] 
			#[-.2, -2] 
			]
		
		self.tmp_pts=[]
		self.drag_lens=-1
		self.hold_shift=False
		self.beam_dist=0
		self.defocus=0
		self.beam_lastxy=[0,0]
		
		#md=[]
		#for t in range(100):
			#self.lens[2][0]=.102009+t*.000001
			
			#df=self.draw_path()
			#print self.lens[2][0], df
			#md.append([self.lens[2][0], df])
		#a=np.argmin(np.array(md)[:,1])
		#print md[a]
			
		
		#print df
		#exit()
		
		if pltwindow:
			pltwindow.show()
		if imgwindow:
			imgwindow.show()
			img=self.draw_wave()
			
			

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
			
		
		self.draw_path()
		glFlush()
		
	def draw_wave(self):
		
		sz=400-1 ### length of x-axis
		xdivz=2e-3 ### x/z scale ratio
		mpix=3e-8 ### pixel to meter
		
		raw=np.zeros(sz) ### input signal
		wid=8 ### width of the slits
		w2=4 ### position of the slits
		
		#raw[sz/2-wid-w2+1:sz/2+wid-w2+1]=1
		#raw[sz/2-wid+w2:sz/2+wid+w2]=1
		
		#raw[sz/2-wid+1:sz/2+wid]=1
		raw0=np.exp(-(np.arange(sz, dtype=float)-sz/2-0)**2/6)
		raw0+=np.exp(-(np.arange(sz, dtype=float)-sz/2-10)**2/3)*.5
		raw=np.exp(-1j*raw0*1)*(np.exp(-raw0*.1))
		wavelen=2e-12/mpix ### wave length
		#print wavelen
		#print raw-raw[::-1]
		#raw=raw*0+1
		#### load lens info 
		mult=400 ### size multiplier from GL lens diagram to wave image
		l=np.array(self.lens)[3:]
		### use relative distance between lens instead of absolute positions
		lens_gl=np.array([[l[i,0]-l[i+1,0], l[i,1]] for i in range(len(l)-1)])
		lens=np.round(mult*lens_gl)
		
		alldata=[]
		for parallel in [False]:
			#### start wave propagation
			imgs=[]
			
			#### from scattering point to the first lens
			pln=np.round((self.lens[2][0]-l[0,0])*mult) ### z position of the first lens
			ix=(np.arange(sz, dtype=float)-sz/2)*xdivz ## indices along x direction
			iz=np.arange(1,pln)[:,None,None] ## indices along z direction
			
			## pln x sz x sz matrix
			dst=np.sqrt((ix-ix[:,None])**2 +iz**2) 
			#print dst
			#print ix.shape, iz.shape, dst.shape, raw.shape, (ix-ix[:,None]).shape
			#pmult=1e-3
			if parallel:
				cpx=np.exp(-1j*2*np.pi*(iz+np.zeros((sz, sz)))/wavelen)*np.mean(raw)#*(1/iz**2)
			else:
				cpx=raw[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)#*(1/dst**2)
			
			
			img=np.sum(cpx, axis=1)
			imgs.append(img)
			shapes=[]
			vz= np.sum(lens[:,0])-len(lens)
			#print img.shape
			for il,ln in enumerate(lens):
				#break
				f=ln[1]
				proj=imgs[-1][-1] ### projection of wave on lens
				if f<0: ### aperture
					proj_ps=proj.copy()
					ap=lens_gl[il][1]+4
					clip=int((1-ap)*sz)/2
					proj_ps[:clip]=0
					proj_ps[-clip:]=0
					#print ln, ap, clip, len(proj_ps[:clip]), len(proj_ps[-clip:])
					shapes.append(EMShape(("rect",.5, .5, 1, 0, vz-2, clip, vz+2, 2)))
					shapes.append(EMShape(("rect",.5, .5, 1, sz-clip, vz-2, sz, vz+2, 2)))
					
				else: ### lens
					ps=((ix)**2)/(f*2)*(2*np.pi/wavelen) ### phase shift
					proj_ps=proj*np.exp(-1j*(-ps)) ### projection after phase shift
					shapes.append(EMShape(("rect",1, .5, .5, 0, vz-2, sz, vz+2, 2)))

				zmax=ln[0]
				iz=np.arange(1,zmax)[:,None,None] ## indices along z direction
				dst=np.sqrt((ix-ix[:,None])**2 +iz**2)
				cpx1=proj_ps[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)#*(1/dst**2)

				img=np.sum(cpx1, axis=1)
				imgs.append(img)
				vz-=zmax-1
			#print [m.shape for m in imgs]
			final=np.vstack(imgs)
			final/=np.sum(abs(final),axis=1)[:,None]
			img=final[::-1,:].copy()
			alldata.append(img)
		#print img.shape
		
		#nrm=(np.sum(abs(alldata[0]), axis=1)+np.sum(abs(alldata[1]), axis=1))/sz
		#for i in [0,1]: alldata[i]/=nrm[:,None]
		
		self.imgwindow.shapes={ i:shapes[i] for i in range(len(shapes)) }
		self.imgwindow.shapechange=1
		#self.imgwindow.set_data([from_numpy(abs(d)*np.sin(np.angle(d))) for d in alldata])
		self.imgwindow.set_data([from_numpy(np.real(d)) for d in alldata])
		
		if self.pltwindow:
			
			a0=alldata[0][0,:]
			#a1=alldata[1][0,:]
			self.pltwindow.set_data([np.arange(sz)-sz/2, abs(a0)/np.max(abs(a0))], "scatter", linetype=0)
			self.pltwindow.set_data([np.arange(sz)-sz/2, raw0], "raw", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, abs(alldata[1][0,:])], "parallel", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, 
			    #abs(a0/abs(a0)+a1/abs(a1))], "contrast", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, np.sin(np.angle(a0))], "phase_scatter", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, np.sin(np.angle(a1))], "phase_parallel", linetype=0)
			
			a0ft=abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift((abs(a0)-np.mean(abs(a0)))))))
			#rawft=np.real(np.fft.fftshift(np.fft.fft(np.fft.fftshift(a0))))
			rawft=np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(raw0-np.mean(raw0)))))
			a0ft/=np.max(a0ft)
			rawft/=np.max(rawft)
			rf=rawft.copy()
			rf[abs(rf)<1e-3]=1
			ctf=a0ft/rf
			ctf[abs(rawft)<1e-3]=0
			
			ftidx=np.fft.fftfreq(len(a0))
			
			#self.pltwindow.set_data([np.arange(sz)-sz/2, a0ft], "scatter_fft", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, rawft], "raw_fft", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, ctf], "ctf", linetype=0)
			#contrast=np.sin(np.angle(a0+a1))
			#contrast=(np.angle(a0)-np.angle(a1))*180/np.pi
			#self.pltwindow.set_data([np.arange(sz)-sz/2, contrast], "phase_contrast", linetype=0)
		return img


	
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
		diff_f=[]
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
				if l[1]==-1 or l[1]==-2 or abs(1./l[1]-1./d0)<1./bign:
					d1=bign
					
				elif l[1]<-3 and l[1]>-4:
					d1=-d0
					
				else:
					d1=1./(1./(l[1])-1./d0)
				
					
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
						if il>0 and s0>self.lens[il-1][0]:
							s1=self.lens[il-1][0]
							x1=w*((s0-s1)/(s0-l[0]))
							self.draw_vertex(sign*x1, s1)
							#print l, self.lens[il-1], s0, x0
						else:
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
				if l[1]==-2: diff_f.append(abs(w))
			
			src=news
			theta=newt
			diverge=newv
			if l[1]>=0:
				scatter=False
		#print ", ".join(["{:.2f}".format(d) for d in dist]), abs(dist[0]*2-dist[1]-dist[2])
		#self.defocus=dist[0]*2-dist[1]-dist[2]
		return np.mean(abs(np.diff(diff_f)))
	
	def draw_lens(self, y=0, focal=.1, scale=1.2):
		
		pts=[]
		
		if focal>=0: ### real lens
			arc=(1./focal)/300.
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
			dy=.04
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
			
		elif focal<=-3 and focal>-4: ### aperture
			asize=1-(focal+4)
			
			dx=scale/2.
			dy=.03
			cc=dx*(.5*asize+.25)
			for lr in [-1, 1]:
				pts=np.array([[dx,dy],[dx-cc,dy],[dx-cc*.9,dy/2.],
						[dx-cc,0], [dx,0], [dx,dy]])
				pts[:,0]*=lr
				pts+=[0,y]
				glColor3f( .8, .8, 1 )
				glEnableClientState(GL_VERTEX_ARRAY)
				glVertexPointerf(pts.tolist())
				glDrawArrays(GL_POLYGON, 0, len(pts))
				glColor3f( .1, .1, 1 )
				glLineWidth(3.)
				glEnableClientState(GL_VERTEX_ARRAY)
				glVertexPointerf(pts.tolist())
				glDrawArrays(GL_LINE_STRIP, 0, len(pts))
			
			
			
		return pts
		
	def scr_to_img(self, pp):
		winsz=np.array(self.win_size, dtype=float)
		p0=np.array([pp.x()-winsz[0]/2., winsz[1]/2.-pp.y()])
		p=p0/winsz*2.
		return p
		
	def mousePressEvent(self, QMouseEvent):
		p=self.scr_to_img(QMouseEvent.pos())
		self.startpy=p[1]
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			self.hold_shift=True
		else:
			self.hold_shift=False
			
		#print p
		self.drag_lens=-1
		for i,l in enumerate(self.lens):
			#print abs(p[1]-l[0])
			if abs(p[1]-l[0])<0.05:
				self.drag_lens=i
				
				break
		
		#if modifiers == QtCore.Qt.ControlModifier:
			#if self.drag_lens>2:
				#print "Removing lens: ", self.drag_lens
				
		if self.drag_lens==-1:
			return
		
		l=self.lens[self.drag_lens][1]
		if l>0 or (l<-3 and l>-4) or l==-1:
			pass
		else:
			self.hold_shift=False
		

	def mouseMoveEvent(self, QMouseEvent):
		if self.drag_lens<0:
			return
		else:
			p=self.scr_to_img(QMouseEvent.pos())
			if self.hold_shift:
				l=self.lens[self.drag_lens][1]
				if l>0:
					### target is a lens
					### adjust strength of lens
					dy=abs(p[1]-self.lens[self.drag_lens][0])
					dy=min(.18, dy)
					dy=0.01/max(.01, dy)
					self.lens[self.drag_lens][1]=dy
					
				elif (l<-3 and l>-4):
					### target is an aperture
					dy=abs(p[1]-self.lens[self.drag_lens][0])
					dy=min(.95, max(.05, dy))
					self.lens[self.drag_lens][1]=-3-dy
					
				elif (l==-1):
					dy=p[1]-self.startpy
					self.lens[self.drag_lens][0]+=dy*.002
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
			print "lens {:d}: py={:.3f}, f={:.3f}".format(
				self.drag_lens, l[0], l[1])
			self.drag_lens=-1
			
			if self.imgwindow:
				img=self.draw_wave()
	
	def closeEvent(self, event):
		print self.lens
		print "Exit.."
		exit()

class MainWindow(QtGui.QMainWindow):
	
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		
		self.imgview = EMImage2DWidget()
		self.pltview = EMPlot2DWidget()
		#img=np.random.rand(512,512)
		#self.imgview.set_data(from_numpy(img))
		
		widget = Microscope(self, self.imgview, self.pltview)
		self.closeEvent=widget.closeEvent
		self.setCentralWidget(widget)

		
		

def main():
	app = EMApp()
#	app=QtGui.QApplication([""])
	
	window = MainWindow()
	window.show()
	
	app.execute()
	
	
	
if __name__ == '__main__':
	main()
	

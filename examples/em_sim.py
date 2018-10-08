#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# Muyuan Chen 2017-08
from past.utils import old_div
from builtins import range
import numpy as np

from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtGui, QtCore, QtOpenGL
from eman2_gui.emimage2d import EMImage2DWidget
from EMAN2 import *
from eman2_gui.emapplication import EMApp
from eman2_gui.emshape import EMShape
from eman2_gui.emplot2d import EMPlot2DWidget



class Microscope(QtOpenGL.QGLWidget):
	
	
	def lens_sets(self, mode=0):
		
		#### C1 and C2
		#### ensure parallel beam from condenser lens
		s=self.source=.9
		d0=.15
		d1=.25
		f0=.085
		#f0=(d0*d1-d0*f1)/(d1-f1+d0)
		f1=old_div((d0*d1-f0*d1-f0*d0),(d0-f0))
		self.lens=[
			[s-d0, f0], ### condenser lens 1
			[s-d0-d1, f1], ### condenser lens 2
			[0., -1], ### f==-1 : specimen stage
			]
		
		#### imaging mode with all lens
		if mode==0:
			self.lens.extend([
				[-0.1, .06], ### objective lens
				#### for aperture, -4<f<-3. aperture size = f+4 (from 0 to 1)
				[-0.16, -3.75],  #### objective aperture
				[-0.35, .06], #### intermediate lens
				[-0.59, .065], #### projector lens
				[-.824, -2], ### f==-2 : detector 
			])
		
		#### diffraction mode
		elif mode==1:
			self.lens.extend([
				[-0.1, .06], ### objective lens
				[-0.264, -3.75], #### objective aperture
				[-0.326, .0756], #### intermediate lens
				[-0.552, .06667], #### projector lens
				[-.824, -2], ### f==-2 : detector 
				])
		
		elif mode==2:
			self.lens.extend([
				[-0.09, old_div(1.,15)], ### objective lens
				[-0.09-old_div(1.,15), -5], #### phase plate
				[-.347, -2], ### f==-2 : detector 
			])
			
		#### imaging mode with a single lens
		else:# mode==2:
			self.lens.extend([
				[-0.1, old_div(1.,15)], ### objective lens
				[-.3, -2], ### f==-2 : detector 
			])
			
				
	
	def __init__(self, parent=None, options=None, imgwindow=None, pltwindow=None, twodwindow=None):
		QtOpenGL.QGLWidget.__init__(self, parent)
		
		self.win_size=[500,1000]
		self.setMinimumSize(self.win_size[0], self.win_size[1])
		self.imgwindow=imgwindow
		self.pltwindow=pltwindow
		self.twodwindow=twodwindow
		self.mag=1
		self.lens_sets(options.mode)
		self.mode=options.mode
		self.cs=options.cs
		self.options=options
		
		
		self.wavesign=1
		
		
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
			self.draw_path()
			imgwindow.show()
			img=self.draw_wave()
		if twodwindow:
			twodwindow.show()
			img = self.draw_wave_twod()
			
			

	def initializeGL(self):
		glClearColor(1.0, 1.0, 1.0, 1.0)
		glClearDepth(1.0)
		glLineStipple(1, 0xAAAA);

		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(30.0, 1.0, 1.0, 30.0)
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
		gluPerspective(30.0, 1.0, 1.0, 30.0)
		self.win_size=[w,h]

	def paintGL(self):

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glLoadIdentity()
		
		for l in self.lens:
			self.draw_lens(l[0], l[1]);
			
		
		self.draw_path()
		glFlush()
		
	def gauss_edge(self, sz, cl=.4, gw=.2, gs=.06):
		raw=np.ones(sz, dtype=np.complex)
		clip=int(sz*cl)
		gwidth=int(sz*gw)
		gwidth=min(clip-1, gwidth)
		gsig=sz*gs
		raw[:clip]=0
		raw[-clip:]=0
		gaus=np.exp(old_div(-np.arange(gwidth, dtype=float)**2,gsig**2))
		raw[clip-gwidth:clip]=gaus[::-1]
		raw[-clip:-clip+gwidth]=gaus
		return raw
	
	def draw_wave(self):
		
		sz=512-1 ### length of x-axis
		xdivz=2e-3 ### x/z scale ratio
		mpix=4e-8 ### pixel to meter
		wavelen=2e-12/mpix ### wave length
		wid=8 ### width of the slits
		w2=4 ### position of the slits
		
		#raw=np.zeros(sz) ### input signal
		#raw[sz/2-wid-w2+1:sz/2+wid-w2+1]=1
		#raw[sz/2-wid+w2:sz/2+wid+w2]=1
		
		#raw[sz/2-wid+1:sz/2+wid]=1
		raw0=np.zeros(sz) ### input signal
		ni=self.options.ninput
		dx=30/ni
		for i in range(ni):
			x=dx*i-15*(ni>1)
			raw0+=np.exp(-(np.arange(sz, dtype="float32")-(sz/2)+x)**2/5)
		
		#raw0+=np.exp(-(np.arange(sz, dtype=float)-sz/2-7.5)**2/5)
		#raw0+=np.exp(-(np.arange(sz, dtype=float)-sz/2+7.5)**2/5)
		#raw=np.exp(-1j*raw0*1)#*(np.exp(-raw0*.1))
		#raw=raw0.copy()
		raw=self.gauss_edge(sz, cl=.2, gw=.1)
		raw*=np.exp(-1j*raw0*1)
		#self.pltwindow.set_data([np.arange(sz)-sz/2, raw], "raw", linetype=0)
		#return
		#print wavelen
		#print raw-raw[::-1]
		#raw=raw*0+1
		#### load lens info 
		mult=400 ### size multiplier from GL lens diagram to wave image
		l=np.array(self.lens)[3:]
		### use relative distance between lens instead of absolute positions
		lens_gl=np.array([[l[i,0]-l[i+1,0], l[i,1]] for i in range(len(l)-1)])
		lens=mult*lens_gl
		#lens[:,0]=np.round(lens[:,0]) ### round z position
		#lens[-1][0]+=10
		
		ix=(np.arange(sz, dtype=float)-old_div(sz,2))*xdivz ## indices along x direction
		ix_mat=(ix-ix[:,None])**2
		cs=self.cs * mult ## spherical abberation
		alldata=[]
		for parallel in [False]: ### compute parallel beam when true
			#### start wave propagation
			imgs=[]
			
			#### from scattering point to the first lens
			zmax=np.round((self.lens[2][0]-l[0,0])*mult) ### z position of the first lens
			#iz=np.arange(1,zmax)[:,None,None] ## indices along z direction
			iz=(np.arange(1, int(zmax)+1, dtype=float)/int(zmax)*zmax)[:,None,None]
			## zmax x sz x sz matrix
			dst=np.sqrt(ix_mat +iz**2) 
			#print dst
			#print ix.shape, iz.shape, dst.shape, raw.shape, (ix-ix[:,None]).shape
			#pmult=1e-3
			if parallel:
				cpx=np.exp(-1j*2*np.pi*(iz+np.zeros((sz, sz)))/wavelen)*np.mean(raw)*(old_div(1,iz**2))
			else:
				cpx=raw[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)*(old_div(1,dst**2))
			
			
			img=np.sum(cpx, axis=1)
			imgs.append(img)
			shapes=[]
			vz= np.sum(lens[:,0])-len(lens) ### to track z position
			#print img.shape
			for il,ln in enumerate(lens):
				#break
				f=ln[1]
				proj=imgs[-1][-1] ### projection of wave on lens
				fv=lens_gl[il][1]
				
				if (fv<-3 and fv>-4): ### aperture
					ap=fv+4
					clip=old_div(int((1-ap)*sz),2)
					msk=self.gauss_edge(sz, cl=old_div((1-ap),2.))
					proj_ps=proj.copy()*msk
					shapes.append(EMShape(("rect",.5, .5, 1, 0, vz-2, clip, vz+2, 2)))
					shapes.append(EMShape(("rect",.5, .5, 1, sz-clip, vz-2, sz, vz+2, 2)))
					
				elif fv==-5: ### phase plate
					msk=self.gauss_edge(sz, cl=old_div(.1,2.))
					phaseplate=self.gauss_edge(sz, cl=.49, gw=.1, gs=.01)
					phaseplate=-phaseplate*np.pi/2.
					proj_ps=proj.copy()*msk*np.exp(-1j*phaseplate)
					
					shapes.append(EMShape(("rect",.8, .5, 1, 0, vz, sz, vz, 2)))
					
				else: ### lens
					ps=old_div(((ix)**2),(f*2))*(2*np.pi/wavelen) ### phase shift
					ps+=cs*(ix**4)/4.
					proj_ps=proj*np.exp(-1j*(-ps)) ### projection after phase 
					shapes.append(EMShape(("rect",1, .5, .5, 0, vz-2, sz, vz+2, 2)))

				zmax=ln[0]
				#iz=np.arange(1,zmax, dtype=float)[:,None,None] ## indices along z direction
				iz=(np.arange(1, int(zmax)+1, dtype=float)/int(zmax)*zmax)[:,None,None]
				dst=np.sqrt(ix_mat +iz**2)
				cpx1=proj_ps[:,None]*np.exp(-1j*2*np.pi*dst/wavelen)*(old_div(1,dst**2))

				img=np.sum(cpx1, axis=1)
				imgs.append(img)
				vz-=zmax-1
			#print [m.shape for m in imgs]
			final=np.vstack(imgs)
			final/=np.sum(abs(final),axis=1)[:,None]
			img=final[::-1,:].copy()
			#img/=np.max(img)
			alldata.append(img)
		#print img.shape
		
		#nrm=(np.sum(abs(alldata[0]), axis=1)+np.sum(abs(alldata[1]), axis=1))/sz
		#for i in [0,1]: alldata[i]/=nrm[:,None]
		
		self.imgwindow.shapes={ i:shapes[i] for i in range(len(shapes)) }
		self.imgwindow.shapechange=1
		#self.imgwindow.set_data([from_numpy(abs(d)*np.sin(np.angle(d))) for d in alldata])
		self.imgwindow.set_data([
			from_numpy(abs(alldata[0]).astype("float32")), 
			from_numpy(np.cos(np.angle(alldata[0])).astype("float32")) ])
		
		if self.pltwindow:
			
			a0=alldata[0][0,:]
			#a1=alldata[1][0,:]
			r0=(np.arange(sz)-old_div(sz,2))*abs(self.mag)
			if self.wavesign<0: r0=r0[::-1]
			clip=np.where(abs(r0)<old_div(sz,2))[0]
			if len(clip)>0:
				c=clip[0]
				rawplot=[r0[c:-c], raw0[c:-c]]
			else:
				rawplot=[r0, raw0]
			
			#rawplot=[np.arange(sz)-sz/2, raw0]
			self.pltwindow.set_data([np.arange(sz)-old_div(sz,2), old_div(abs(a0),np.max(abs(a0)))], "scatter", linetype=0)
			self.pltwindow.set_data(rawplot, "raw", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, abs(alldata[1][0,:])], "parallel", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, 
			    #abs(a0/abs(a0)+a1/abs(a1))], "contrast", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, np.sin(np.angle(a0))], "phase_scatter", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, np.sin(np.angle(a1))], "phase_parallel", linetype=0)
			aa=abs(a0)
			if self.wavesign>0:
				bd=np.where(raw0>.01)[0]
			else:
				bd=np.where(raw0[::-1]>.01)[0]
			bd=[bd[0], bd[-1]]
			pad=int((bd[1]-bd[0])*1.)
			bd=[bd[0]-pad, bd[1]+pad]
			bd=[int((b-old_div(sz,2))*abs(self.mag)+old_div(sz,2)) for b in bd]
			
			bd[0]=max(0,bd[0])
			bd[1]=min(sz-1, bd[1])
			#print bd
			aa-=old_div((aa[bd[0]]+aa[bd[1]]),2.)
			#print bd
			aa[:bd[0]]=0
			aa[bd[1]-1:]=0
			aa/=np.max(abs(aa))
			a0ft=np.real(np.fft.fftshift(np.fft.fft(np.fft.fftshift((aa)))))
			#rawft=np.real(np.fft.fftshift(np.fft.fft(np.fft.fftshift(a0))))
			rawft=np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(-raw0))))
			a0ft/=np.max(a0ft)
			#a0ft[sz/2]=1
			rawft/=np.max(rawft)
			rf=rawft.copy()
			rf[abs(rf)<1e-3]=1
			ctf=old_div(a0ft,rf)
			ctf[abs(rawft)<1e-3]=0
			
			ftidx=np.fft.fftfreq(len(a0))
			#print self.mag
			#self.pltwindow.set_data([np.arange(sz)-sz/2, aa], "scatter_msk", linetype=0)
			#self.pltwindow.set_data([r0[c:-c], a0ft[c:-c]], "scatter_fft", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, rawft], "raw_fft", linetype=0)
			#self.pltwindow.set_data([np.arange(sz)-sz/2, ctf], "ctf", linetype=0)
			#contrast=np.sin(np.angle(a0+a1))
			#contrast=(np.angle(a0)-np.angle(a1))*180/np.pi
			#self.pltwindow.set_data([np.arange(sz)-sz/2, contrast], "phase_contrast", linetype=0)
		return img

	def specimen_propagation(self):
		sz = self.sz
		xydivz=2e-3 ### xy/z scale ratio
		mpix=3e-8 ### pixel to meter
		wavelen=old_div(2e-12,mpix) ### wave length
		mult=floor(sz*0.5) ### size multiplier from GL lens diagram to wave image

		yap,xap = np.ogrid[old_div(-sz,2):old_div(sz,2), old_div(-sz,2):old_div(sz,2)]

		#### from top to bottom of specimen
		spec = []

		layers = []
		layers.append( [[ -10,  -10, 0.15, 4., 8.]] ) # t = 3-0 = 3
		layers.append( [[  0,  0,  0.3, 8., 8.]] ) # t = 3-1 = 2
		layers.append( [[  10, 10, 0.45, 4., 8.]] ) # t = 3-2 = 1

		scale = 5
		t = len(layers)*scale # specimen thickness (in pixels)

		raw1 = np.zeros((sz,sz))

		ix2d, iy2d=(np.array(np.indices((sz,sz)), dtype=np.float32)-old_div(sz,2))*xydivz
		rr=ix2d**2+iy2d**2
		ixy = (ix2d[:,:,None,None]-ix2d[None,None,:])**2+ (iy2d[:,:,None,None]-iy2d[None,None,:])**2
		for t0,coords in enumerate(layers):
			for x0,y0,s,sigx,sigy in coords:
				x0 = float(x0)
				y0 = float(y0)
				s = float(s)
				sigx = float(sigx)
				sigy = float(sigy)
				raw1 += s*np.exp(-(old_div((xap-x0)**2,sigx)+old_div((yap-y0)**2,sigy)),dtype=np.float32)
			raw1 /= len(coords)
			raw1 /= np.max(raw1)
			phaseobj = np.exp(-1j*raw1*2*np.pi)

			z = t0*scale
			dst = np.sqrt(ixy+(t-z)**2)
			cpx=phaseobj[:,:,None,None]*np.exp(-1j*2.*np.pi*dst/wavelen)*(old_div(1,dst**2))
			img = np.sum(cpx, axis=(1,0))
			img/=np.max(abs(img))
			spec.append(img)

		return np.sum(spec,axis=0)

	def draw_wave_twod(self):
		self.sz = 64
		sz=self.sz ### length of x-axis

		xydivz=2e-3 ### xy/z scale ratio
		mpix=3e-8 ### pixel to meter
		wavelen=old_div(2e-12,mpix) ### wave length
		mult=floor(sz*0.78) ### size multiplier from GL lens diagram to wave image

		#### load lens info
		l=np.array(self.lens)[3:]
		lens_gl=np.array([[l[i,0]-l[i+1,0], l[i,1]] for i in range(len(l)-1)]) ### use relative distance between lens instead of absolute positions
		lens=np.round(mult*lens_gl)

		#### start wave propagation
		imgs=[]

		### specimen exit waveform
		raw = self.specimen_propagation()
		imgs.append(raw)
		
		#### from exiting specimen to the first lens
		yap,xap = np.ogrid[old_div(-sz,2):old_div(sz,2), old_div(-sz,2):old_div(sz,2)]
		iz=np.round((self.lens[2][0]-l[0,0])*mult) ### z position of the first lens
		ix2d, iy2d=(np.array(np.indices((sz,sz)), dtype=np.float32)-old_div(sz,2))*xydivz
		rr=ix2d**2+iy2d**2
		ixy = (ix2d[:,:,None,None]-ix2d[None,None,:])**2+ (iy2d[:,:,None,None]-iy2d[None,None,:])**2
		dst = np.sqrt(ixy+iz**2)
		cpx=raw[:,:,None,None]*np.exp(-1j*2.*np.pi*dst/wavelen)*(old_div(1,dst**2))
		img = np.sum(cpx, axis=(1,0))
		img/=np.max(abs(img))
		imgs.append(img)
		vz= np.sum(lens[:,0])-len(lens)
		for il,ln in enumerate(lens):
			zmax=ln[0]
			f=ln[1]
			proj=imgs[-1] ### projection of wave on lens
			if f<0: ### aperture
				proj_ps=proj.copy()
				ap=lens_gl[il][1]+4
				clip=old_div(int((1-ap)*sz),2)
				proj_ps *= (xap*xap+yap*yap<clip**2)
			else: ### lens
				ps=old_div(rr,(f*2))*(2*np.pi/wavelen)  ### phase shift
				ps+=self.cs*(rr**4)/4.
				proj_ps=proj*np.exp(-1j*(-ps)) ### projection after phase shift
			dst = np.sqrt(ixy+zmax**2)
			cpx1=proj_ps[:,:,None,None]*np.exp(-1j*2*np.pi*dst/wavelen)*(old_div(1,dst**2))
			img=np.sum(cpx1,axis=(1,0))
			img/=np.max(abs(img))
			imgs.append(img)
			vz -= zmax-1

		self.twodwindow.set_data([from_numpy(np.real(d).copy().astype("float32")) for d in imgs])
		return imgs[-1]
	
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
		self.wavesign=1
		src=[self.source] ## input position
		theta=[np.tan(old_div(np.pi,3.))] ## input angle
		diverge=[False for s in src] ### whether the beam at previous layer diverges
		
		bign=1e10 ### a big number. so we do not reach inf when the beam is paralleled
		
		
		beamstop=False ### stop the beam at some point..
		scatter=False ### scattered beam from last stage
		ymax=self.lens[-1][0]
		
		self.beam_dist=0
		dist=[0] ### track beam distance
		passsample=False
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
				if l[1]==-1 or l[1]==-2 or abs(old_div(1.,l[1])-old_div(1.,d0))<old_div(1.,bign):
					d1=bign
					
				elif l[1]<-3 and l[1]>-4:
					d1=-d0
					
				else:
					d1=old_div(1.,(old_div(1.,(l[1]))-old_div(1.,d0)))
				
					
				if abs(d0)>bign:
					d0=np.sign(d0)*bign
					
				w=tt*d0 ### width of beam when it hit the lens
				tt=old_div(w,d1) ### output angle
				
				
				#### deal with beam scattering
				if not scatter:
					s0=s
					x0=0
				else:
					## beam cannot be scattered at the first stage.
					x0=s-self.lens[il-1][0]
					s0=self.lens[il-1][0]
					#print si, x0, s0, d0
				
				ym=s-old_div(d0,abs(w))
				xm=1
				if ym<ymax:
					xm=old_div(xm,ym)*ymax
					ym=ymax

				glColor3f( .5, 1., .2 + (si>0)*.7 );
				glBegin(GL_LINES)
				if d1>0: ### output beam converge to next focal point
					
					if (scatter and abs(x0)<bign) or (d0>0 and dvg==False): 
						#### lens below last focal point, draw input beam
						if si==0: 
							self.wavesign*=-1
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
						mult=old_div((l1-l[0]),d1)-1
						if w*(2+mult)>maxwidth:
							#mult=abs(1./w)
							mult=old_div((ymax-l[0]),d1)-1
							#print mult
							beamstop=True
					else:
						mult=abs(old_div(1.,w))
					
					
					for sign in [-1,1]:
						if il>0 and s0>self.lens[il-1][0]:
							s1=self.lens[il-1][0]
							x1=w*(old_div((s0-s1),(s0-l[0])))
							self.draw_vertex(sign*x1, s1)
							#print l, self.lens[il-1], s0, x0
						else:
							if si==0 and x0==0 and sign==-1: 
								self.wavesign*=-1
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
					if passsample:
						print("Multiple sample stage??")
					diff_f.append(abs(w))
					news.append(s)
					newt.append(tt)
					newv.append(False)
					
					news.append(l[0]+abs(w))
					newt.append(np.tan(old_div(np.pi,4.)))
					newv.append(False)
					dist.append(dist[si])
					
					news.append(l[0]-abs(w))
					newt.append(np.tan(old_div(np.pi,4.)))
					newv.append(False)
					dist.append(dist[si])
					scatter=True
					passsample=True
				
				else:
					if not beamstop:
						news.append(s)
						newt.append(tt)
						newv.append(dvg)
					else:
						news.append(None)
						newt.append(0)
						newv.append(False)
				
				#if passsample and si==0: 
					#print il, self.wavesign
					#self.wavesign*=np.sign(w)*-1
				if l[1]==-2 and si==0: diff_f.append(w)
			
			src=news
			theta=newt
			diverge=newv
			if l[1]>=0:
				scatter=False
		#print ", ".join(["{:.2f}".format(d) for d in dist]), abs(dist[0]*2-dist[1]-dist[2])
		#self.defocus=dist[0]*2-dist[1]-dist[2]
		self.mag= old_div(diff_f[1],diff_f[0])
		#print diff_f, self.mag, self.wavesign
		return self.mag
	
	def draw_lens(self, y=0, focal=.1, scale=1.2):
		
		pts=[]
		
		if focal>=0: ### real lens
			arc=old_div((old_div(1.,focal)),400.)
			t=np.arange(-np.pi*arc, np.pi*arc+1e-5, np.pi*arc/20.)+old_div(np.pi,2.)
			p0=np.vstack([np.cos(t), np.sin(t)]).T
			p0[:,1]-=p0[0,1]
			l=p0[0,0]-p0[-1,0]
			n=len(p0)
			p1=np.vstack([np.cos(t[::-1]), -np.sin(t[::-1])]).T
			p1[:,1]-=p1[0,1]
			pts=np.vstack([p0,p1])
			pts=old_div(pts,l)*scale
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
			dx=old_div(scale,2.)
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
			dx=old_div(scale,2.)
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
			
		elif focal<=-3 and focal>-4: ### aperture
			asize=1-(focal+4)
			
			dx=old_div(scale,2.)
			dy=.03
			cc=dx*(.5*asize+.25)
			for lr in [-1, 1]:
				pts=np.array([[dx,dy],[dx-cc,dy],[dx-cc*.9,old_div(dy,2.)],
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
			
		
		elif focal==-5: ### phase plate
			
			dx=old_div(scale,2.)
			
			for lr in [-1, 1]:
				dy=.01
				cc=dx*.95
				pts=np.array([[dx,old_div(dy,2.)],[dx-cc,old_div(dy,2.)],
						[dx-cc,old_div(-dy,2.)], [dx,old_div(-dy,2.)], [dx,old_div(dy,2.)]])
				pts[:,0]*=lr
				pts+=[0,y]
				glColor3f( .9, .8, 1 )
				glEnableClientState(GL_VERTEX_ARRAY)
				glVertexPointerf(pts.tolist())
				glDrawArrays(GL_POLYGON, 0, len(pts))
				#####
				dy=.03
				cc=dx*.3
				pts=np.array([[dx,old_div(dy,2.)],[dx-cc,old_div(dy,2.)],
						[dx-cc,old_div(-dy,2.)], [dx,old_div(-dy,2.)], [dx,old_div(dy,2.)]])
				pts[:,0]*=lr
				pts+=[0,y]
				glColor3f( .9, .8, 1 )
				glEnableClientState(GL_VERTEX_ARRAY)
				glVertexPointerf(pts.tolist())
				glDrawArrays(GL_POLYGON, 0, len(pts))
				
				glColor3f( .5, .1, 1 )
				glLineWidth(3.)
				glEnableClientState(GL_VERTEX_ARRAY)
				glVertexPointerf(pts.tolist())
				glDrawArrays(GL_LINE_STRIP, 0, len(pts))
			
		return pts
		
	def scr_to_img(self, pp):
		winsz=np.array(self.win_size, dtype=float)
		p0=np.array([pp.x()-old_div(winsz[0],2.), old_div(winsz[1],2.)-pp.y()])
		p=old_div(p0,winsz)*2.
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
					dy=min(.15, dy)
					dy=old_div(0.01,max(.01, dy))
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
				f=old_div((d0*d1-f0*d1-f0*d0),(d0-f0))
				self.lens[1][1]=max(0.02, f)
			elif self.drag_lens==1:
				f=old_div((d0*d1-d0*f1),(d1-f1+d0))
				self.lens[0][1]=max(0.02, f)
				
			self.updateGL()
		
	
	
	def mouseReleaseEvent(self, QMouseEvent):
		if self.drag_lens>=0:
			l=self.lens[self.drag_lens]
			if self.imgwindow:
				img=self.draw_wave()
			if self.twodwindow:
				img=self.draw_wave_twod()
			print("lens {:d}: py={:.3f}, f={:.3f}".format(
				self.drag_lens, l[0], l[1]))
			self.drag_lens=-1
			
	
	def closeEvent(self, event):
		print(self.lens)
		print("Exit..")
		exit()

class MainWindow(QtGui.QMainWindow):
	
	def __init__(self, options):
		QtGui.QMainWindow.__init__(self)
		
		if options.twod:
			self.twodview = EMImage2DWidget()
			widget = Microscope(self, options, None, None, self.twodview)
		else:
			self.imgview = EMImage2DWidget()
			self.pltview = EMPlot2DWidget()
			widget = Microscope(self, options, self.imgview, self.pltview, None)
			#widget = Microscope(self, None, None, None, mode, cs)

		#widget = Microscope(self, None, None)
		self.closeEvent=widget.closeEvent
		self.setCentralWidget(widget)

		
		

def main():
	
	
	usage="Electron Microscope simulation"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--mode", type=int,help="operating mode. 0: imaging mode with all lens; 1:diffraction mode with all lens; 2: imaging mode with single lens and phase plate; else: imaging mode with one lens.", default=0)
	parser.add_argument("--cs", type=float,help="Cs of microscope.", default=0.0)
	parser.add_argument("--ninput", type=int,help="number of peaks in sampel.", default=4)
	parser.add_argument("--twod", action="store_true", help="Show twod image instead of 1D plot and beam propagation in the microscope.", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	
	app = EMApp()
#	app=QtGui.QApplication([""])
	
	window = MainWindow(options)
	window.show()
	window.raise_()
	
	app.execute()
	
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	

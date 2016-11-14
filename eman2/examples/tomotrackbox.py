#!/usr/bin/env python

#
# Author: Steven Ludtke 10/05/2009 (sludtke@bcm.edu)
# Copyright (c) 2000-2009 Baylor College of Medicine
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

from optparse import OptionParser
import sys
import os
from EMAN2 import *
from emapplication import EMApp
from emimage2d import EMImage2DWidget
from emimage3d import EMImage3DModule
from valslider import ValSlider
import weakref
from emshape import EMShape
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
#import EMAN2db

reconmodes=["gauss_2","gauss_3","gauss_5"]

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input ali>
	
This program will allow the user to select a specific feature within a tomogram and track it, boxing out the
feature from all slices. Generally best for uniform objects like vesicles."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--prefix",type="string",help="Output file prefix",default=None)
	parser.add_option("--tiltstep",type="float",help="Step between tilts in degrees, default=2",default=2.0)
	parser.add_option("--maxshift",type="int",help="Maximum shift in autoalign, default=8",default=8)
	parser.add_option("--seqali",action="store_true",default=False,help="Average particles in previous slices for better alignments of near-spherical objects")
	parser.add_option("--invert",action="store_true",default=False,help="Inverts image data contrast on the fly")

	#parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=64)
	#parser.add_option("--writeoutput",action="store_true",default=False,help="Uses coordinates stored in the tomoboxer database to write output")
	#parser.add_option("--stack",action="store_true",default=False,help="Causes the output images to be written to a stack")
	#parser.add_option("--force",action="store_true",default=False,help="Force overwrite output")
	#parser.add_option("--normproc", help="Normalization processor to apply to particle images. Should be normalize, normalize.edgemean or None", default="normalize.edgemean")
	#parser.add_option("--invertoutput",action="store_true",help="If writing output only, this will invert the pixel intensities of the boxed files",default=False)
	#parser.add_option("--dbls",type="string",help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	#parser.add_option("--outformat", help="Format of the output particles images, should be bdb,img, or hdf", default="bdb")
	
	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)
	
	app = EMApp()

	
	track=TrackerControl(app,options.maxshift,options.invert,options.seqali,options.tiltstep)
	track.set_image(args[0])

	app.execute()
	
	E2end(logid)
	
class TrackerControl(QtGui.QWidget):
	def __init__(self,app,maxshift,invert=False,seqali=False,tiltstep=2.0):
		self.app=app
		self.maxshift=maxshift
		self.seqali=seqali
		self.invert=invert
		self.tiltstep=tiltstep
		
		# the control panel
		QtGui.QWidget.__init__(self,None)

		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(0)
		self.gbl.setSpacing(6)
		self.gbl.setObjectName("hbl")
		
		# action buttons
		self.bcenalign=QtGui.QPushButton("Center Align")
		self.bprojalign=QtGui.QPushButton("Proj. Realign")
		self.btiltaxis=QtGui.QPushButton("Tilt Axis")
		self.btiltaxisval=QtGui.QLineEdit("90.0")
		self.bsavedata=QtGui.QPushButton("Save Data")
		self.breconst=QtGui.QPushButton("3D Normal")
		self.sbmode=QtGui.QSpinBox(self)
		self.sbmode.setRange(0,2)
		self.sbmode.setValue(0)
		self.bmagict=QtGui.QPushButton("3D Tomofill")
		self.bmagics=QtGui.QPushButton("3D Sph")
		self.bmagicc=QtGui.QPushButton("3D Cyl")
		self.vslpfilt=ValSlider(self,(0,.5),"Filter",0.5,50)
		
		self.gbl.addWidget(self.bcenalign,0,0)
		self.gbl.addWidget(self.bprojalign,0,1)
		self.gbl.addWidget(self.btiltaxis,0,2)
		self.gbl.addWidget(self.btiltaxisval,0,3)
#		self.gbl.addWidget(self.bsavedata,0,3)
		self.gbl.addWidget(self.breconst,1,0)
		self.gbl.addWidget(self.sbmode,2,0,1,1)
		self.gbl.addWidget(self.vslpfilt,3,0,1,4)
		self.gbl.addWidget(self.bmagict,1,1)
		self.gbl.addWidget(self.bmagics,1,2)
		self.gbl.addWidget(self.bmagicc,1,3)
		
		QtCore.QObject.connect(self.bcenalign,QtCore.SIGNAL("clicked(bool)"),self.do_cenalign)
		QtCore.QObject.connect(self.bprojalign,QtCore.SIGNAL("clicked(bool)"),self.do_projalign)
		QtCore.QObject.connect(self.btiltaxis,QtCore.SIGNAL("clicked(bool)"),self.do_tiltaxis)
		QtCore.QObject.connect(self.bsavedata,QtCore.SIGNAL("clicked(bool)"),self.do_savedata)
		QtCore.QObject.connect(self.breconst,QtCore.SIGNAL("clicked(bool)"),self.do_reconst)
		QtCore.QObject.connect(self.bmagict,QtCore.SIGNAL("clicked(bool)"),self.do_magict)
		QtCore.QObject.connect(self.bmagics,QtCore.SIGNAL("clicked(bool)"),self.do_magics)
		QtCore.QObject.connect(self.bmagicc,QtCore.SIGNAL("clicked(bool)"),self.do_magicc)
		QtCore.QObject.connect(self.vslpfilt,QtCore.SIGNAL("valueChanged"),self.do_filter)

		# the single image display widget
		self.im2d =    EMImage2DWidget(application=app,winid="tomotrackbox.big")
		self.imboxed = EMImage2DWidget(application=app,winid="tomotrackbox.small")
		self.improj =  EMImage2DWidget(application=app,winid="tomotrackbox.proj")
		self.imslice = EMImage2DWidget(application=app,winid="tomotrackbox.3dslice")
		self.imvol =   EMImage3DModule(application=app,winid="tomotrackbox.3d")
	
		# get some signals from the window. 
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("mousedown"),self.down)
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("mousedrag"),self.drag)
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("mouseup"),self.up)
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("increment_list_data"),self.change_tilt)
	
		self.imagefile=None
		self.imageparm=None
		self.tiltshapes=None
		self.curtilt=0
		self.oldtilt=self.curtilt
		self.map3d=None
		self.downloc=None
		self.downadjloc=None
		
		self.show()
		self.im2d.show()
		
	def closeEvent(self,event):
		self.im2d.closeEvent(QtGui.QCloseEvent())
		self.imboxed.closeEvent(QtGui.QCloseEvent())
		self.improj.closeEvent(QtGui.QCloseEvent())
		self.imslice.closeEvent(QtGui.QCloseEvent())
		self.imvol.closeEvent(QtGui.QCloseEvent())
		event.accept()
		
	def do_cenalign(self,x=0):
		"""In response to the center align button. Just a wrapper"""
		self.cenalign_stack()
		self.update_stack()
		
	def do_projalign(self,x=0):
		"""In response to the projection align button. Just a wrapper"""
		self.projection_align(self.tiltstep)
		self.update_stack()
#		self.do_reconst()

	def do_tiltaxis(self):
		"""In response to the tilt axis button. Just a wrapper"""
		self.tilt_axis()
		
	def do_reconst(self,x=0):
		"""In response to the normal reconstruction button. Just a wrapper"""
		stack=self.get_boxed_stack()
		mode=self.sbmode.value()
		self.map3d=self.reconstruct(stack,self.tiltstep,mode)
		self.update_3d()

		
	def do_magict(self,x):
		"""In response to tomographic filling reconstruction button. Just a wrapper"""
		stack=self.get_boxed_stack()
#		self.map3d=self.reconstruct_ca(stack[5:-4],0.5)
#		init=self.reconstruct_ca(stack[5:-4],0.5)
		mode=self.sbmode.value()
		self.map3d=self.reconstruct_wedgefill(stack,self.tiltstep,mode)
		self.update_3d()
		
	def do_magics(self,x):
		"""In response to the 3D Sph button. Just a wrapper"""
		return
		
	def do_magicc(self,x):
		"""In response to the 3D cyl button. Just a wrapper"""
		return
		
	def do_filter(self,v):
		"""In response to the filter ValSlider"""
		if self.map3d==None : return
		self.lpfilt=v
		self.update_3d()
		
	def do_savedata(self):
		""
	
	def update_3d(self):
		if self.map3d==None : return
		
		self.filt3d=self.map3d.process("filter.lowpass.gauss",{"cutoff_abs":self.vslpfilt.getValue()})
		
		self.imvol.set_data(self.filt3d)
		self.imvol.show()
		self.imvol.updateGL()

		sz=self.map3d["nx"]
		xsum=self.filt3d.process("misc.directional_sum",{"axis":"x"})
		xsum.set_size(sz,sz,1)
		ysum=self.filt3d.process("misc.directional_sum",{"axis":"y"})
		ysum.set_size(sz,sz,1)
		zsum=self.filt3d.process("misc.directional_sum",{"axis":"z"})
		zsum.set_size(sz,sz,1)
		
		self.improj.set_data([zsum,ysum,xsum])
		self.improj.show()
		self.improj.updateGL()

		self.imslice.set_data(self.filt3d)
		self.imslice.show()
		self.imslice.updateGL()
	
	def update_stack(self):
		stack=self.get_boxed_stack()
		self.imboxed.set_data(stack)
		self.imboxed.show()
		self.imboxed.updateGL()
		

	def set_image(self,fsp):
		"""Takes an ali file to process"""
		self.imageparm=EMData(fsp,0,True).get_attr_dict()
		print "%d slices at %d x %d"%(self.imageparm["nz"],self.imageparm["nx"],self.imageparm["ny"])

		self.imagefile=fsp

		self.curtilt=self.imageparm["nz"]/2
		self.tiltshapes=[None for i in range(self.imageparm["nz"])]
		self.update_tilt()

	def update_tilt(self):
		if self.imagefile==None : return
		
		self.curimg=EMData(self.imagefile,0,False,Region(0,0,self.curtilt,self.imageparm["nx"],self.imageparm["ny"],1))
		if self.invert : self.curimg.mult(-1.0)
		self.im2d.set_data(self.curimg)

		s=EMShape(["scrlabel",.7,.3,0,20.0,20.0,"%d"%self.curtilt,200.0,1])
		self.im2d.add_shape("tilt",s)
		
		if self.tiltshapes[self.curtilt]!=None :
			self.im2d.add_shape("finalbox",self.tiltshapes[self.curtilt])
			
			s0=self.tiltshapes[self.oldtilt].getShape()
			s1=self.tiltshapes[self.curtilt].getShape()
			dx=s0[4]-s1[4]
			dy=s0[5]-s1[5]
			
			self.im2d.set_origin(self.im2d.origin[0]-dx,self.im2d.origin[1]-dy)
			self.oldtilt=self.curtilt
			
		self.im2d.updateGL()

	def change_tilt(self,direc):
		"""When the user presses the up or down arrow"""
		self.oldtilt=self.curtilt
		self.curtilt+=direc
		if self.curtilt<0 : self.curtilt=0
		if self.curtilt>=self.imageparm["nz"] : self.curtilt=self.imageparm["nz"]-1
		
		self.update_tilt()

	def down(self,event,lc):
		"""The event contains the x,y coordinates in window space, lc are the coordinates in image space"""

		if event.buttons()&Qt.LeftButton:
			if event.modifiers()&Qt.ShiftModifier : 
				self.downadjloc=(lc,self.tiltshapes[self.curtilt].getShape()[4:8])
			else :
				self.downloc=lc

	def drag(self,event,lc):
		if self.downloc!=None:
			dx=abs(lc[0]-self.downloc[0])
			dy=abs(lc[1]-self.downloc[1])
			dx=max(dx,dy)	# Make box square
			dx=good_size(dx*2)/2	# use only good sizes
			dy=dx
			s=EMShape(["rectpoint",0,.7,0,self.downloc[0]-dx,self.downloc[1]-dy,self.downloc[0]+dx,self.downloc[1]+dy,1])
			self.im2d.add_shape("box",s)
			s=EMShape(["scrlabel",.7,.7,0,20.0,20.0,"%d (%d x %d)"%(self.curtilt,dx*2,dy*2),200.0,1])
			self.im2d.add_shape("tilt",s)
		elif self.downadjloc!=None:
			dx=(lc[0]-self.downadjloc[0][0])
			dy=(lc[1]-self.downadjloc[0][1])
			s=self.tiltshapes[self.curtilt].getShape()[:]
			s[4]=self.downadjloc[1][0]+dx
			s[5]=self.downadjloc[1][1]+dy
			s[6]=self.downadjloc[1][2]+dx
			s[7]=self.downadjloc[1][3]+dy
			self.im2d.add_shape("box",EMShape(s))

		self.im2d.updateGL()
	
	def up(self,event,lc):
		if self.downloc!=None :
			dx=abs(lc[0]-self.downloc[0])
			dy=abs(lc[1]-self.downloc[1])
			dx=max(dx,dy)	# Make box square
			dx=good_size(dx*2)/2	# use only good sizes
			dy=dx
			s=EMShape(["rectpoint",.7,.2,0,self.downloc[0]-dx,self.downloc[1]-dy,self.downloc[0]+dx,self.downloc[1]+dy,1])
			self.im2d.del_shape("box")
			if hypot(lc[0]-self.downloc[0],lc[1]-self.downloc[1])>5 : 
				self.tiltshapes=[None for i in range(self.imageparm["nz"])]
				self.find_boxes(s)
			
			self.update_tilt()
			self.downloc=None
		elif self.downadjloc!=None :
			dx=(lc[0]-self.downadjloc[0][0])
			dy=(lc[1]-self.downadjloc[0][1])
			s=self.tiltshapes[self.curtilt].getShape()[:]
			s[4]=self.downadjloc[1][0]+dx
			s[5]=self.downadjloc[1][1]+dy
			s[6]=self.downadjloc[1][2]+dx
			s[7]=self.downadjloc[1][3]+dy
			self.tiltshapes[self.curtilt]=EMShape(s)
			self.im2d.add_shape("finalbox",self.tiltshapes[self.curtilt])
			self.im2d.del_shape("box")

			self.update_tilt()
			self.update_stack()
			self.downadjloc=None
			

	def get_boxed_stack(self):
		stack=[]
		for i in range(self.imageparm["nz"]):
			refshape=self.tiltshapes[i].getShape()
			img=EMData(self.imagefile,0,False,Region(refshape[4],refshape[5],i,refshape[6]-refshape[4],refshape[7]-refshape[5],1))
			img["ptcl_source_coord"]=(int((refshape[6]+refshape[4])/2.0),int((refshape[7]+refshape[5])/2.0),i)
			img["ptcl_source_image"]=str(self.imagefile)
			if self.invert : img.mult(-1.0)
			img.process_inplace("normalize.edgemean")
			stack.append(img)
			
		return stack

	def cenalign_stack(self):
		"""This will perform an iterative centering process on a stack of particle images, centering each on the average.
	It will modify the current stack of boxing parameters in-place"""
		
		for it in range(5):
			stack=self.get_boxed_stack()
			
			# Average the stack, and center it
			av=stack[0].copy()
			for im in stack[1:]:
				av.add(im)
			av.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			av.process_inplace("filter.highpass.gauss",{"cutoff_abs":.02})
			av.process_inplace("xform.centeracf")
			#display((av,av2))
			
			# align to the average
			for i,im in enumerate(stack):
				im.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
				im.process_inplace("filter.highpass.gauss",{"cutoff_abs":.02})
				ali=im.align("translational",av)
				trans=ali["xform.align2d"].get_trans()
				shape=self.tiltshapes[i]
				shape.translate(-trans[0],-trans[1])

		# Update the stack display
		stack=self.get_boxed_stack()
		self.imboxed.set_data(stack)
	
	def reconstruct_wedgefill(self,stack,angstep,mode=2):
		"""Fills the missing wedge with the average of the slices"""
		print "Making 3D tomofill"
		
		taxis=float(self.btiltaxisval.text())
		boxsize=stack[0]["nx"]
		pad=Util.calc_best_fft_size(int(boxsize*1.5))

		# average all of the slices together
		av=stack[0].copy()
		for p in stack[1:]: av+=p
		av.del_attr("xform.projection")
		av.mult(1.0/(len(stack)))
		av=av.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
		
		for i,p in enumerate(stack) : 
			p["alt"]=(i-len(stack)/2)*angstep
		
		# Determine a good angular step for filling Fourier space
		fullsamp=360.0/(boxsize*pi)
		if angstep/fullsamp>2.0 :
			samp=1.0/(floor(angstep/fullsamp))
		else :
			samp=angstep
		
		print "Subsampling = %1.2f"%samp

		# Now the reconstruction
		recon=Reconstructors.get("fourier", {"sym":"c1","size":(pad,pad,pad),"mode":reconmodes[mode],"verbose":True})
		recon.setup()
		
		for ri in range(5):
			print "Iteration ",ri
			for a in [i*samp for i in range(-int(90.0/samp),int(90.0/samp)+1)]:
				for ii in range(len(stack)-1):
					if stack[ii]["alt"]<=a and stack[ii+1]["alt"]>a : break
				else: ii=-1
				
				if a<stack[0]["alt"] :
					p=av
					#frac=0.5*(a-stack[0]["alt"])/(-90.0-stack[0]["alt"])
					## a bit wierd. At the ends (missing wedge) we use the average over all tilts. This could be improved
					#p=stack[0].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*(1.0-frac)+stack[-1].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*frac
#					print a," avg ",frac,stack[0]["alt"]
				elif ii==-1 :
					p=av
					#frac=0.5*(a-stack[-1]["alt"])/(90.0-stack[-1]["alt"])+.5
					## a bit wierd. At the ends (missing wedge) we use the average over all tilts. This could be improved
					#p=stack[-1].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*(1.0-frac)+stack[0].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*frac
#					print a," avg ",frac
				else:
					# We average slices in real space, producing a rotational 'smearing' effect
					frac=(a-stack[ii]["alt"])/angstep
					p=stack[ii].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*(1.0-frac)+stack[ii+1].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))*frac
#					print a,ii,ii+1,stack[ii]["alt"],frac
				
				xf=Transform({"type":"eman","alt":a,"az":-taxis,"phi":taxis})
				p["xform.projection"]=xf
						
				if ri%2==1:	
					recon.determine_slice_agreement(p,xf,1)
				else :
					recon.insert_slice(p,xf)
		
		ret=recon.finish()
		print "Done"
		ret=ret.get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
		ret.process_inplace("normalize.edgemean")
#		ret=ret.get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
		
		return ret
	
	def reconstruct_ca(self,stack,angstep,mode=2):
		"""Cylindrically averaged tomographic model, generally used for filling empty spaces. Returned volume is padded."""
		print "Making CA"
		
		taxis=float(self.btiltaxisval.text())
		boxsize=stack[0]["nx"]
		pad=Util.calc_best_fft_size(int(boxsize*1.5))

		# average all of the slices together
		av=stack[0].copy()
		for p in stack[1:]: av+=p
		av.del_attr("xform.projection")
		p.mult(1.0/len(stack))
		av=av.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))

		recon=Reconstructors.get("fourier", {"quiet":True,"sym":"c1","x_in":pad,"y_in":pad})
		recon.setup()
		
		for ri in range(3):
			if ri>0 :
				alt=-180.0
				while (alt<180.0):
					recon.determine_slice_agreement(av,Transform({"type":"eman","alt":alt,"az":-taxis,"phi":taxis}),1)
					alt+=angstep
			alt=-180.0
			while (alt<180.0) :
				recon.insert_slice(av,Transform({"type":"eman","alt":alt,"az":-taxis,"phi":taxis}))
				alt+=angstep
		
		ret=recon.finish()
		ret.process_inplace("normalize.edgemean")
#		ret=ret.get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
		
		return ret
		
	
	def reconstruct(self,stack,angstep,mode=0,initmodel=None):
		""" Tomographic reconstruction of the current stack """
		if initmodel!=None : print "Using initial model"
		
		taxis=float(self.btiltaxisval.text())
		
		boxsize=stack[0]["nx"]
		pad=good_size(int(boxsize*1.5))
		
		for i,p in enumerate(stack) : p["xform.projection"]=Transform({"type":"eman","alt":(i-len(stack)/2)*angstep,"az":-taxis,"phi":taxis})
		
		recon=Reconstructors.get("fourier", {"sym":"c1","size":(pad,pad,pad),"mode":reconmodes[mode],"verbose":True})
		if initmodel!=None : recon.setup(initmodel,.01)
		else : recon.setup()
		scores=[]
		
		# First pass to assess qualities and normalizations
		for i,p in enumerate(stack):
			p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
			p2=recon.preprocess_slice(p2,p["xform.projection"])
			recon.insert_slice(p2,p["xform.projection"],1.0)
			print " %d    \r"%i
		print ""

		# after building the model once we can assess how well everything agrees
		for p in stack:
			p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
			p2=recon.preprocess_slice(p2,p["xform.projection"])
			recon.determine_slice_agreement(p2,p["xform.projection"],1.0,True)
			scores.append((p2["reconstruct_absqual"],p2["reconstruct_norm"]))
			print " %d\t%1.3f    \r"%(i,scores[-1][0])
		print ""

		# clear out the first reconstruction (probably no longer necessary)
#		ret=recon.finish(True)
#		ret=None
		
		# setup for the second run
		if initmodel!=None : recon.setup(initmodel,.01)
		else : recon.setup()

		thr=0.7*(scores[len(scores)/2][0]+scores[len(scores)/2-1][0]+scores[len(scores)/2+1][0])/3;		# this is rather arbitrary
		# First pass to assess qualities and normalizations
		for i,p in enumerate(stack):
			if scores[i][0]<thr : 
				print "%d. %1.3f *"%(i,scores[i][0])
				continue
			
			print "%d. %1.2f \t%1.3f\t%1.3f"%(i,p["xform.projection"].get_rotation("eman")["alt"],scores[i][0],scores[i][1])
			p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
			p2=recon.preprocess_slice(p2,p["xform.projection"])
			p2.mult(scores[i][1])
			recon.insert_slice(p2,p["xform.projection"],1.0)

#		plot(scores)
		
		recon.set_param("savenorm","norm.mrc")
		ret=recon.finish(True)
		ret=ret.get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
#		print "Quality: ",qual
		
		return ret
	
	def tilt_axis(self):
		ntilt=self.imageparm["nz"]
		sz=good_size(self.imageparm["nx"]/2)
		while 1:
			av=None
			n=0
			for i in range(ntilt):
				refshape=self.tiltshapes[i].getShape()
				if refshape[4]<=sz/2 or refshape[5]<=sz/2 or self.imageparm["nx"]-refshape[4]<=sz/2 or self.imageparm["ny"]-refshape[5]<=sz/2 : break
				img=EMData(self.imagefile,0,False,Region(refshape[4]-sz/2,refshape[5]-sz/2,i,sz,sz,1))
				if self.invert : img.mult(-1.0)
				img.process_inplace("normalize.edgemean")
				
				if av==None: av=img
				else : av.add(img)
				n+=1
				
			if n==ntilt : break
			sz/=2
			if sz<32: return
			print "You may wish to center on a feature closer to the center of the image next time -> ",sz
		
		sz2=good_size(sz+128)
		av2=av.get_clip(Region((sz-sz2)/2,(sz-sz2)/2,sz2,sz2))
		av2.process_inplace("mask.zeroedgefill")
		av2.process_inplace("filter.flattenbackground",{"radius":64})
		av=av2.get_clip(Region((sz2-sz)/2,(sz2-sz)/2,sz,sz))
		av.process_inplace("normalize.edgemean")
		av.process_inplace("mask.sharp",{"outer_radius":sz/2-1})

#		display(av)
		f=av.do_fft()
		d=f.calc_az_dist(360,-90.25,0.5,10.0,sz/2-1)
		d=[(i,j*0.5-90) for j,i in enumerate(d)]
		self.btiltaxisval.setText(str(max(d)[1]))
#		print max(d)
#		print min(d)
#		plot(d)

	
	def projection_align(self,angstep=2.0):
		"""realign the current set of boxes using iterative projection matching"""
		
		taxis=float(self.btiltaxisval.text())

		stack=self.get_boxed_stack()
		for i,p in enumerate(stack) : 
			ort=Transform({"type":"eman","alt":(i-len(stack)/2)*angstep,"az":-taxis,"phi":taxis})		# is this right ?
			curshape=self.tiltshapes[i].getShape()
			
			# Read the reference at the user specified size, then pad it a bit
			ref=self.map3d.project("standard",ort)
			ref.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			ref.process_inplace("normalize.edgemean")
			ref=ref.get_clip(Region(-self.maxshift,-self.maxshift,ref["nx"]+self.maxshift*2,ref["ny"]+self.maxshift*2))
			
			# when we read the alignment target, we pad with actual image data since the object will have moved
			trg=EMData(self.imagefile,0,False,Region(curshape[4]-self.maxshift,curshape[5]-self.maxshift,i,curshape[6]-curshape[4]+self.maxshift*2,curshape[7]-curshape[5]+self.maxshift*2,1))
			trg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			trg.process_inplace("normalize.edgemean")
			if self.invert : trg.mult(-1.0)
			
			aln=ref.align("translational",trg,{"intonly":1,"maxshift":self.maxshift*4/5})
			trans=aln["xform.align2d"].get_trans()
			print i,trans[0],trans[1]
			if i>len(stack)-4 : display([ref,trg,aln])
#			if i==self.curtilt+3 : display((ref,trg,aln,ref.calc_ccf(trg)))

			self.tiltshapes[i].translate(trans[0],trans[1])

		
	
	def find_boxes(self,mainshape):
		"""Starting with a user selected box at the current tilt, search for the same shape in the entire
	tilt series"""
	
		if self.imagefile==None: return 
		
		self.tiltshapes[self.curtilt]=mainshape
		
		lref=None
		for i in range(self.curtilt+1,self.imageparm["nz"]):
			refshape=self.tiltshapes[i-1].getShape()
			
			# Read the reference at the user specified size, then pad it a bit
			ref=EMData(self.imagefile,0,False,Region(refshape[4],refshape[5],i-1,refshape[6]-refshape[4],refshape[7]-refshape[5],1))
			ref.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0})
			ref.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			ref.process_inplace("normalize.edgemean")
			ref=ref.get_clip(Region(-self.maxshift,-self.maxshift,ref["nx"]+self.maxshift*2,ref["ny"]+self.maxshift*2))
			if lref!=None and self.seqali : ref.add(lref)
			ref.process_inplace("normalize.edgemean")		# older images contribute less
			lref=ref
			
			# when we read the alignment target, we pad with actual image data since the object will have moved
			trg=EMData(self.imagefile,0,False,Region(refshape[4]-self.maxshift,refshape[5]-self.maxshift,i,refshape[6]-refshape[4]+self.maxshift*2,refshape[7]-refshape[5]+self.maxshift*2,1))
			trg.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0})
			trg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			trg.process_inplace("normalize.edgemean")
			
			aln=ref.align("translational",trg,{"intonly":1,"maxshift":self.maxshift*4/5,"masked":1})
			ref.write_image("dbug.hdf",-1)
			trg.write_image("dbug.hdf",-1)
			aln.write_image("dbug.hdf",-1)
			trans=aln["xform.align2d"].get_trans()
#			if i==self.curtilt+3 : display((ref,trg,aln,ref.calc_ccf(trg)))

			self.tiltshapes[i]=EMShape(["rectpoint",.7,.2,0,refshape[4]+trans[0],refshape[5]+trans[1],refshape[6]+trans[0],refshape[7]+trans[1],1])
			print i,trans[0],trans[1]
			
		lref=None
		for i in range(self.curtilt-1,-1,-1):
			refshape=self.tiltshapes[i+1].getShape()
			
			# Read the reference at the user specified size, then pad it a bit
			ref=EMData(self.imagefile,0,False,Region(refshape[4],refshape[5],i+1,refshape[6]-refshape[4],refshape[7]-refshape[5],1))
			ref.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			ref.process_inplace("normalize.edgemean")
			ref=ref.get_clip(Region(-self.maxshift,-self.maxshift,ref["nx"]+self.maxshift*2,ref["ny"]+self.maxshift*2))
			if lref!=None and self.seqali : ref.add(lref)
			ref.process_inplace("normalize.edgemean")
			lref=ref
			
			# when we read the alignment target, we pad with actual image data since the object will have moved
			trg=EMData(self.imagefile,0,False,Region(refshape[4]-self.maxshift,refshape[5]-self.maxshift,i,refshape[6]-refshape[4]+self.maxshift*2,refshape[7]-refshape[5]+self.maxshift*2,1))
			trg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			trg.process_inplace("normalize.edgemean")
			
			aln=ref.align("translational",trg,{"intonly":1,"maxshift":self.maxshift*4/5,"masked":1})
			trans=aln["xform.align2d"].get_trans()
			if i==self.curtilt+3 : display((ref,trg,aln,ref.calc_ccf(trg)))

			self.tiltshapes[i]=EMShape(["rectpoint",.7,.2,0,refshape[4]+trans[0],refshape[5]+trans[1],refshape[6]+trans[0],refshape[7]+trans[1],1])
			print i,trans[0],trans[1]
		
		self.update_stack()


if __name__ == "__main__":
	main()

		
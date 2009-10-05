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
from emapplication import get_application, EMStandAloneApplication,EMQtWidgetModule
from emimage2d import EMImage2DModule
import weakref
from emshape import EMShape
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import math
#import EMAN2db
import time

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input ali>
	
This program will allow the user to select a specific feature within a tomogram and track it, boxing out the
feature from all slices. Generally best for uniform objects like vesicles."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--prefix",type="string",help="Output file prefix",default=None)
	parser.add_option("--maxshift",type="int",help="Maximum shift in autoalign, default=8",default=8)
	parser.add_option("--seqali",action="store_true",default=False,help="Average particles in previous slices for better alignments of near-spherical objects")

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
	
	app = EMStandAloneApplication()

	
	track=TrackerControl(app,options.maxshift,options.seqali)
	track.set_image(args[0])

	app.execute()
	
class TrackerControl():
	def __init__(self,app,maxshift,seqali):
		self.maxshift=maxshift
		self.seqali=seqali
		
		# the single image display widget
		self.im2d = EMImage2DModule(application=app)
	
		# get some signals from the window. With a regular Qt object this would just
		# be the window, but with 'Module' objects we must use the emitter
		QtCore.QObject.connect(self.im2d.emitter(),QtCore.SIGNAL("mousedown"),self.down)
		QtCore.QObject.connect(self.im2d.emitter(),QtCore.SIGNAL("mousedrag"),self.drag)
		QtCore.QObject.connect(self.im2d.emitter(),QtCore.SIGNAL("mouseup"),self.up)
		QtCore.QObject.connect(self.im2d.emitter(),QtCore.SIGNAL("increment_list_data"),self.change_tilt)
	
		self.imagefile=None
		self.imageparm=None
		self.tiltshapes=None
		self.curtilt=0
		
		self.im2d.show()

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
		self.im2d.set_data(self.curimg)

		s=EMShape(["scrlabel",.7,.7,0,20.0,20.0,"%d"%self.curtilt,200.0,1])
		self.im2d.add_shape("tilt",s)
		
		if self.tiltshapes[self.curtilt]!=None :
			self.im2d.add_shape("box",self.tiltshapes[self.curtilt])
		self.im2d.updateGL()

	def change_tilt(self,direc):
		"""When the user presses the up or down arrow"""
		self.curtilt+=direc
		if self.curtilt<0 : self.curtilt=0
		if self.curtilt>=self.imageparm["nz"] : self.curtilt=self.imageparm["nz"]-1
		
		self.update_tilt()

	def down(self,event,lc):
		"""The event contains the x,y coordinates in window space, lc are the coordinates in image space"""

		self.downloc=lc	
		self.tiltshapes=[None for i in range(self.imageparm["nz"])]

	def drag(self,event,lc):
		dx=abs(lc[0]-self.downloc[0])
		dy=abs(lc[1]-self.downloc[1])
		s=EMShape(["rectpoint",0,.7,0,self.downloc[0]-dx,self.downloc[1]-dy,self.downloc[0]+dx,self.downloc[1]+dy,1])
		self.im2d.add_shape("box",s)
		s=EMShape(["scrlabel",.7,.7,0,20.0,20.0,"%d (%d x %d)"%(self.curtilt,dx*2,dy*2),200.0,1])
		self.im2d.add_shape("tilt",s)

		self.im2d.updateGL()
	
	def up(self,event,lc):
		dx=abs(lc[0]-self.downloc[0])
		dy=abs(lc[1]-self.downloc[1])
		s=EMShape(["rectpoint",.7,.2,0,self.downloc[0]-dx,self.downloc[1]-dy,self.downloc[0]+dx,self.downloc[1]+dy,1])
		self.im2d.del_shape("box")
		if hypot(lc[0]-self.downloc[0],lc[1]-self.downloc[1])>5 : self.find_boxes(s)
		
		self.update_tilt()

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
			ref.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			ref.process_inplace("normalize.edgemean")
			ref=ref.get_clip(Region(-self.maxshift,-self.maxshift,ref["nx"]+self.maxshift*2,ref["ny"]+self.maxshift*2))
			if lref!=None and self.seqali : ref.add(lref)
			lref=ref
			
			# when we read the alignment target, we pad with actual image data since the object will have moved
			trg=EMData(self.imagefile,0,False,Region(refshape[4]-self.maxshift,refshape[5]-self.maxshift,i,refshape[6]-refshape[4]+self.maxshift*2,refshape[7]-refshape[5]+self.maxshift*2,1))
			trg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			trg.process_inplace("normalize.edgemean")
			
			aln=ref.align("translational",trg,{"intonly":1,"maxshift":self.maxshift*4/5})
			trans=aln["xform.align2d"].get_trans()
			if i==self.curtilt+3 : display((ref,trg,aln,ref.calc_ccf(trg)))

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
			lref=ref
			
			# when we read the alignment target, we pad with actual image data since the object will have moved
			trg=EMData(self.imagefile,0,False,Region(refshape[4]-self.maxshift,refshape[5]-self.maxshift,i,refshape[6]-refshape[4]+self.maxshift*2,refshape[7]-refshape[5]+self.maxshift*2,1))
			trg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			trg.process_inplace("normalize.edgemean")
			
			aln=ref.align("translational",trg,{"intonly":1,"maxshift":self.maxshift*4/5})
			trans=aln["xform.align2d"].get_trans()
			if i==self.curtilt+3 : display((ref,trg,aln,ref.calc_ccf(trg)))

			self.tiltshapes[i]=EMShape(["rectpoint",.7,.2,0,refshape[4]+trans[0],refshape[5]+trans[1],refshape[6]+trans[0],refshape[7]+trans[1],1])
			print i,trans[0],trans[1]
			
if __name__ == "__main__":
	main()

		
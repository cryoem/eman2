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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# e2aligntest.py  09/21/2004  Steven Ludtke
# This program is used to generate various alignment test images


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
#import pdb
#from bisect import insort
import numpy as np
from numpy import arange

BOXSIZE=1024

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <movie stack>"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--axes",type="string",help="String list 3 axes from xyzaqp",default="xaq")
	
	(options, args) = parser.parse_args()
#	if len(args)<2 : parser.error("Input1 and Input2 required")

	n=EMUtil.get_image_count(args[0])
	h=EMData(args[0],0,True)
	nx=h["nx"]
	ny=h["ny"]
	
	N=0
	clpav=[]
	for im in xrange(n):
		print im,"/",n
		img=EMData(args[0],im)
		rgn=0
		for x in xrange(50,nx-BOXSIZE-50,BOXSIZE):
			for y in xrange(50,ny-BOXSIZE-50,BOXSIZE):
				clp=img.get_clip(Region(x,y,BOXSIZE,BOXSIZE))
				clp.process_inplace("normalize.unitlen");
				clpf=clp.do_fft()
				
				# coherent average
				try: clpav[rgn].add(clp)
				except: 
					clpav.append(clp.copy())
				
				# incoherent average
				clpf.ri2inten()
				try: pws.add(clpf)
				except: 
					pws=clpf
				N+=1
				rgn+=1
	
#	display(clpav,True)
	for im in xrange(len(clpav)):
		clpav[im]=clpav[im].do_fft()
		clpav[im].ri2inten()
		
		try: pwsc.add(clpav[im])
		except: 
			pwsc=clpav[im]
		
		clpav[im]=None
					
	pws.mult(1.0/N)
	pws.process_inplace("math.sqrt")
	pws["is_intensity"]=0			# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
	pws["is_complex_ri"]=0
	print pws.get_attr_dict()
#	pws.sub(pws["mean"]-pws["sigma"]/2.0)
	pws["is_complex_ri"]=1
	pws[0,0]=0
	pws.do_ift().write_image("pws.hdf",0)

	pwsc.mult(1.0/N)
	pwsc[0,0]=0
	pwsc.process_inplace("math.sqrt")
	pwsc["is_intensity"]=0			# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
	pwsc.do_ift().write_image("pws.hdf",1)
	
#	ds=1.0/(pws["apix_x"]*ny)
#	pws1=pws.calc_radial_dist(ny/2,0,1.0,1)
#	pws1bg=ctf.low_bg_curve(pws1,ds)

	df,flipim=fit_defocus(pws)
	
	img1=EMData.read_images(args[0],(5,6,7,8))
	for i in img1:
		try: img1c.add(i.get_clip(Region(2000,2000,BOXSIZE,BOXSIZE)))
		except: img1c=i.get_clip(Region(2000,2000,BOXSIZE,BOXSIZE))
		
	img2=EMData.read_images(args[0],(12,13,14,15))
	for i in img2:
		try: img2c.add(i.get_clip(Region(2000,2000,BOXSIZE,BOXSIZE)))
		except: img2c=i.get_clip(Region(2000,2000,BOXSIZE,BOXSIZE))

	ccf1=img1c.calc_ccf(img2c)
	ccf1.process_inplace("normalize")
	ccf1.process_inplace("xform.phaseorigin.tocenter")
	img2cf=img2c.do_fft()
	img2cf.mult(flipim)
	img2cp=img2cf.do_ift()
	ccf2=img1c.calc_ccf(img2cp)
	ccf2.process_inplace("normalize")
	ccf1.process_inplace("xform.phaseorigin.tocenter")
	
	display((ccf1,ccf2),True)
	
def fit_defocus(img):
	ds=1.0/(img["apix_x"]*img["ny"])
	ns=min(int(floor(.25/ds)),img["ny"]/2)

	# the data curve we are trying to fit
	oned=np.array(img.calc_radial_dist(ns,0,1.0,1)[1:])
	oned-=min(oned)	# get rid of bulk background
	oned/=max(oned)	# normalize a bit for convienience

	# so we can generate simulated curves
	ctf=EMAN2Ctf()
	ctf.voltage=300.0
	ctf.cs=4.7
	ctf.ampcont=10
	ctf.bfactor=0
	ctf.dfdiff=0
	s=[ds*(i+1) for i,j in enumerate(oned)]
	
	dfl=[]
	ql=[]
	for df in arange(0.6,5,.01):
		ctf.defocus=df
		curve=np.array(ctf.compute_1d(ns*2,ds,Ctf.CtfType.CTF_AMP)[1:])
		# we square the curve (no B-factor), then "normalize" it so constant background won't enter the fit other than edge effects
		curve*=curve
		curve-=curve.mean()
#		plot((s,list(oned)),(s,list(curve)))
		
		qual=curve.dot(oned)

		print df,qual
		dfl.append(df)
		ql.append(qual)
	
	qls=np.convolve(ql,[.2,.2,.2,.2,.2,.2,.2],mode="same")
	qls[0]=qls[1]=qls[2]=qls[-3]=qls[-2]=qls[-1]=-2
	a=np.argmax(qls)
	df=dfl[a]
	
	print "Best defocus ",df
	plot((dfl,ql),(dfl,qls))
	
	ctf.defocus=df
	aliimg=img.copy()
	ctf.compute_2d_complex(aliimg,Ctf.CtfType.CTF_ALIFILT,None)
	
	display(aliimg)
	return df,aliimg

if __name__ == "__main__":  main()

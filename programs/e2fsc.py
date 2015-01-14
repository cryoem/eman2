#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/20/2012 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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



from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
import time
from numpy import *


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2fsc.py [options] input1 input2

Simple 2 volume FSCs can be computed with e2proc3d.py. In addition to the overall fsc (saved to fsc.txt), 
it also computes a "local resolution" through the volume. This method is of HIGHLY questionable usefulness,
and this program should be regarded as experimental.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--input",type=str,help="Similarity matrix to analyze",default=None)
#	parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
	parser.add_argument("--output",type=str,help="Output text file",default="zvssim.txt")
	parser.add_argument("--localsize", type=int, help="Size in pixels of the local region to compute the resolution in",default=-1)
	parser.add_argument("--overlap", type=int, help="Amount of oversampling to use in local resolution windows. Larger value -> larger output map",default=12)
	parser.add_argument("--apix", type=float, help="A/pix to use for the comparison (default uses Vol1 apix)",default=0)
	#parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
	#parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
	#parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
	#parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbosity [0-9]. Larger values produce more output.")

	print "WARNING: This program is considered highly experimental, and there are mathematical \narguments that local estimation techniques will not produce reliable values.\n"
	print "Having said that, the fsc.txt file is a normal FSC between the two volumes, and IS \nreliable, though e2proc3d.py could compute it far more easily"

	(options, args) = parser.parse_args()
		
	if len(args)<2 : 
		print "Please specify 2 input files"
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	v1=EMData(args[0],0)
	v2=EMData(args[1],0)
	
	if options.apix>0 : apix=options.apix
	else :
		apix=v1["apix_x"]
		print "Using %1.2f A/pix"%apix
	
	nx,ny,nz=v1["nx"],v1["ny"],v1["nz"]
	print "%d x %d x %d"%(nx,ny,nz)
	if nx!=ny or nx!=nz : print "Warning: non-cubic volumes may produce unexpected results"
	
	if options.localsize==-1 : lnx=nx/10*2
	else: lnx=options.localsize
	if apix*lnx/2.0<15.0 :
		print "WARNING: Local sampling box is <15 A. Adjusting to 20 A."
		lnx=int(floor(40.0/apix))
	print "Local region is %d pixels"%lnx
	
	thresh1=v1["mean"]+v1["sigma"]
	thresh2=v2["mean"]+v2["sigma"]
	print "Thresholds : ",thresh1,thresh2
	
	# overall fsc
	fsc=v1.calc_fourier_shell_correlation(v2)
	fx=array(fsc[0:len(fsc)/3])/apix
	fy=fsc[len(fsc)/3:len(fsc)*2/3]

	out=file("fsc.txt","w")
	for i,x in enumerate(fx):
		out.write("%1.5f\t%1.4f\n"%(x,fy[i]))
	out.close()
	
	# Create a centered Gaussian mask with a size ~1/10th of the box size
	cenmask=EMData(lnx,lnx,lnx)
	cenmask.to_one()
	cenmask.process_inplace("mask.gaussian",{"inner_radius":lnx/6,"outer_radius":lnx/6})
	print "Approx feature size for assessment = %1.1f A"%(apix*lnx/2.0)
	cenmask.write_image("cenmask.hdf")
	#display(cenmask)
	
	overlap=options.overlap		# This is the fraction of the window size to use as a step size in sampling
	if overlap<1 or overlap>lnx :
		print "Invalid overlap specified, using default"
	
	xr=xrange(0,nx-lnx,lnx/overlap)
	yr=xrange(0,ny-lnx,lnx/overlap)
	zr=xrange(0,nz-lnx,lnx/overlap)
	resvol=EMData(len(xr),len(yr),len(zr))
	resvol["apix_x"]=apix*lnx/overlap
	resvol["apix_y"]=apix*lnx/overlap
	resvol["apix_z"]=apix*lnx/overlap
	resvol143=resvol.copy()
	
	fys=[]
	funny=[]		# list of funny curves
	t=time.time()
	for oz,z in enumerate(zr):
		for oy,y in enumerate(yr):
			for ox,x in enumerate(xr):
				if options.verbose>1 : print "%d, %d, %d :"%(x,y,z),
				elif options.verbose:
					if time.time()-t>.5 :
						print "  %3d,%3d,%3d / %d,%d,%d\r"%(ox,oy,oz,len(xr),len(yr),len(zr)),
						sys.stdout.flush()
						t=time.time()
						
				v1m=v1.get_clip(Region(x,y,z,lnx,lnx,lnx))
				v1m.mult(cenmask)
				
				v2m=v2.get_clip(Region(x,y,z,lnx,lnx,lnx))
				v2m.mult(cenmask)
				
				
				## make a copy of the mask centered at the current position
				#mask=cenmask.process("xform.translate.int",{"trans":(x-nx/2,y-ny/2,z-nz/2)})		

				#v1m=v1*mask
				#v2m=v2*mask

				if options.verbose>1 : print v1m["sigma_nonzero"], v2m["sigma_nonzero"],
				if v1m["maximum"]<thresh1 or v2m["maximum"]<thresh2 :
					if options.verbose>1 : print " "
					resvol[ox,oy,oz]=0.0
					continue
				if options.verbose>1 : print " ***"
#				display(v1m)
				
				fsc=v1m.calc_fourier_shell_correlation(v2m)
				fx=array(fsc[1:len(fsc)/3])/apix
				fy=fsc[len(fsc)/3+1:len(fsc)*2/3]
				
				# 0.5 resolution
				for i,xx in enumerate(fx[:-1]):
					if fy[i]>0.5 and fy[i+1]<0.5 : break
				res=(0.5-fy[i])*(fx[i+1]-fx[i])/(fy[i+1]-fy[i])+fx[i]
				if res<0 : res=0.0
				if res>fx[-1]: 
					res=fx[-1]		# This makes the resolution at Nyquist, which is not a good thing
					funny.append(len(fys))
#				if res>0 and res<0.04 : funny.append(len(fys))
				resvol[ox,oy,oz]=res
				
				fys.append(fy)
				if isnan(fy[0]): print "NAN"

				# 0.143 resolution
				for i,xx in enumerate(fx[:-1]):
					if fy[i]>0.143 and fy[i+1]<0.143 : break
				res143=(0.143-fy[i])*(fx[i+1]-fx[i])/(fy[i+1]-fy[i])+fx[i]
				if res143<0 : res143=0.0
				if res143>fx[-1]: 
					res143=fx[-1]		# This makes the resolution at Nyquist, which is not a good thing
#				if res>0 and res<0.04 : funny.append(len(fys))
				resvol143[ox,oy,oz]=res143
				
				fys.append(fy)


	resvol.write_image("resvol.hdf")
	resvol143.write_image("resvol143.hdf")
	
	out=file("fsc.curves.txt","w")
	out.write("# This file contains individual FSC curves from e2fsc.py. Only a fraction of computed curves are included.\n")
	if len(fys)>100 : 
		step=len(fys)/100
		print "Saving 1/%d of curves to fsc.curves.txt + %d"%(step,len(funny))
	else: 
		step=1
		print "Saving all curves to fsc.curves.txt"
	
	for i,x in enumerate(fx):
		out.write( "%f\t"%x)
		for j in range(0,len(fys),step):
			out.write( "%f\t"%fys[j][i])
		
		# Also save any particularly low resolution curves
		for j in funny:
			out.write( "%f\t"%fys[j][i])
			
		out.write("\n")
		
	if len(funny)>1 :
		print "WARNING: %d/%d curves were evaluated as being >0.5 AT Nyquist. While these values have been set to \
		Nyquist (the maximum resolution for your sampling), these values are not meaningful, and could imply \
		insufficient sampling, or bias in the underlying reconstructions."%(len(funny),len(fys))

	E2end(logid)

if __name__ == "__main__":
        main()


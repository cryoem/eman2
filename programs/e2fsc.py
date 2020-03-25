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



from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
from math import *
import os
import sys
import time
from numpy import *
import queue

def procthread(jsd,vals,lnx,thresh1,thresh2,apix,v1,v2,cenmask,avgmask,options,ttl):
	for ox,x,oy,y,oz,z in vals:
		if options.verbose>2 : print("%d, %d, %d :"%(x,y,z), end=' ')
		elif options.verbose>3:
			if time.time()-t>.5 :
				print("  %3d,%3d,%3d / %d,%d,%d\r"%(ox,oy,oz,len(xr),len(yr),len(zr)), end=' ')
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

		if options.verbose>3 : print(v1m["sigma_nonzero"], v2m["sigma_nonzero"], end=' ')
		if v1m["maximum"]<thresh1 or v2m["maximum"]<thresh2 :
			if options.verbose>3 : print(" ")
			jsd.put((0,0,[],ox,x,oy,y,oz,z,0,None,None))
			continue
		
		if options.verbose>3 : print(" ***")
	#				display(v1m)
		
		fsc=v1m.calc_fourier_shell_correlation(v2m)
		fx=old_div(array(fsc[1:old_div(len(fsc),3)]),apix)
		fy=fsc[old_div(len(fsc),3)+1:old_div(len(fsc)*2,3)]
		
		# 0.5 resolution
		if fy[0]<0.5 and fy[1]<0.5 : i,xx,res=1,fx[1],fx[1]
		else:
			for i,xx in enumerate(fx[:-1]):
				if fy[i]>0.5 and fy[i+1]<0.5 : break
			res=(0.5-fy[i])*(fx[i+1]-fx[i])/(fy[i+1]-fy[i])+fx[i]
		if res<0 : res=0.0
		if res>fx[-1]: 
			res=fx[-1]		# This makes the resolution at Nyquist, which is not a good thing
	#		funny.append(len(fys))
	#				if res>0 and res<0.04 : funny.append(len(fys))
		
		if isnan(fy[0]): print("NAN",x,y,z)

		cutoff=options.cutoff
		# 0.143 resolution
		if fy[0]<cutoff and fy[1]<cutoff : si,xx,res143=1,fx[1],fx[1]
		else:
			for si,xx in enumerate(fx[:-1]):
				if fy[si]>cutoff and fy[si+1]<cutoff : break
			res143=(cutoff-fy[si])*(fx[si+1]-fx[si])/(fy[si+1]-fy[si])+fx[si]
		if res143<0 : res143=0.0
		if res143>fx[-1]: 
			res143=fx[-1]		# This makes the resolution at Nyquist, which is not a good thing
	#				if res>0 and res<0.04 : funny.append(len(fys))

		# now we build the locally filtered volume
		v1m=v1.get_clip(Region(x,y,z,lnx,lnx,lnx))
		v2m=v2.get_clip(Region(x,y,z,lnx,lnx,lnx))
	#				if res143>.23 : v1m.write_image("zones.hdf",-1)
	
		
		if options.gauss:
			filtername="filter.lowpass.gauss"
		else:
			filtername="filter.lowpass.tophat"
		v1m.process_inplace(filtername,{"cutoff_pixels":si+1})	# sharp low-pass at 0.143 cutoff
		v2m.process_inplace(filtername,{"cutoff_pixels":si+1})	# sharp low-pass at 0.143 cutoff
	#				if res143>.23 : v1m.write_image("zones.hdf",-1)
		v1m.mult(avgmask)
		v2m.mult(avgmask)

	#				if res143>.2 : print x,y,z,si,lnx,fx[si],res143

		jsd.put((res,res143,fy,ox,x,oy,y,oz,z,si,v1m,v2m))
	jsd.put(ttl)		# joins the thread

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2fsc.py [options] input1 input2

Simple 2 volume FSCs can be computed with e2proc3d.py. In addition to the overall fsc (saved to fsc.txt), 
it also computes a "local resolution" through the volume. These local resolutions are based on poor statistics.
The smaller the box size, the worse they are, and this is not taken into account when computing actual
resolution values. Regardless, it will give a reasonable comparison of how much variation exists in different
regions of the map, and will produce locally filtered maps with a reasonable level of detail, given the two
input volumes.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--input",type=str,help="Similarity matrix to analyze",default=None)
#	parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
	parser.add_argument("--output",type=str,help="Output .143 resolution volume",default="resvol143.hdf")
	parser.add_argument("--outfilt",type=str,help="Output locally filtered average volume",default="res143_filtered.hdf")
	parser.add_argument("--outfilte",type=str,help="Apply the local filter to the even map as well and write to specified file",default=None)
	parser.add_argument("--outfilto",type=str,help="Apply the local filter to the odd map as well and write to specified file",default=None)
	parser.add_argument("--localsize", type=int, help="Size in pixels of the local region to compute the resolution in",default=-1)
	parser.add_argument("--overlap", type=int, help="Amount of oversampling to use in local resolution windows. Larger value -> larger output map",default=4)
	parser.add_argument("--apix", type=float, help="A/pix to use for the comparison (default uses Vol1 apix)",default=0)
	parser.add_argument("--cutoff", type=float, help="fsc cutoff. default is 0.143",default=0.143)
	parser.add_argument("--gauss", action="store_true", help="use gaussian filter instead of tophat",default=False)
	parser.add_argument("--mask",type=str,help="Mask to apply to both input images before calculation",default=None)
	#parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
	#parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
	#parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
	#parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on the local computer")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbosity [0-9]. Larger values produce more output.")

	#print "WARNING: This program is considered highly experimental, and there are mathematical \narguments that local estimation techniques will not produce reliable values.\n"
	#print "Having said that, the fsc.txt file is a normal FSC between the two volumes, and IS \nreliable, though e2proc3d.py could compute it far more easily"

	(options, args) = parser.parse_args()
		
	if len(args)<2 : 
		print("Please specify 2 input files")
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	if options.threads<2 : options.threads=2

	v1=EMData(args[0],0)
	v2=EMData(args[1],0)
	if options.mask!=None:
		mask=EMData(options.mask)
		v1.mult(mask)
		v2.mult(mask)
	
	if options.apix>0 : apix=options.apix
	else :
		apix=v1["apix_x"]
		print("Using %1.2f A/pix"%apix)
	
	nx,ny,nz=v1["nx"],v1["ny"],v1["nz"]
	print("%d x %d x %d"%(nx,ny,nz))
	if nx!=ny or nx!=nz : print("Warning: non-cubic volumes may produce unexpected results")
	
	overlap=options.overlap		# This is the fraction of the window size to use as a step size in sampling
	if overlap<1:
		print("Invalid overlap specified, using default")
		overlap=6
		
	if options.localsize==-1 : 
		lnx=int(old_div(32,apix))
		if lnx<16: lnx=16
		lnx=(((lnx-1)//overlap)+1)*overlap
	else: lnx=options.localsize
	if apix*lnx/2.0<10.0 :
		print("WARNING: Local sampling box is <10 A. Adjusting to 16 A.")
		lnx=int(floor(32.0/apix))
	print("Local region is %d pixels"%lnx)
	if overlap>lnx : overlap=lnx
	
	thresh1=v1["mean"]+v1["sigma"]
	thresh2=v2["mean"]+v2["sigma"]
	print("Thresholds : ",thresh1,thresh2)
	
	if options.verbose: print("Computing overall FSC")
	# overall fsc
	fsc=v1.calc_fourier_shell_correlation(v2)
	fx=old_div(array(fsc[0:old_div(len(fsc),3)]),apix)
	fy=fsc[old_div(len(fsc),3):old_div(len(fsc)*2,3)]

	out=open("fsc.txt","w")
	for i,x in enumerate(fx):
		out.write("%1.5f\t%1.4f\n"%(x,fy[i]))
	out.close()
	
	if options.verbose: print("Preparing for local calculation")
	# Create a centered Gaussian mask with a size ~1/10th of the box size
	cenmask=EMData(lnx,lnx,lnx)
	cenmask.to_one()
	cenmask.process_inplace("mask.gaussian",{"inner_radius":old_div(lnx,6),"outer_radius":old_div(lnx,6)})
	print("Approx feature size for assessment = %1.1f A"%(apix*lnx/2.0))
#	cenmask.write_image("cenmask.hdf")
	#display(cenmask)
	
	# Create a Gaussian with the correct size to produce a flat average in 3-D
	avgmask=EMData(lnx,lnx,lnx)
	avgmask.to_one()
	d=float(lnx//overlap)
#	avgmask.process_inplace("mask.gaussian",{"outer_radius":2.0*d/log(8.0) })	# this mask is adjusted to the precise width necessary so a sum of tiled overlapping Gaussians will be flat
	avgmask.process_inplace("mask.gaussian",{"outer_radius":3.0*d/log(8.0) })	# make it a bit wider since we are weighting anyway, this should produce smoother surfaces
	
	off=old_div((nx%(lnx//overlap)),2)
	xr=list(range(off,nx-lnx,lnx//overlap))
	yr=list(range(off,ny-lnx,lnx//overlap))
	zr=list(range(off,nz-lnx,lnx//overlap))
	resvol=EMData(len(xr),len(yr),len(zr))
	resvol["apix_x"]=apix*lnx//overlap
	resvol["apix_y"]=apix*lnx//overlap
	resvol["apix_z"]=apix*lnx//overlap
	resvol143=resvol.copy()
	
	print("Local region: ",lnx," with step ",lnx//overlap)
	
	# volfilt will contain the locally filtered version of the map
	volfilt=v1.copy()
	volfilt.to_zero()
	volnorm=v1.copy()
	volnorm.to_zero()

	if options.outfilte!=None : 
		volfilte=v1.copy()
		volfilte.to_zero()
	
	if options.outfilto!=None : 
		volfilto=v1.copy()
		volfilto.to_zero()
	
	# now do all of the tiled calculations
	fys=[]
	funny=[]		# list of funny curves
	t=time.time()
	jsd=queue.Queue(0)
	thrds=[]
	for oz,z in enumerate(zr):
		for oy,y in enumerate(yr):
			thrd=[]
			for ox,x in enumerate(xr):
				thrd.append((ox,x,oy,y,oz,z))		# this is different from the normal parallelism pattern. We make the actual Threads on demand
			thrds.append(thrd)						# one row of x values per thread
	

	if options.verbose: print(len(thrds)," threads")
	
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose>1 : 
				print("\rStarting thread {}/{}                                                  ".format(thrtolaunch,len(thrds)), end=' ')
				sys.stdout.flush()
			vals=thrds[thrtolaunch]
			thr=threading.Thread(target=procthread,args=(jsd,vals,lnx,thresh1,thresh2,apix,v1,v2,cenmask,avgmask,options,thrtolaunch))
			thrds[thrtolaunch]=thr		# replace the values with the actual thread
			thr.start()
			thrtolaunch+=1
		else: time.sleep(1)
#		if options.verbose>1: print "{}% complete".format(100.0*frac),
	
		while not jsd.empty():
			vals=jsd.get()
			if isinstance(vals,int) : 
				thrds[vals].join()
				thrds[vals]=None
			else :
				res,res143,fy,ox,x,oy,y,oz,z,si,v1m,v2m=vals

				resvol[ox,oy,oz]=res
				resvol143[ox,oy,oz]=res143
				fys.append(fy)
				if res==0 : continue


				volfilt.insert_scaled_sum(v1m,(x+old_div(lnx,2),y+old_div(lnx,2),z+old_div(lnx,2)))
				volfilt.insert_scaled_sum(v2m,(x+old_div(lnx,2),y+old_div(lnx,2),z+old_div(lnx,2)))
				if options.outfilte!=None : 
					volfilte.insert_scaled_sum(v1m,(x+old_div(lnx,2),y+old_div(lnx,2),z+old_div(lnx,2)))
				if options.outfilto!=None : 
					volfilto.insert_scaled_sum(v2m,(x+old_div(lnx,2),y+old_div(lnx,2),z+old_div(lnx,2)))
					
				volnorm.insert_scaled_sum(avgmask,(x+old_div(lnx,2),y+old_div(lnx,2),z+old_div(lnx,2)))
			
	if options.verbose>1: print("\nAll threads complete")

	# while the size of avgmask was selected to produce a nearly normalized image without further work
	# there were minor artifacts. The normalization deals with this.
	volnorm.process_inplace("math.reciprocal")
	volfilt.mult(volnorm)
	
	resvol.write_image("resvol.hdf")
	resvol143.write_image(options.output)
	volfilt.write_image(options.outfilt)

	if options.outfilte!=None : 
		volfilte.mult(volnorm)
		volfilte.write_image(options.outfilte)
	
	if options.outfilto!=None : 
		volfilto.mult(volnorm)
		volfilto.write_image(options.outfilto)
	
	#out=open("fsc.curves.txt","w")
	#out.write("# This file contains individual FSC curves from e2fsc.py. Only a fraction of computed curves are included.\n")
	#for y in fys:
		#print(len(y),y)
	#if len(fys)>100 : 
		#step=old_div(len(fys),100)
		#print("Saving 1/%d of curves to fsc.curves.txt + %d"%(step,len(funny)))
	#else: 
		#step=1
		#print("Saving all curves to fsc.curves.txt")
	
	#for i,x in enumerate(fx):
		#out.write( "%f\t"%x)
		#for j in range(0,len(fys),step):
			#out.write( "%f\t"%fys[j][i])
		
		## Also save any particularly low resolution curves
		#for j in funny:
			#out.write( "%f\t"%fys[j][i])
			
		#out.write("\n")
		
	#if len(funny)>1 :
		#print("WARNING: %d/%d curves were evaluated as being >0.5 AT Nyquist. While these values have been set to \
		#Nyquist (the maximum resolution for your sampling), these values are not meaningful, and could imply \
		#insufficient sampling, or bias in the underlying reconstructions."%(len(funny),len(fys)))

	E2end(logid)

if __name__ == "__main__":
        main()


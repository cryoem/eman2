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

def calc_oneres(jsd,vol1f,vol2f,apix,freq,ftsize,tophat=False,cutoff=0.143,rmask=None):
	"""Calculates a local correlation map and locally weighted filtered volume centered on a single spatial frequency
	this would normally be called multiple times (for each spatial frequency) and the results merged. rmask, if provided
	will also compute a correlation value under a single masked region. This could be done separately, but it is more 
	efficient when done here"""
	global curvol,minvol
	nx,ny,nz=vol1f["nx"],vol1f["ny"],vol1f["nz"]
	band=EMData(nx,ny,nz)
	band.set_complex(1)
	band.set_ri(1)
	
	# bandpass volume we apply as a filter to each image FFT
#	band.process_inplace("testimage.fourier.gaussianband",{"center":freq,"width":sqrt(2.0)})
	band.process_inplace("testimage.fourier.gaussianband",{"center":freq,"width":3})
	vol1b=(vol1f*band).do_ift()
	vol2b=(vol2f*band).do_ift()
	
	# We are computing a normalized dot product over a region, so we start with the squared and cross map
	volcor=vol1b*vol2b		# A*B
	vol1b2=vol1b.process("math.squared")
	vol2b2=vol2b.process("math.squared")
	
	corunmask=volcor["mean"]/sqrt(vol1b2["mean"]*vol2b2["mean"])
	if rmask!=None:
		rdot=(volcor*rmask)["mean"]
		rv1=(vol1b2*rmask)["mean"]
		rv2=(vol2b2*rmask)["mean"]
		cormask=rdot/sqrt(rv1*rv2)
	else: cormask=None
	
	# then we low-pass filter (local convolution) to select the "region" for each vector
	volcor.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	vol1b2.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	vol2b2.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	vol1b2.process_inplace("threshold.belowtozero",{"minval":0})   # obviously this should not be necessary, but (likely due to FFT roundoff error) it is
	vol2b2.process_inplace("threshold.belowtozero",{"minval":0})

	
	# put it all together to get a normalized local correlation at this one frequency
	vol1b2.mult(vol2b2)
	vol1b2.process_inplace("math.sqrt")
	vol1b2.process_inplace("math.reciprocal",{"zero_to":0.0})
	volcor.mult(vol1b2)
	
	if tophat:
		# Now let's turn that into a binary filter (tophat)
		filt=volcor.process("threshold.binary",{"value":cutoff})
	else:
		# Now let's turn that into a Wiener filter
		if freq==3: 
			minvol=volcor.copy()
			curvol=3
		elif freq>3:
			while curvol<freq-1: time.sleep(0.1)		# we wait until we have all lower frequency calculations done to insure monotonic decrease
			minvol.update_min(volcor)
#			volcor.update_min(minvol)	# could just copy, but timing is likely about the same
			# if we find a low resolution patch where the FSC falls below zero, but then rises above 0.5, we reset our below zero threshold
			msk=volcor.process("threshold.binary",{"value":.5})
			minvol.add(msk)
			# Then we zero the FSC in any regions which have previously fallen below zero to avoid spurious high resolution info in 
			msk=minvol.process("threshold.binary",{"value":.01})
			volcor.mult(msk)
			msk=None
			
			curvol=max(freq,curvol)		# should always be freq, but just to be safe...
		
		# Setting scalesnr here is slightly under-filtering the Wiener filter. This is emperical to help permit convergence
		# of higher frequency information
		filt=volcor.process("math.ccc_snr_wiener",{"wiener":1,"scalesnr":3.0})

	filt1=vol1b*filt
	filt2=vol2b*filt
	
	jsd.put((freq,filt1,filt2,volcor,corunmask,cormask))
	
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
	parser.add_argument("--compressbits", type=int,help="Bits to keep when writing volumes with compression. 0->lossless floating point. Default 10 (3 significant figures)", default=10)
	parser.add_argument("--localsizea", type=int, help="Size in Angstroms of the local region to compute the resolution in",default=50)
	parser.add_argument("--apix", type=float, help="A/pix to use for the comparison (default uses Vol1 apix)",default=0)
	parser.add_argument("--normin",type=str,help="Apply a real space normalization to each input before FSC. Default normalize.edgemean. Use 'none' to disable.",default="normalize.edgemean")
	parser.add_argument("--cutoff", type=float, help="fsc cutoff. default is 0.143",default=0.143)
	parser.add_argument("--mask",type=str,help="Optional mask to produce masked overall FSC curve. Must have the same dimensions as the input volumes.",default=None)
	parser.add_argument("--tophat",action="store_true",help="If set, the local filter is a tophat filter, otherwise a local Wiener filter is applied",default=False)
	parser.add_argument("--sampfscs",action="store_true",help="If set, full fsc curves are stored for a range of specific locations within the volume",default=False)
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
	
	if options.normin!=None and options.normin.lower()!="none" and len(options.normin)>3 :
		v1.process_inplace(options.normin)
		v2.process_inplace(options.normin)
		
	#if options.mask!=None:
		#mask=EMData(options.mask)
		#v1.mult(mask)
		#v2.mult(mask)
	
	if options.apix>0 : apix=options.apix
	else :
		apix=v1["apix_x"]
		print("Using %1.2f A/pix"%apix)
	
	# make sure apix set for all volumes
	v1["apix_x"],v1["apix_y"],v1["apix_z"]=apix,apix,apix
	v2["apix_x"],v2["apix_y"],v2["apix_z"]=apix,apix,apix
	
	nx,ny,nz=v1["nx"],v1["ny"],v1["nz"]
	print("%d x %d x %d"%(nx,ny,nz))
	if nx!=ny or nx!=nz : print("Warning: non-cubic volumes may produce unexpected results")

	box=good_size(ny+options.localsizea//apix)
	print("Using box-size: ",box)
	v1f=v1.get_clip(Region((nx-box)/2,(ny-box)/2,(nz-box)/2,box,box,box)).do_fft()
	v2f=v2.get_clip(Region((nx-box)/2,(ny-box)/2,(nz-box)/2,box,box,box)).do_fft()
	if options.mask!=None:
		mask=EMData(options.mask).get_clip(Region((nx-box)/2,(ny-box)/2,(nz-box)/2,box,box,box))
	else: mask=None
	
	if options.verbose: print("Preparing for local calculation")

	# This averager will contain the resolution map when done
	resvola=Averagers.get("minmax",{"max":1})
	
	# These averagers will contain the final filtered volumes
	filtvol1=Averagers.get("mean")
	filtvol2=Averagers.get("mean")
	
	# These are used to guarantee sequence and monotonic reduction after a certain threshold in Wiener filtration
	global curvol,minvol
	curvol=0
	minvol=None

	# used for overall FSC curves
	fscum=[0]*(box//2)
	fscm=[0]*(box//2)
	fscs=[r/(box*apix) for r in range(1,box//2)]

	NTHREADS=max(options.threads,2)		# we have one thread just writing results
	jsd=queue.Queue(0)
	thrds=[(jsd,v1f,v2f,apix,f,options.localsizea,options.tophat,options.cutoff,mask) for f in range(1,box//2)]

	if options.sampfscs:
		#list of locations to sample full FSC curves, relative to middle of box
		samplocs=((0,0,0),(10,0,0),(0,0,10),(50,0,0),(0,0,50),(10,25,25),(0,box//4,0),(-box//6,-box//6,-box//6))
		samplocs=[(i[0]+box//2,i[1]+box//2,i[2]+box//2) for i in samplocs]
		samps=[[0]*(box//2) for i in range(len(samplocs))]

	# here we run the threads and save the results, no actual alignment done here
	if options.verbose: print(len(thrds)," threads")
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose>1 : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch]=threading.Thread(target=calc_oneres,args=thrds[thrtolaunch])
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
		if options.verbose>1 and thrtolaunch>0:
			frac=thrtolaunch/float(len(thrds))
			print("{:1.1f}% complete".format(100.0*frac))
	
		while not jsd.empty():
			freq,filt1,filt2,volcor,corunmask,cormask=jsd.get()
			filt1=filt1.get_clip(Region((box-nx)/2,(box-ny)/2,(box-nz)/2,nx,ny,nz))
			filt2=filt2.get_clip(Region((box-nx)/2,(box-ny)/2,(box-nz)/2,nx,ny,nz))
			fscum[freq-1]=corunmask
			fscm[freq-1]=cormask
			volcor=volcor.get_clip(Region((box-nx)/2,(box-ny)/2,(box-nz)/2,nx,ny,nz))
			filtvol1.add_image(filt1)
			filtvol2.add_image(filt2)
			# Threshold the map representing a single spatial frequency, then 
			# replace the value with spatial freq. Max of all of the individual resolution maps
			# should represent the local resolution map
			#volcor.get_clip(Region(0,0,nz*2//3,nx,ny,1)).write_image("res_slice.hdf",freq-1)
			if options.sampfscs:
				for i,s in enumerate(samplocs):
					samps[i][freq]=volcor[s]
			
			volcor.process_inplace("threshold.binary",{"value":0.143})
			volcor.mult(freq/(ny*apix))
			resvola.add_image(volcor)

	for t in thrds:
		t.join()

	if options.sampfscs:
		for i,s in enumerate(samplocs):
			out=open(f"fsc-{i}_{s[0]}_{s[1]}_{s[2]}.txt","w")
			for j in range(len(samps[i])):
				out.write(f"{j}\t{samps[i][j]}\n")
	
	av1=filtvol1.finish()
	if options.outfilte!=None: av1.write_compressed(options.outfilte,0,options.compressbits,erase=True)
	av2=filtvol2.finish()
	if options.outfilto!=None: av2.write_compressed(options.outfilto,0,options.compressbits,erase=True)
	((av1+av2)/2).write_compressed(options.outfilt,0,options.compressbits,erase=True)
	resvol=resvola.finish()
	resvol.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":2/options.localsizea})
	# compression mode 0 here not because of need for precision, but so Chimera gets the correct values
	resvol.write_compressed(options.output,0,0,erase=True)

	Util.save_data(1/(box*apix),1/(box*apix),fscum,"fsc_unmasked.txt")
	if mask!=None: Util.save_data(1/(box*apix),1/(box*apix),fscm,"fsc_masked.txt")
	
	E2end(logid)

if __name__ == "__main__":
        main()


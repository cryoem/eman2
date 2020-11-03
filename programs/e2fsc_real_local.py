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

def calc_oneres(jsd,vol1f,vol2f,freq,ftsize):
	nx,ny,nz=vol1f["nx"],vol1f["ny"],vol1f["nz"]
	band=EMData(nx,ny,nz)
	band.set_complex(1)
	band.set_ri(1)
	
	# bandpass volume we apply as a filter to each image FFT
	band.process_inplace("testimage.fourier.gaussianband",{"center":freq,"width":sqrt(2.0)})
	vol1b=(vol1f*band).do_ift()
	vol2b=(vol2f*band).do_ift()
	
	# We are computing a normalized dot product over a region, so we start with the squared and cross map
	volcor=vol1b*vol2b		# A*B
	vol1b.process_inplace("math.squared")
	vol2b.process_inplace("math.squared")
	
	# then we low-pass filter (local convolution) to select the "region" for each vector
	volcor.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	vol1b.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	vol2b.process_inplace("filter.lowpass.gauss",{"cutoff_resolv":1.0/ftsize})
	
	# put it all together to get a normalized local correlation at this one frequency
	vol1b.mult(vol2b)
	vol1b.process_inplace("math.sqrt")
	vol1b.process_inplace("math.reciprocal")
	volcor.mult(vol1b)
	
	# Now let's turn that into a Wiener filter
	filt=volcor.process("math.ccc_snr_wiener",{"wiener":1})
	filtav=((vol1f+vol2f)*filt).do_ift()
	
	jsd.put((filtav,volcor))
	
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
	parser.add_argument("--localsize", type=int, help="Size in pixels of the local region to compute the resolution in",default=16)
	parser.add_argument("--apix", type=float, help="A/pix to use for the comparison (default uses Vol1 apix)",default=0)
	parser.add_argument("--cutoff", type=float, help="fsc cutoff. default is 0.143",default=0.143)
#	parser.add_argument("--mask",type=str,help="Mask to apply to both input images before calculation",default=None)
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
	#if options.mask!=None:
		#mask=EMData(options.mask)
		#v1.mult(mask)
		#v2.mult(mask)
	
	if options.apix>0 : apix=options.apix
	else :
		apix=v1["apix_x"]
		print("Using %1.2f A/pix"%apix)
	
	nx,ny,nz=v1["nx"],v1["ny"],v1["nz"]
	print("%d x %d x %d"%(nx,ny,nz))
	if nx!=ny or nx!=nz : print("Warning: non-cubic volumes may produce unexpected results")

	

	if apix*lnx/2.0<10.0 :
		print("WARNING: Local sampling box is <10 A. Adjusting to 16 A.")
		lnx=int(floor(32.0/apix))
	print("Local region is %d pixels"%lnx)
	
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
	
	resvol=EMData(nx,ny,nz)
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

	NTHREADS=max(options.threads,2)		# we have one thread just writing results
	jsd=queue.Queue(0)
	thrds=[(jsd,) for i in range(NTHREADS-1)]
	
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
			print("{}% complete".format(100.0*frac))
	
		while not jsd.empty():
			fsp,nds=jsd.get()
			for n,d in nds:
				angs[(fsp,n)]=d
				if options.saveali:
					v=EMData(fsp,n)
					v.transform(d["xform.align2d"])
					v.write_image("{}/aliptcls_{:02d}.hdf".format(options.path,options.iter),n)


	for t in thrds:
		t.join()

	
	resvol.write_image("resvol.hdf")
	resvol143.write_image(options.output)
	volfilt.write_image(options.outfilt)

	E2end(logid)

if __name__ == "__main__":
        main()


#!/usr/bin/env python
#
# Author: Steven Ludtke  4/28/2024
# Copyright (c) 2023- Baylor College of Medicine
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

from EMAN3 import *
from EMAN3tensor import *
import numpy as np
import sys
import time

def main():

	usage="""e3make3d_gauss.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--volout", type=str,help="Volume output file", default="threed.hdf")
	parser.add_argument("--gaussout", type=str,help="Gaussian list output file",default=None)
	parser.add_argument("--volfilt", type=float, help="Lowpass filter to apply to output volume, absolute, Nyquist=0.5", default=0.3)
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--savesteps", action="store_true",help="Save the gaussian parameters for each refinement step, for debugging and demos")
	parser.add_argument("--tomo", action="store_true",help="tomogram mode, changes optimization steps")

	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	tf_set_device(dev=0,maxmem=options.gpuram)

	llo=E3init(sys.argv)

	nptcl=EMUtil.get_image_count(args[0])
	nxraw=EMData(args[0],0,True)["nx"]

	if options.savesteps: 
		try: os.unlink("steps.hdf")
		except: pass

	if options.verbose: print(f"{nptcl} particles at {nxraw}^3")

	# definition of downsampling sequence for stages of refinement
	# #ptcl, downsample, iter, frc weight, amp threshold, replicate, step coef
	# replication skipped in final stage
	if options.tomo:
		stages=[
			[256,32,  16,1.8,  0,1,.01, 3.0],
			[256,32,  16,1.8,0.3,2,.01, 1.0],
			[256,64,  24,1.5,  0,1,.005,1.0],
			[256,64,  48,1.5,0.2,2,.005,0.5],
			[256,128, 32,1.2,0.2,3,.003,0.5],
			[256,512, 32,1.0,0.2,3,.001,0.25],
			[256,1024,24,0.8,0.2,1,.001,0.1]
		]
	else:
		stages=[
			[500,   16,16,1.8,-3  ,1,.01, 2.0],
			[500,   16,16,1.8, 0  ,4,.01, 1.0],
			[1000,  32,16,1.5,-1  ,1,.005,1.5],
			[1000,  32,16,1.5,-1.5,3,.005,1.0],
			[2500,  64,24,1.2,-1.5,3,.003,1.0],
			[10000,256,24,1.0,-2  ,2,.001,0.5],
			[25000,512,12,0.8,-2  ,1,.001,0.2]
		]

	times=[time.time()]

	# Cache initialization
	if options.verbose: print("Caching particle data")
	downs=sorted(set([s[1] for s in stages]))
	caches={down:StackCache(f"tmp_{down}.cache",nptcl) for down in downs} 	# dictionary keyed by box size
	for i in range(0,nptcl,2500):
		if options.verbose>1: print(f"Caching {i}/{nptcl}")
		stk=EMStack2D(EMData.read_images(args[0],range(i,min(i+2500,nptcl))))
		orts,tytx=stk.orientations
		tytx/=nxraw
		stkf=stk.do_fft()
		for down in downs:
			stkfds=stkf.downsample(min(down,nxraw))
			caches[down].write(stkfds,i,orts,tytx)

	gaus=Gaussians()
	#Initialize Gaussians to random values with amplitudes over a narrow range
	rnd=tf.random.uniform((options.initgauss,4))     # specify the number of Gaussians to start with here
	rnd+=(-.5,-.5,-.5,10.0)
	if options.tomo: gaus._data=rnd/(.9,.9,.9,10.0)	# amplitudes set to ~1.0, positions random within 2/3 box size
	else: gaus._data=rnd/(1.5,1.5,1.5,100.0)	# amplitudes set to ~1.0, positions random within 2/3 box size

	times.append(time.time())
	ptcls=[]
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")
#
# 		# stage 1 - limit to ~1000 particles for initial low resolution work
# 		if options.verbose: print(f"\tReading Files {min(stage[0],nptcl)} ptcl")
# 		ptcls=EMStack2D(EMData.read_images(args[0],range(0,nptcl,max(1,nptcl//stage[0]))))
# 		orts,tytx=ptcls.orientations
# 		tytx/=nxraw
# 		ptclsf=ptcls.do_fft()
#
# 		if options.verbose: print(f"\tDownsampling {min(nxraw,stage[1])} px")
# 		if stage[1]<nxraw: ptclsfds=ptclsf.downsample(stage[1])    # downsample specifies the final size, not the amount of downsampling
# 		else: ptclsfds=ptclsf			# if true size is smaller than stage size, don't downsample, obviously
# 		ny=stage[1]
# 		ptcls=None		# free resouces since we're committed to re-reading files for now
# 		ptclsf=None
# #		ny=ptclsfds.shape[1]

		nliststg=range(sn,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity

#	print(ptclsfds.shape,tytx.shape)
		
		if options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tshift_grad\tamp_grad")
		for i in range(stage[2]):		# training epochs
			for j in range(0,len(nliststg),500):	# compute the gradient step piecewise due to memory limitations, 1000 particles at a time
				ptclsfds,orts,tytx=caches[stage[1]].read(nliststg[j:j+500])
				step0,qual0,shift0,sca0=gradient_step(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
				if j==0:
					step,qual,shift,sca=step0,qual0,shift0,sca0
				else:
					step+=step0
					qual+=qual0
					shift+=shift0
					sca+=sca0
			norm=len(nliststg)//500+1
			step/=norm
			qual/=norm
			shift/=norm
			sca/=norm
			gaus.add_tensor(step)
			if options.savesteps: from_numpy(gaus.numpy).write_image("steps.hdf",-1)

			print(f"{i}: {qual:1.4f}\t{shift:1.4f}\t\t{sca:1.4f}")

		# if options.savesteps:
		# 	vol=gaus.volume(nxraw)
		# 	vol.emdata[0].process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt})
		# 	vol.write_images(f"A_vol_opt_{sn}.hdf")

		# filter results and prepare for stage 2
		g0=len(gaus)
		gaus.norm_filter(sig=stage[4])
		g1=len(gaus)
		if stage[5]>0: gaus.replicate(stage[5],stage[6])
		g2=len(gaus)
		print(f"Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians  {local_datetime()}")
		times.append(time.time())
	
		# do this at the end of each stage in case of early termination
		if options.gaussout is not None:
			out=open(options.gaussout,"w")
			for x,y,z,a in gaus.tensor: out.write(f"{x:1.5f}\t{y:1.5f}\t{z:1.5f}\t{a:1.3f}\n")

		vol=gaus.volume(nxraw)
		#vol.emdata[0].process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt})
		vol.write_images(options.volout)


	times=np.array(times)
	times-=times[0]
	if options.verbose>1 : print(times)

	E3end(llo)

#@tf.function
def gradient_step(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0):
	"""Computes one gradient step on the Gaussian coordinates given a set of particle FFTs at the appropriate scale,
	computing FRC to axial Nyquist, with specified linear weighting factor (def 1.0). Linear weight goes from
	0-2. 1 is unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns step, qual, shift, scale
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]

	with tf.GradientTape() as gt:
		gt.watch(gaus.tensor)
		projs=gaus.project_simple(orts,ny,tytx=tytx)
		projsf=projs.do_fft()
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,3)	# specifying ny/2 radius explicitly so weight functions

	grad=gt.gradient(frcs,gaus._data)
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	shift=tf.math.reduce_std(grad[:,:3])	# translational std
	sca=tf.math.reduce_std(grad[:,3])		# amplitude std
	xyzs=relstep/(shift*500)   				# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*(xyzs,xyzs,xyzs,relstep/(sca*250))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))
#	print(f"{i}) {float(qual)}\t{float(shift)}\t{float(sca)}")


if __name__ == '__main__':
	main()

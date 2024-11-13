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
from EMAN3jax import *
import numpy as np
import sys
import time
import os

def main():

	usage="""e3make3d_gauss.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--volout", type=str,help="Volume output file", default="threed.hdf")
	parser.add_argument("--gaussout", type=str,help="Gaussian list output file",default=None)
	parser.add_argument("--volfiltlp", type=float, help="Lowpass filter to apply to output volume in A, 0 disables, default=40", default=40)
	parser.add_argument("--volfilthp", type=float, help="Highpass filter to apply to output volume in A, 0 disables, default=2500", default=2500)
	parser.add_argument("--apix", type=float, help="A/pix override for raw data", default=-1)
	parser.add_argument("--thickness", type=float, help="For tomographic data specify the Z thickness in A to limit the reconstruction domain", default=-1)
	parser.add_argument("--preclip",type=int,help="Trim the input images to the specified (square) box size in pixels", default=-1)
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--savesteps", action="store_true",help="Save the gaussian parameters for each refinement step, for debugging and demos")
	parser.add_argument("--tomo", action="store_true",help="tomogram mode, changes optimization steps")
#	parser.add_argument("--ctf", action="store_true",help="Includes ctf in the projections")
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--dfmin", type=float, help="The minimum defocus appearing in the project, for use with --ctf",default=0.5)
	parser.add_argument("--dfmax", type=float, help="The maximum defocus appearing in the project, for use with --ctf",default=2.0)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	jax_set_device(dev=0,maxmem=options.gpuram)

	llo=E3init(sys.argv)

	nptcl=EMUtil.get_image_count(args[0])
	nxraw=EMData(args[0],0,True)["nx"]
	if options.preclip>0: nxraw=options.preclip
	nxrawm2=good_size_small(nxraw-2)
	if options.apix>0: apix=options.apix
	else: apix=EMData(args[0],0,True)["apix_x"]
	if options.thickness>0: zmax=options.thickness/(apix*nxraw*2.0)		# instead of +- 0.5 Z range, +- zmax range
	else: zmax=0.5
	if options.ctf>0:
		if options.tomo:
			ctf=EMData(args[0],0,True)["ctf"].to_dict() # Assuming tomo uses the file from particles, created by extract particles
		else:
			js=js_open_dict(info_name(EMData(args[0],0,True)["ptcl_source_image"])) # Assuming SPR uses lst file ptcls_XX.lst created by spt refinement
			ctf=js["ctf"][0].to_dict()
			js.close()
		cs=ctf["cs"]
		voltage=ctf["voltage"]
		ampcont=ctf["ampcont"]
		dfrange=(options.dfmin,options.dfmax)
		# Create the ctf stack
		ctf_stack,dfstep=create_ctf_stack(dfrange,voltage,cs,ampcont,nxrawm2,apix)

	if options.verbose: print(f"Input data box size {nxraw}x{nxraw} at {apix} A/pix. Maximum downsampled size for refinement {nxrawm2}. Thickness limit {zmax}. {nptcl} input images")

	if options.savesteps: 
		try: os.unlink("steps.hdf")
		except: pass

	if options.verbose: print(f"{nptcl} particles at {nxraw}^3")

	# definition of downsampling sequence for stages of refinement
	# 0) #ptcl, 1) downsample, 2) iter, 3) frc weight, 4) amp threshold, 5) replicate, 6) repl. spread, 7) step coef
	# replication skipped in final stage
	if options.tomo:
		stages=[
			[256,32,  32,1.8, -1,1,.05, 3.0],
			[256,32,  32,1.8, -1,2,.05, 1.0],
			[256,64,  48,1.5, -2,1,.04,1.0],
			[256,64,  48,1.5, -2,2,.02,0.5],
			[256,128, 32,1.2, -3,4,.01,2],
			[256,256, 32,1.2, -2,1,.01,3],
			[256,512, 48,1.2, -2,1,.01,3],
			[256,1024,48,1.2, -3,1,.004,5]
		]
	else:
		stages=[
			[512,   16,16,1.8,-3  ,1,.03, 2.0],
			[512,   16,16,1.8, 0  ,4,.03, 1.0],
			[1024,  32,16,1.5, 0  ,4,.02,1.5],
			[1024,  32,16,1.5,-1  ,3,.02,1.0],
			[4096,  64,24,1.2,-1.5,3,.01,1.0],
			[8192, 256,24,1.0,-2  ,3,.005,1.0],
			[32768,512,12,0.8,-2  ,1,.002,0.75]
		]

	times=[time.time()]

	# Cache initialization
	if options.verbose: print("Caching particle data")
	downs=sorted(set([s[1] for s in stages]))
	caches={down:StackCache(f"tmp_{os.getpid()}_{down}.cache",nptcl) for down in downs} 	# dictionary keyed by box size
	for i in range(0,nptcl,1000):
		if options.verbose>1:
			print(f" Caching {i}/{nptcl}",end="\r",flush=True)
			sys.stdout.flush()
		stk=EMStack2D(EMData.read_images(args[0],range(i,min(i+1000,nptcl))))
		if options.preclip>0 : stk=stk.center_clip(options.preclip)
		orts,tytx=stk.orientations
		tytx/=(nxraw,nxraw,1) # Don't divide the defocus
		for im in stk.emdata: im.process_inplace("normalize.edgemean")
		stkf=stk.do_fft()
		for down in downs:
			stkfds=stkf.downsample(min(down,nxrawm2))
			caches[down].write(stkfds,i,orts,tytx)

	# Forces all of the caches to share the same orientation information so we can update them simultaneously below (FRCs not jointly cached!)
	for down in downs[1:]:
		caches[down].orts=caches[downs[0]].orts
		caches[down].tytx=caches[downs[0]].tytx

	if options.verbose>1: print("")

	gaus=Gaussians()
	#Initialize Gaussians to random values with amplitudes over a narrow range
	rng = np.random.default_rng()
	rnd=rng.uniform(0.0,1.0,(options.initgauss,4))		# start with completely random Gaussian parameters
#	rnd=tf.random.uniform((options.initgauss,4))     # specify the number of Gaussians to start with here
	rnd+=(-.5,-.5,-.5,2.0)
	if options.tomo: gaus._data=rnd/(.9,.9,1.0/zmax,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size
	else: gaus._data=rnd/(1.5,1.5,1.5,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size

	times.append(time.time())
	ptcls=[]
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")

#		nliststg=range(sn,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity

		if options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tshift_grad\tamp_grad\timshift\tgrad_scale")
		lqual=-1.0
		rstep=1.0
		for i in range(stage[2]):		# training epochs
			if nptcl>stage[0]: idx0=sn+i
			else: idx0=0
			nliststg=range(idx0,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current epoch in the current stage, sn+i provides stochasticity
			imshift=0.0
			for j in range(0,len(nliststg),512):	# compute the gradient step piecewise due to memory limitations, 512 particles at a time
				ptclsfds,orts,tytx=caches[stage[1]].read(nliststg[j:j+512])
				# standard mode, optimize gaussian parms only
#				if not options.tomo or sn<2:
				if options.ctf==0:
					step0,qual0,shift0,sca0=gradient_step(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
					if j==0:
						step,qual,shift,sca=step0,qual0,shift0,sca0
					else:
						step+=step0
						qual+=qual0
						shift+=shift0
						sca+=sca0
				elif options.ctf==2:
					dsapix=apix*nxraw/ptclsfds.shape[1]
					step0,qual0,shift0,sca0=gradient_step_layered_ctf(gaus,ptclsfds,orts,ctf_stack.downsample(ptclsfds.shape[1]),tytx,dfrange,dfstep,dsapix,stage[3],stage[7])
					if j==0:
						step,qual,shift,sca=step0,qual0,shift0,sca0
					else:
						step+=step0
						qual+=qual0
						shift+=shift0
						sca+=sca
				elif options.ctf==1:
					step0,qual0,shift0,sca0=gradient_step_ctf(gaus,ptclsfds,orts,ctf_stack.downsample(ptclsfds.shape[1]),tytx,dfrange,dfstep,stage[3],stage[7])
					if j==0:
						step,qual,shift,sca=step0,qual0,shift0,sca0
					else:
						step+=step0
						qual+=qual0
						shift+=shift0
						sca+=sca
				# optimize gaussians and image shifts
				else:
					step0,stept0,qual0,shift0,sca0,imshift0=gradient_step_tytx(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
					if j==0:
						step,stept,qual,shift,sca,imshift=step0,stept0,qual0,shift0,sca0,imshift0
						caches[stage[1]].add_orts(nliststg[j:j+512],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
					else:
						step+=step0
						caches[stage[1]].add_orts(nliststg[j:j+512],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
						qual+=qual0
						shift+=shift0
						sca+=sca0
						imshift+=imshift0
			norm=len(nliststg)//512+1
			qual/=norm
			if qual<lqual: rstep/=2.0	# if we start falling or oscillating we reduce the step within the epoch
			step*=rstep/norm
			shift/=norm
			sca/=norm
			imshift/=norm
			gaus.add_array(step)
			lqual=qual
			if options.savesteps: from_numpy(gaus.numpy).write_image("steps.hdf",-1)

			print(f"{i}: {qual:1.5f}\t{shift:1.5f}\t\t{sca:1.5f}\t{imshift:1.5f}\t{rstep:1.5f}")
			if qual>0.99: break

		# end of epoch, save images and projections for comparison
		if options.verbose>3:
			dsapix=apix*nxraw/ptclsfds.shape[1]
			ctf_stackds=ctf_stack.downsample(ptclsfds.shape[1])
			projs=gaus.project_simple(orts,ptclsfds.shape[1],tytx=tytx)
			ctf_projs=gaus.project_layered_ctf(orts,ctf_stackds,ptclsfds.shape[1],dsapix,dfrange,dfstep,tytx=tytx)
			single_ctf_projs=gaus.project_ctf(orts,ctf_stackds,ptclsfds.shape[1],dfrange,dfstep,tytx=tytx)
			transforms=orts.transforms(tytx)
#			# Need to calculate the ctf corrected projection then write 1. particle 2. simple projection 3. corrected simple projection 4.ctf projection
			ptclds=ptclsfds.do_ift()
			for i in range(len(projs)):
				a=ptclds.emdata[i]
				b=projs.emdata[i]
				c=single_ctf_projs.emdata[i]
				d=ctf_projs.emdata[i]
				a["apix_x"]=dsapix
				a["apix_y"]=dsapix
				b["apix_x"]=dsapix
				b["apix_y"]=dsapix
				c["apix_x"]=dsapix
				c["apix_y"]=dsapix
				d["apix_x"]=dsapix
				d["apix_y"]=dsapix
				a.process_inplace("normalize")
				b.process_inplace("filter.matchto",{"to":a})
#				a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2)
#				b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2+1)
				c.process_inplace("filter.matchto",{"to":a})
				d.process_inplace("filter.matchto",{"to":a})
				a["xform.projection"]=transforms[i]
				b["xform.projection"]=transforms[i]
				c["xform.projection"]=transforms[i]
				d["xform.projection"]=transforms[i]
				a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4)
				b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+1)
				c.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+2)
				d.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+3)


		# if options.savesteps:
		# 	vol=gaus.volume(nxraw,zmax)
		# 	vol.emdata[0].process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt})
		# 	vol.write_images(f"A_vol_opt_{sn}.hdf")

		# filter results and prepare for stage 2
		g0=len(gaus)
		if options.tomo: gaus.norm_filter(sig=stage[4])		# gaussians outside the box may be important!
		else: gaus.norm_filter(sig=stage[4],rad_downweight=0.33)
		g1=len(gaus)
		if stage[5]>0: gaus.replicate(stage[5],stage[6])
		g2=len(gaus)
		print(f"Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians  {local_datetime()}")
		times.append(time.time())
	
		# do this at the end of each stage in case of early termination
		if options.gaussout is not None and g2 != 0:
			np.savetxt(options.gaussout,gaus.numpy,fmt="%0.4f",delimiter="\t")
			# out=open(options.gaussout,"w")
			# for x,y,z,a in gaus.tensor: out.write(f"{x:1.5f}\t{y:1.5f}\t{z:1.5f}\t{a:1.3f}\n")

		# show individual shifts at high verbosity
		if options.verbose>2:
			print("TYTX: ",(caches[stage[1]].tytx*nxraw).astype(np.int32))

	outsz=min(1024,nxraw)
	times.append(time.time())
	vol=gaus.volume(outsz,zmax).emdata[0]
	times.append(time.time())
	vol["apix_x"]=apix*nxraw/outsz
	vol["apix_y"]=apix*nxraw/outsz
	vol["apix_z"]=apix*nxraw/outsz
	vol.write_image(options.volout.replace(".hdf","_unfilt.hdf"),0)
	if options.volfilthp>0: vol.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.volfilthp})
	if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
	times.append(time.time())
	vol.write_image(options.volout,0)


	times=np.array(times)
	#times-=times[0]
	times=times[1:]-times[:-1]
	if options.verbose>1 : print(times.astype(np.int32))

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
	mx=orts.to_mx2d(swapxy=True)
	gausary=gaus.jax
	ptcls=ptclsfds.jax
#	print("mx ",mx.shape)

	frcs,grad=gradvalfn(gausary,mx,tytx,ptcls,weight)

	qual=frcs.mean()			# this is the average over all projections, not the average over frequency
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*500)   	# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*jnp.array((xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))
#	print(f"{i}) {float(qual)}\t{float(shift)}\t{float(sca)}")

def prj_frc(gausary,mx2d,tytx,ptcls,weight):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_simple_fn(gausary,mx2d,ny,tytx)
	return jax_frc_jit(jax_fft2d(prj),ptcls,weight,2)

gradvalfn=jax.value_and_grad(prj_frc)

def gradient_step_tytxccf(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0):
	"""Computes one gradient step on the Gaussian coordinates and image shifts given a set of particle FFTs at the appropriate scale,
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
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,2)	# specifying ny/2 radius explicitly so weight functions

	grad,gradtytx=gt.gradient(frcs,(gaus._data,tytx))
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	shift=tf.math.reduce_std(grad[:,:3])	# translational std
	imshift=tf.math.reduce_std(gradtytx)	# image shift std
	sca=tf.math.reduce_std(grad[:,3])		# amplitude std
	xyzs=relstep/(shift*500)   				# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*(xyzs,xyzs,xyzs,relstep/(sca*250))	# amplitude scale, 500 heuristic, TODO: may change
	tytxstep=gradtytx*relstep/(imshift*2000)
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,tytxstep,float(qual),float(shift),float(sca),float(imshift))


def gradient_step_tytx(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0):
	"""Computes one gradient step on the Gaussian coordinates and image shifts given a set of particle FFTs at the appropriate scale,
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
		gt.watch(tytx)
		projs=gaus.project_simple(orts,ny,tytx=tytx)
		projsf=projs.do_fft()
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,2)	# specifying ny/2 radius explicitly so weight functions

	grad,gradtytx=gt.gradient(frcs,(gaus._data,tytx))
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	shift=tf.math.reduce_std(grad[:,:3])	# translational std
	imshift=tf.math.reduce_std(gradtytx)	# image shift std
	sca=tf.math.reduce_std(grad[:,3])		# amplitude std
	xyzs=relstep/(shift*500)   				# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*(xyzs,xyzs,xyzs,relstep/(sca*250))	# amplitude scale, 500 heuristic, TODO: may change
	tytxstep=gradtytx*relstep/(imshift*2000)
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,tytxstep,float(qual),float(shift),float(sca),float(imshift))
#	print(f"{i}) {float(qual)}\t{float(shift)}\t{float(sca)}")

def gradient_step_ctf(gaus,ptclsfds,orts,ctf_stackds,tytx,dfrange,dfstep,weight=1.0,relstep=1.0):
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
		projs=gaus.project_ctf(orts,ctf_stackds,ny,dfrange,dfstep,tytx=tytx)
		projsf=projs.do_fft()
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,2)	# specifying ny/2 radius explicitly so weight functions

	grad=gt.gradient(frcs,gaus._data)
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	shift=tf.math.reduce_std(grad[:,:3])	# translational std
	sca=tf.math.reduce_std(grad[:,3])		# amplitude std
	xyzs=relstep/(shift*500)   				# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*(xyzs,xyzs,xyzs,relstep/(sca*250))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))

def gradient_step_layered_ctf(gaus,ptclsfds,orts,ctf_stackds,tytx,dfrange,dfstep,dsapix,weight=1.0,relstep=1.0):
	"""Computes one gradient step on the Gaussian coordinates given a set of particle FFTs at the appropriate scale,
	computing FRC to axial Nyquist, with specified linear weighting factor (def 1.0). Linear weight goes from
	0-2. 1 is unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns step, qual, shift, scale
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]

	with tf.GradientTape(persistent=True) as gt:
		gt.watch(gaus.tensor)
		projs=gaus.project_layered_ctf(orts,ctf_stackds,ny,dsapix,dfrange,dfstep,tytx=tytx)
		projsf=projs.do_fft() # TODO: Remove and have projection return fourier transform? It is already in fourier space there...
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,2)	# specifying ny/2 radius explicitly so weight functions

	grad=gt.gradient(frcs,gaus._data)
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	shift=tf.math.reduce_std(grad[:,:3])	# translational std
	sca=tf.math.reduce_std(grad[:,3])		# amplitude std
	xyzs=relstep/(shift*500)   				# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*(xyzs,xyzs,xyzs,relstep/(sca*250))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))



if __name__ == '__main__':
	main()

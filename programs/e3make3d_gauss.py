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
import jax
import optax
import numpy as np
import sys
import time
import os

# used for profiling. This can be commented out as long as the @profile line is also commented out
#from line_profiler import profile

try: os.mkdir(".jaxcache")
except: pass

# We cache the JIT compilation results to speed up future runs
jax.config.update("jax_compilation_cache_dir", "./.jaxcache")
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 2)
jax.config.update("jax_persistent_cache_enable_xla_caches", "xla_gpu_per_fusion_autotune_cache_dir")

jax.config.update("jax_default_matmul_precision", "float32")

# @profile
def main():

	usage="""e3make3d_gauss.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--volout", type=str,help="Volume output file. Note that volumes will be appended to an existing file", default="threed.hdf")
	parser.add_argument("--gaussout", type=str,help="Gaussian list output file",default=None)
	parser.add_argument("--volfiltlp", type=float, help="Lowpass filter to apply to output volume in A, 0 disables, default=40", default=40)
	parser.add_argument("--volfilthp", type=float, help="Highpass filter to apply to output volume in A, 0 disables, default=2500", default=2500)
	parser.add_argument("--frc_z", type=float, help="FRC Z threshold (mean-sigma*Z)", default=3.0)
	parser.add_argument("--apix", type=float, help="A/pix override for raw data", default=-1)
	parser.add_argument("--thickness", type=float, help="For tomographic data specify the Z thickness in A to limit the reconstruction domain", default=-1)
	parser.add_argument("--outbox",type=int,help="output boxsize, permitting over/undersampling (impacts A/pix)", default=-1)
	parser.add_argument("--preclip",type=int,help="Trim the input images to the specified (square) box size in pixels", default=-1)
	parser.add_argument("--postclip",type=int,help="Trim the output volumes to the specified (square) box size in pixels (no impact on A/pix)", default=-1)
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--savesteps", action="store_true",help="Save the gaussian parameters for each refinement step, for debugging and demos")
	parser.add_argument("--combineiters", type=int, help="Specify an additional number of iterations to add to the end of refinement, volume will use all Gaussian positions during these iterations", default=-1)
	parser.add_argument("--tomo", action="store_true",help="tomogram mode, changes optimization steps")
	parser.add_argument("--tomo_seqali", type=int,default=0,help="align each image in the tilt series to the adjacent image, starting with the center image and working outward. Specify region size in pixels in image center for alignment.")
	parser.add_argument("--cttomo", action="store_true",help="Continous tilt tomogram mode, changes optimization steps")
	parser.add_argument("--spt", action="store_true",help="subtomogram averaging mode, changes optimization steps")
	parser.add_argument("--quick", action="store_true",help="single particle mode with less thorough refinement, but faster results")
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--ptcl3d_id", type=str, help="only use 2-D particles with matching ptcl3d_id parameter (lst file/header, use + for range)",default=None)
	parser.add_argument("--class", dest="classid", type=int, help="only use 2-D particles with matching class parameter (lst file/header)",default=-1)
	parser.add_argument("--dfmin", type=float, help="Minimum defocus override, for use with --ctf",default=-1)
	parser.add_argument("--dfmax", type=float, help="Maximum defocus override, for use with --ctf",default=-1)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--fscdebug", type=str,help="Compute the FSC of the final map with a reference volume for debugging",default=None)
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--profile", action="store_true",help="Used for code development only, not routine use")
	parser.add_argument("--cachepath",type=str,help="path for storing the cached images, ideally on a high speed drive. Default='.'",default=".")
#	parser.add_argument("--precache",type=str,help="Rather than perform a reconstruction, only perform caching on the input file for later use. String is the folder to put the cache files in.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	jax_set_device(dev=0,maxmem=options.gpuram)
	llo=E3init(sys.argv,options.ppid)

	nptcl=EMUtil.get_image_count(args[0])
	if options.quick : nptcl=min(nptcl,8192+4096)

	if options.ptcl3d_id is not None:
		if args[0][-4:]!=".lst" : error_exit("--ptcl3d_id only works with .lst input files")
		try:
			options.ptcl3d_id=int(options.ptcl3d_id)
			rng=1
		except:
			rng=int(options.ptcl3d_id.split("+")[1])
			options.ptcl3d_id=int(options.ptcl3d_id.split("+")[0])
		lsx=LSXFile(args[0])
		selimg=[i for i in range(len(lsx)) if lsx[i][2]["ptcl3d_id"] in range(options.ptcl3d_id,options.ptcl3d_id+rng)]
		nptcl=len(selimg)
		lsx=None

	if options.classid>=0:
		if args[0][-4:]!=".lst" : error_exit("--class only works with .lst input files")
		lsx=LSXFile(args[0])
		selimg=[i for i in range(len(lsx)) if lsx[i][2]["class"]==options.classid]
		nptcl=len(selimg)
		lsx=None

	if options.profile:
		selimg=tuple(range(0,min(2050,nptcl)))
		nptcl=len(selimg)
		print("WARNING: profiling mode enabled. Actual gaussian map results will not be useful, used for development only!")

	nxraw=EMData(args[0],0,True)["nx"]
	if options.preclip>0: nxraw=options.preclip
	nxrawm2=good_size_small(nxraw-2)
	frc_Z=options.frc_z
	if options.apix>0: apix=options.apix
	else: apix=EMData(args[0],0,True)["apix_x"]
	if options.thickness>0: zmax=options.thickness/(apix*nxraw*2.0)		# instead of +- 0.5 Z range, +- zmax range
	else: zmax=0.5

	if options.verbose: print(f"Input data box size {nxraw}x{nxraw} at {apix} A/pix. Maximum downsampled size for refinement {nxrawm2}. Thickness limit +-{zmax}. {nptcl} input images")

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
			[256,32,  32,1.8, -1,4,.05, 1.0],
			[256,64,  48,1.5, -2,4,.04,1.0],
			[256,64,  48,1.5, -2,16,.02,0.5],
			[256,128, 32,1.2, -3,16,.01,2.0],
			[256,256, 32,1.2, -2,64,.01,3.0],
			[256,512, 48,1.2, -2,128,.01,3.0],
			[256,1024,48,1.2, -3,0,.004,5.0]
		]
	elif options.cttomo:
		stages=[
			[4096,32,  32,1.8, -1,1,.05, 3.0],
			[4096,32,  32,1.8, -1,4,.05, 1.0],
			[4096,64,  48,1.5, -2,4,.04,1.0],
			[4096,64,  48,1.5, -2,16,.02,0.5],
			[4096,128, 32,1.2, -3,16,.01,2.0],
			[4096,256, 32,1.2, -2,64,.01,3.0],
			[4096,512, 48,1.2, -2,128,.01,3.0],
			[4096,1024,48,1.2, -3,0,.004,5.0]
		]
	elif options.spt:
		stages=[
			[2**10,32,  32,1.8, -1,1,.05, 3.0],
			[2**10,32,  32,1.8, -1,2,.05, 1.0],
			[2**12,64,  48,1.5, -2,2,.04,1.0],
			[2**12,64,  48,1.5, -2,8,.02,0.5],
			[2**14,128, 32,1.2, -3,16,.01,2],
			[2**16,256, 16,1.2, -2,32,.01,3],
			[2**17,512, 8,1.2, -2,0,.01,3]
#			[2**20,1024,48,1.2, -3,,.004,5]
		]
	elif options.profile:
		stages=[
			[512,   16,16,1.8,-3  ,1,.01, 2.0],
			[1024,  32,16,1.5, 0  ,4,.005,1.5],
			[2048, 128,4 ,1.2,-1.5,0,.003,1.0]
#			[8192, 256,32,1.0,-2  ,3,.003,1.0],
#			[32768,512,32,0.8,-2  ,1,.001,0.75]
		]
	elif options.quick:
		stages=[
			[512,   16,24,1.8,-3  ,1,.01, 2.0],
			[512,   16,24,1.8, 0  ,4,.01, 1.0],
			[1024,  32,24,1.5, 0  ,8,.005,1.5],
			[4096,  64,24,1.2,-1.5,16,.003,1.0],
			[8192, 256,24,1.0,-2  ,0,.003,1.0],
		]
	else:
		stages=[
			[512,   16,32,1.8,-3  ,1,.01, 2.0],
			[512,   16,32,1.8, 0  ,4,.01, 1.0],
			[1024,  32,32,1.5, 0  ,4,.005,1.5],
			[1024,  32,32,1.5,-1  ,8,.005,1.0],
			[4096,  64,32,1.2,-1.5,16,.003,1.0],
			[16384, 256,32,1.0,-2 ,32,.003,1.0],
			[65536, 512,32,0.8,-2 ,0,.001,0.75]
		]

	# limit sampling to (at most) the box size of the raw data
	# we do this by altering stages to help with jit compilation
	for i in range(len(stages)):
		stages[i][1]=min(stages[i][1],nxrawm2)

	batchsize=192
	if options.combineiters>0:
		stages[-1][2]+=options.combineiters*24 # Increase the number of iterations so can save Gaussians
		final_gaus = Gaussians()
		final_gaus._data = []
	times=[time.time()]

	# Cache initialization
	if options.verbose: print(f"{local_datetime()}: Caching particle data")
	downs=sorted(set([s[1] for s in stages]))
	if options.ctf > 0:
		mindf = float('inf')
		maxdf = 0

	# critical for later in the program, this initializes the radius images for all of the samplings we will use
	for d in downs: rad_img_int(d)

	caches={down:StackCache(f"{options.cachepath}/tmp_{os.getpid()}_{down}.cache",nptcl) for down in downs} 	# dictionary keyed by box size
	if options.ptcl3d_id is not None and options.ptcl3d_id>=0 :
		if options.verbose>1:
			print(f" Caching {nptcl}")
		stk=EMStack2D(EMData.read_images(args[0],selimg))
		if options.preclip>0 : stk=stk.center_clip(options.preclip)
		orts,tytx=stk.orientations
		tytx/=jnp.array((nxraw,nxraw,1)) # Don't divide the defocus
		for im in stk.emdata: im-=im["mean"]
			#im.process_inplace("normalize.edgemean")
		stkf=stk.do_fft()
		for down in downs:
			stkfds=stkf.downsample(down)
			caches[down].write(stkfds,0,orts,tytx)
	else:
		for i in range(0,nptcl,1000):
			if options.verbose>1:
				print(f" Caching {i}/{nptcl}",end="\r",flush=True)
				sys.stdout.flush()
			stk=EMStack2D(EMData.read_images(args[0],range(i,min(i+1000,nptcl))))
			if options.preclip>0 : stk=stk.center_clip(options.preclip)
			if options.tomo and options.tomo_seqali!=0 :
				stk.center_align_seq(options.tomo_seqali)
				if options.verbose>3: stk.write_images("dbg_ali.hdf")
			orts,tytx=stk.orientations
			tytx/= jnp.array((nxraw,nxraw, 1))
			if options.ctf>0:
				mindf = min(mindf, float(jnp.min(tytx[:, 2])))
				maxdf = max(maxdf, float(jnp.max(tytx[:, 2])))
			for im in stk.emdata: im.process_inplace("normalize.edgemean")
			stkf=stk.do_fft()
			for down in downs:
				stkfds=stkf.downsample(down)
				caches[down].write(stkfds,i,orts,tytx)

	# Forces all of the caches to share the same orientation information so we can update them simultaneously below (FRCs not jointly cached!)
	for down in downs[1:]:
		caches[down].orts=caches[downs[0]].orts
		caches[down].tytx=caches[downs[0]].tytx

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
#		dfrange=(options.dfmin,options.dfmax)
		dfstep = apix*apix/100 # Rough approximation of the correct step that doesn't involve calculating the wavelength
		boxlen = apix*stages[-1][1]*sqrt(3) # stages[-1][1] is the largest downsampling for the particle
		df_buffer = (boxlen/2)/(dfstep*10000) + dfstep
		dfrange=(mindf - df_buffer, maxdf + df_buffer)
		if options.dfmin > 0 and options.dfmax > 0:
			dfrange=(options.dfmin, options.dfmax)
		# Create the ctf stack
		ctf_stack,dfstep=create_ctf_stack(dfrange,voltage,cs,ampcont,nxrawm2,apix)

	if options.verbose>1: print(f"\n{local_datetime()}: Refining")

	gaus=Gaussians()
	#Initialize Gaussians to random values with amplitudes over a narrow range
	rng = np.random.default_rng()
	rnd=rng.uniform(0.0,1.0,(options.initgauss*9//10,4))		# start with completely random Gaussian parameters
#	rnd=tf.random.uniform((options.initgauss,4))     # specify the number of Gaussians to start with here
	neg = rng.uniform(0.0, 1.0, (options.initgauss//10, 4))		# 10% of gaussians are negative WRT the background (zero)
	rnd+=(-.5,-.5,-.5,2.0)
	neg+=(-.5,-.5,-.5,-3.0)
	rnd = np.concatenate((rnd, neg))
	if options.tomo: gaus._data=rnd/(.9,.9,1.0/zmax,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size
	else: gaus._data=rnd/(1.5,1.5,1.5,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size

	times.append(time.time())
	ptcls=[]
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")
		if options.profile and sn==2 : jax.profiler.start_trace("jax_trace")

#		nliststg=range(sn,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity

		if options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tshift_grad\tamp_grad\timshift\tgrad_scale")
		lqual=-1.0
		rstep=1.0
		# TODO: Ok, this should really use one of the proper optimization algorithms available from the deep learning toolkits
		# this basic conjugate gradient gets the job done, but not very efficiently I suspect...
		optim = optax.adam(.003)		# parm is learning rate
#		optim = optax.lion(.003)		# tried, seems not quite as good as Adam in test, but maybe worth another try
#		optim = optax.lamb(.005)		# tried, slightly better than adam, worse than lion
#		optim = optax.fromage(.01)		# tried, not as good
		optim_state=optim.init(gaus._data)		# initialize with data
		for i in range(stage[2]):		# training epochs
			if rstep<.01: break		# don't continue if we've optimized well at this level
			if nptcl>stage[0]*2: idx0=sn+i
			else: idx0=0
			nliststg=range(idx0,nptcl,max(1,nptcl//stage[0]+1))		# all of the particles to use in the current epoch in the current stage, sn+i provides stochasticity
			imshift=0.0
			for j in range(0,len(nliststg)-10,batchsize):	# compute the gradient step piecewise due to memory limitations, batchsize particles at a time. The "-10" prevents very small residual batches from being computed
				ptclsfds,orts,tytx=caches[stage[1]].read(nliststg[j:j+batchsize])
				if len(orts)<5 :
					print("Abort tiny batch: ",len(nliststg),j,batchsize)
					continue
				# standard mode, optimize gaussian parms only
#				if not options.tomo or sn<2:
				if options.ctf==0:
					step0,qual0,shift0,sca0=gradient_step_optax(gaus,ptclsfds,orts,tytx,stage[3],stage[7],frc_Z)
					# TODO: These nan_to_num shouldn't be necessary. Not sure what is causing nans
					step0=jnp.nan_to_num(step0)
					shift0=jnp.nan_to_num(shift0)
					sca0=jnp.nan_to_num(sca0)
					if j==0:
						step,qual,shift,sca=step0,-qual0,shift0,sca0
					else:
						step+=step0
						qual-=qual0
						shift+=shift0
						sca+=sca0

					# update alignments of data to gaussian projections
					# if i>3:
					# 	dtytx=align_2d(gaus,orts,tytx,ptclsfds)
					# 	caches[stage[1]].add_orts(nliststg[j:j+batchsize],None,-dtytx)

					# if i==stage[2]-1:
					# # 	fscs0=jax_frc(jax_fft2d(gauss_project_simple_fn(gaus.jax,orts.to_mx2d(swapxy=True),ptclsfds.jax.shape[1],tytx)),ptclsfds.jax,-1,1.0,2)
					# 	out=open(f"fscs_{sn}.txt","a" if j>0 else "w")
					# 	for ii,fsc in enumerate(np.array(fscs0)): out.write(f"{nliststg[j+ii]:d}\t{fsc:0.5f}\n")
					# 	out.close()
				elif options.ctf==2:
					dsapix=apix*nxraw/ptclsfds.shape[1]
					step0,qual0,shift0,sca0=gradient_step_layered_ctf_optax(gaus,ptclsfds,orts,jax_downsample_2d(ctf_stack.jax,ptclsfds.shape[1]),tytx,dfrange,dfstep,dsapix,stage[3],stage[7],frc_Z)
					step0=jnp.nan_to_num(step0)
					if j==0:
						step,qual,shift,sca=step0,-qual0,shift0,sca0
					else:
						step+=step0
						qual-=qual0
						shift+=shift0
						sca+=sca
				elif options.ctf==1:
					step0,qual0,shift0,sca0=gradient_step_ctf_optax(gaus,ptclsfds,orts,jax_downsample_2d(ctf_stack.jax,ptclsfds.shape[1]),tytx,dfrange,dfstep,stage[3],stage[7],frc_Z)
					step0=jnp.nan_to_num(step0)
					if j==0:
						step,qual,shift,sca=step0,-qual0,shift0,sca0
					else:
						step+=step0
						qual-=qual0
						shift+=shift0
						sca+=sca
				# optimize gaussians and image shifts
				else:
					step0,stept0,qual0,shift0,sca0,imshift0=gradient_step_tytx(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
					step0=jnp.nan_to_num(step0)
					if j==0:
						step,stept,qual,shift,sca,imshift=step0,stept0,qual0,shift0,sca0,imshift0
						caches[stage[1]].add_orts(nliststg[j:j+batchsize],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
					else:
						step+=step0
						caches[stage[1]].add_orts(nliststg[j:j+batchsize],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
						qual+=qual0
						shift+=shift0
						sca+=sca0
						imshift+=imshift0
			norm=len(nliststg)//batchsize+1
			qual/=norm
			# # if the quality got worse, we take smaller steps, starting by stepping back almost to the last good step
			# if qual<lqual:
			# 	rstep/=2.0			# if we start falling or oscillating we reduce the step within the epoch
			# 	step=-lstep*.95		# new gradient doesn't matter, first we want to mostly undo the previous step
			# 	lstep*=.05
			# 	gaus.add_array(step)
			# 	if options.savesteps: from_numpy(gaus.numpy).write_image("steps.hdf",-1)
			# 	print(f"{i}: {qual:1.5f}\t     \t\t     \t     \t{rstep:1.5f} reverse")
			# 	continue
			# step*=rstep/norm
			# lstep=step
			# gaus.add_array(step)
			# lqual=qual
			shift/=norm
			sca/=norm
			imshift/=norm

			if isnan(shift) or isnan(sca) :
				if i==0:
					print("ERROR: nan on gradient descent, saving crash images and exiting")
					ptclsfds.do_ift().write_images("crash_lastb_images.hdf",0)
					out=open("crash_lastb_ortdydx.txt","w")
					for i in range(len(orts)):
						out.write(f"{orts[i][0]:1.6f}\t{orts[i][1]*1000:1.6f}\t{orts[i][2]*1000:1.6f}\t{tytx[i][0]*1000:1.2f}\t{tytx[i][1]*1000:1.2f} (/1000)\n")
					sys.exit(1)
				else:
					print("ERROR: encountered nan on gradient descent, skipping epoch. Image numbers saved to crash_img_S_E.lst")
					try: os.unlinK("crash_img.lst")
					except: pass
					out=LSXFile(f"crash_img_{sn}_{i}.lst")
					for ii in nliststg: out.write(-1,ii,args[0])
					out=None
					continue

			update, optim_state = optim.update(step, optim_state, gaus._data)
			gaus._data = optax.apply_updates(gaus._data, update)

			if options.savesteps: from_numpy(gaus.numpy).write_image("steps.hdf",-1)

			print(f"{i}\t{qual:1.5f}\t{shift*1000:1.6f}\t\t{sca*1000:1.6f}\t{imshift*1000:1.6f}  # /1000")

			if qual>0.99: break

			# Combine final few iterations to give the final volume
			if options.combineiters>0 and sn == len(stages)-1 and stages[sn][2] - i <= options.combineiters*24 and (i+1-(stages[sn][2]-options.combineiters*24))%24 == 0:
				if len(final_gaus) == 0: final_gaus._data = np.array(gaus._data)
				else: final_gaus._data = np.concatenate([final_gaus.numpy, np.array(gaus.jax)], axis=0)
				rng=np.random.default_rng()
				gaus.coerce_jax()
				std=1/(5*min(1024,nxraw))
				gaus._data+=rng.normal(0,std,gaus._data.shape)
#				outsz=min(1024,nxraw)
#				vol=final_gaus.volume_np(outsz,zmax).center_clip(outsz).emdata[0]
#				vol["apix_x"]=apix*nxraw/outsz
#				vol["apix_y"]=apix*nxraw/outsz
#				vol["apix_z"]=apix*nxraw/outsz
#				if options.volfilthp>0: vol.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.volfilthp})
#				if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
#				vol.process_inplace("normalize.edgemean")
#				vol.write_image("testing_combineiters_vol.hdf",-1)

		# end of epoch, save images and projections for comparison
		if options.verbose>3:
			dsapix=apix*nxraw/ptclsfds.shape[1]
			mx2d=orts.to_mx2d(swapxy=True)
			gausary=gaus.jax
			ny=ptclsfds.shape[1]
			projs=EMStack2D(gauss_project_simple_fn(gausary,mx2d,ny,tytx))
			if options.ctf>0:
				mx3d=orts.to_mx3d()
				ctfaryds=jax_downsample_2d(ctf_stack.jax,ny)
				ctf_projs=EMStack2D(gauss_project_ctf_fn(gausary,mx2d,ctfaryds,ny,dfrange[0],dfrange[1],dfstep,tytx))
				layered_ctf_projs=EMStack2D(gauss_project_layered_ctf_fn(gausary,mx3d,ctfaryds,ny,dfrange[0],dfrange[1],dfstep,dsapix,tytx))
			transforms=orts.transforms(tytx)
#			# Need to calculate the ctf corrected projection then write 1. particle 2. simple projection 3. corrected simple projection 4.ctf projection
			ptclds=ptclsfds.do_ift()
			for i in range(len(projs)):
				a=ptclds.emdata[i]
				b=projs.emdata[i]
				a["apix_x"]=dsapix
				a["apix_y"]=dsapix
				b["apix_x"]=dsapix
				b["apix_y"]=dsapix
				a["xform.projection"]=transforms[i]
				b["xform.projection"]=transforms[i]
				a.process_inplace("normalize")
				b.process_inplace("filter.matchto",{"to":a})
				if options.ctf>0:
					c=ctf_projs.emdata[i]
					d=layered_ctf_projs.emdata[i]
					c["apix_x"]=dsapix
					c["apix_y"]=dsapix
					d["apix_x"]=dsapix
					d["apix_y"]=dsapix
					c.process_inplace("filter.matchto",{"to":a})
					d.process_inplace("filter.matchto",{"to":a})
					c["xform.projection"]=transforms[i]
					d["xform.projection"]=transforms[i]
					a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4)
					b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+1)
					c.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+2)
					d.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+3)
				else:
					a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2)
					b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2+1)

		# if options.savesteps:
		# 	vol=gaus.volume(nxraw,zmax)
		# 	vol.emdata[0].process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt})
		# 	vol.write_images(f"A_vol_opt_{sn}.hdf")

		# filter results and prepare for stage 2
		if stage[5]>0:			# no filter/replicate in the last stage
			g0=len(gaus)
			if options.tomo: gaus.norm_filter(sig=stage[4], cyl_mask=0.45)		# gaussians outside the box may be important!
			else: gaus.norm_filter(sig=stage[4],rad_downweight=0.33)
			g1=len(gaus)
			# Replicate gaussians to produce a specified total number for each stage. Critical that these numbers
			# fall in a small set of specific N for JIT compilation
			gaus.replicate_abs(stage[5]*options.initgauss,stage[6])
			g2=len(gaus)
		else: g0=g1=g2=len(gaus)
		print(f"{local_datetime()}: Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians")
		times.append(time.time())

		# do this at the end of each stage in case of early termination
		if options.gaussout is not None and g2 != 0:
			np.savetxt(options.gaussout,gaus.numpy,fmt="%0.4f",delimiter="\t")
			# out=open(options.gaussout,"w")
			# for x,y,z,a in gaus.tensor: out.write(f"{x:1.5f}\t{y:1.5f}\t{z:1.5f}\t{a:1.3f}\n")

		# show individual shifts at high verbosity
#		if options.verbose>2:
#			print("TYTX: ",(caches[stage[1]].tytx*nxraw).astype(np.int32))

<<<<<<< HEAD
	if options.profile : jax.profiler.stop_trace()
=======
>>>>>>> projection

	if options.outbox>0: outsz=options.outbox
	else: outsz=min(1024,nxraw)
	times.append(time.time())
#	if options.combineiters>0:np.savetxt("testing_combine_iters.hdf", final_gaus.numpy, fmt="%0.4f", delimiter="\t") # For testing
	if options.combineiters>0 and options.postclip>0: vol = final_gaus.volume_np(outsz,zmax).center_clip(options.postclip)
	elif options.combineiters>0: vol=final_gaus.volume_np(outsz,zmax).center_clip(outsz)
	elif options.postclip>0 : vol=gaus.volume(outsz,zmax).center_clip(options.postclip)
	else : vol=gaus.volume(outsz,zmax).center_clip(outsz)
	vol=vol.emdata[0]
	if options.sym not in ("c1","C1","I","i"):
		if options.verbose>0 : print(f"Apply {options.sym} symmetry to map (not gaussians)")
		vol.process_inplace("xform.applysym",{"sym":options.sym})
	times.append(time.time())
	vol["apix_x"]=apix*nxraw/outsz
	vol["apix_y"]=apix*nxraw/outsz
	vol["apix_z"]=apix*nxraw/outsz
	if options.ptcl3d_id is not None and options.ptcl3d_id>=0 : vol["ptcl3d_id"]=options.ptcl3d_id
	vol.write_image(options.volout.replace(".hdf","_unfilt.hdf"),-1)
	if options.volfilthp>0: vol.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.volfilthp})
	if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
	vol.process_inplace("normalize.edgemean")
	times.append(time.time())
	vol.write_image(options.volout,-1)

	# this is just to save some extra processing steps
	if options.fscdebug is not None: 
		os.system(f'e2proc3d.py {options.volout.split(":")[0]} {options.volout.rsplit(".",1)[0]}_fsc.txt --calcfsc {options.fscdebug}')

	times=np.array(times)
	#times-=times[0]
	times=times[1:]-times[:-1]
	if options.verbose>1 : print(times.astype(np.int32))

	E3end(llo)

#@tf.function
def gradient_step(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0,frc_Z=3.0):
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

	frcs,grad=gradvalfn(gausary,mx,tytx,ptcls,weight,frc_Z)

#	qual=frcs.mean()			# this is the average over all projections, not the average over frequency
	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*100)   	# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*jnp.array((xyzs,xyzs,xyzs,relstep/(sca*100)))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))
#	print(f"{i}) {float(qual)}\t{float(shift)}\t{float(sca)}")

# @profile
def gradient_step_optax(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0,frc_Z=3.0):
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

	frcs,grad=gradvalfnl(gausary,mx,tytx,ptcls,weight,frc_Z)

	qual=frcs			# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational (gauss) std
	sca=grad[:,3].std()			# amplitude std

	return (grad,float(qual),float(shift),float(sca))

# @profile
def prj_frc_loss(gausary,mx2d,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations. Returns -frc since optax wants to minimize, not maximize"""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_simple_fn(gausary,mx2d,ny,tytx)
#	print(prj.shape,ptcls.shape,weight,frc_Z)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfnl=jax.jit(jax.value_and_grad(prj_frc_loss))

def prj_frc(gausary,mx2d,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_simple_fn(gausary,mx2d,ny,tytx)
	return jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfn=jax.value_and_grad(prj_frc)

def align_2d(gaus,orts,tytx,ptclsfds):
	ny=ptclsfds.shape[1]
	#ptcls=ptclsfds.jax
	mx=orts.to_mx2d(swapxy=True)
	prj=jax_fft2d(gauss_project_simple_fn(gaus.jax,mx,ny,tytx))	# FFT of gaussian projection for each particle
	return ptclsfds.align_translate(prj)/ny			# ccf between each particle and its projection


# TODO: This function is not updated to jax
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

# TODO: This function is not updated to jax
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

def gradient_step_ctf(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,weight=1.0,relstep=1.0,frc_Z=3.0):
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

	frcs,grad=gradvalfn_ctf(gausary,mx,ctfaryds,dfrange[0],dfrange[1],dfstep,tytx,ptcls,weight,frc_Z)

#	qual=frcs.mean()			# this is the average over all projections, not the average over frequency
	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*500)   	# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*jnp.array((xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))

def prj_frc_ctf(gausary,mx2d,ctfary,dfmin,dfmax,dfstep,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_ctf_fn(gausary,mx2d,ctfary,ny,dfmin,dfmax,dfstep,tytx)
	return jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfn_ctf=jax.value_and_grad(prj_frc_ctf)


def gradient_step_ctf_optax(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,weight=1.0,relstep=1.0,frc_Z=3.0):
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

	frcs,grad=gradvalfnl_ctf(gausary,mx,ctfaryds,dfrange[0],dfrange[1],dfstep,tytx,ptcls,weight,frc_Z)

	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*500)   	# xyz scale factor, 1000 heuristic, TODO: may change

	return (grad,float(qual),float(shift),float(sca))

def prj_frc_loss_ctf(gausary,mx2d,ctfary,dfmin,dfmax,dfstep,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	prj=gauss_project_ctf_fn(gausary,mx2d,ctfary,ny,dfmin,dfmax,dfstep,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfnl_ctf=jax.value_and_grad(prj_frc_loss_ctf)


def gradient_step_layered_ctf(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,dsapix,weight=1.0,relstep=1.0,frc_Z=3.0):
	"""Computes one gradient step on the Gaussian coordinates given a set of particle FFTs at the appropriate scale,
	computing FRC to axial Nyquist, with specified linear weighting factor (def 1.0). Linear weight goes from
	0-2. 1 is unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns step, qual, shift, scale
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]
	mx=orts.to_mx3d()
	gausary=gaus.jax
	ptcls=ptclsfds.jax

	frcs,grad=gradvalfn_layered_ctf(gausary,mx,ctfaryds,dfrange[0],dfrange[1],dfstep,dsapix,tytx,ptcls,weight, frc_Z)

#	qual=frcs.mean()			# this is the average over all projections, not the average over frequency
	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*500)   	# xyz scale factor, 1000 heuristic, TODO: may change
#	gaus.add_tensor(grad*(xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	step=grad*jnp.array((xyzs,xyzs,xyzs,relstep/(sca*250)))	# amplitude scale, 500 heuristic, TODO: may change
	#print(f"{qual}\t{shift}\t{sca}")

	return (step,float(qual),float(shift),float(sca))

def prj_frc_layered_ctf(gausary,mx3d,ctfary,dfmin,dfmax,dfstep,apix,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_layered_ctf_fn(gausary,mx3d,ctfary,ny,dfmin,dfmax,dfstep,apix,tytx)
	return jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfn_layered_ctf=jax.value_and_grad(prj_frc_layered_ctf)

def gradient_step_layered_ctf_optax(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,dsapix,weight=1.0,relstep=1.0,frc_Z=3.0):
	"""Computes one gradient step on the Gaussian coordinates given a set of particle FFTs at the appropriate scale,
	computing FRC to axial Nyquist, with specified linear weighting factor (def 1.0). Linear weight goes from
	0-2. 1 is unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns step, qual, shift, scale
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]
	mx=orts.to_mx3d()
	gausary=gaus.jax
	ptcls=ptclsfds.jax

	frcs,grad=gradvalfnl_layered_ctf(gausary,mx,ctfaryds,dfrange[0],dfrange[1],dfstep,dsapix,tytx,ptcls,weight, frc_Z)

	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std
	xyzs=relstep/(shift*500)   	# xyz scale factor, 1000 heuristic, TODO: may change

	return (grad,float(qual),float(shift),float(sca))

def prj_frc_layered_ctf_loss(gausary,mx3d,ctfary,dfmin,dfmax,dfstep,apix,tytx,ptcls,weight,frc_Z):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	prj=gauss_project_layered_ctf_fn(gausary,mx3d,ctfary,ny,dfmin,dfmax,dfstep,apix,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)

gradvalfnl_layered_ctf=jax.value_and_grad(prj_frc_layered_ctf_loss)



if __name__ == '__main__':
	main()

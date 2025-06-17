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

#@profile
def main():

	usage="""e3spa_refine_gauss.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--fromscratch", action="store_true",help="Ignore orientations from input file and refine from scratch")
	parser.add_argument("--volfiltlp", type=float, help="Lowpass filter to apply to output volume in A, 0 disables, default=40", default=40)
	parser.add_argument("--volfilthp", type=float, help="Highpass filter to apply to output volume in A, 0 disables, default=2500", default=2500)
	parser.add_argument("--apix", type=float, help="A/pix override for raw data", default=-1)
	parser.add_argument("--preclip", type=int, help="Trim the input images to the specified (square) box size in pixels", default=-1)
	parser.add_argument("--outbox",type=int,help="output boxsize, permitting over/undersampling (impacts A/pix)", default=-1)
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--spt", action="store_true",help="subtomogram averaging mode, changes optimization steps")
	parser.add_argument("--quick", action="store_true",help="single particle mode with less thorough refinement, but faster results")
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--ptcl3d_id", type=str, help="only use 2-D particles with matching ptcl3d_id parameter (lst file/header, use + for range)",default=None)
	parser.add_argument("--class", dest="classid", type=int, help="only use 2-D particles with matching class parameter (lst file/header)",default=-1)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--cachepath",type=str,help="path for storing the cached images, ideally on a high speed drive. Default='.'",default=".")
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
#	parser.add_argument("--precache",type=str,help="Rather than perform a reconstruction, only perform caching on the input file for later use. String is the folder to put the cache files in.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	jax_set_device(dev=0,maxmem=options.gpuram)
	llo=E3init(sys.argv,options.ppid)
	rng = np.random.default_rng()

	nptcl=EMUtil.get_image_count(args[0])
	if options.quick : nptcl=min(nptcl,8192+4096)

	if options.path == None:
		options.path=num_path_new("r3dgauss_")
		print("Writing to ", options.path)

	# Store the input arguments in the refinement directory
	db = js_open_dict(f"{options.path}/0_refine_parms.json")
	db.update(vars(options))
	db["commandline"]=" ".join(sys.argv)
	db["timestamp"]=str(time.ctime())
	db.close()
	lsxin=LSXFile(args[0])

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

	nxraw=EMData(args[0],0,True)["nx"]
	if options.preclip>0: nxraw=options.preclip
	nxrawm2=good_size_small(nxraw-2)
	if options.apix>0: apix=options.apix
	else: apix=EMData(args[0],0,True)["apix_x"]
	zmax=0.5
	if options.outbox>0: outsz=options.outbox
	else: outsz=min(1024,nxraw)

	if options.verbose: print(f"Input data box size {nxraw}x{nxraw} at {apix} A/pix. Maximum downsampled size for refinement {nxrawm2}. {nptcl} input images")

	# definition of downsampling sequence for stages of refinement
	# 0) #ptcl, 1) downsample, 2) iter, 3) frc weight, 4) amp threshold, 5) replicate, 6) repl. spread, 7) frc loc threshold (9 disables gradient)
	# replication skipped in final stage
	if options.fromscratch:
		stages=[
			[200,   24,16,1.8, 0  ,  1,.05,9],
			[200,   24,16,1.8, 0  ,  1,.05,9],
			[200,   24,16,1.8,-1  ,  2,.01,-1],
			[5000,  24,16,1.8,-1  ,  2,.1, -1],
			[5000,  24,16,1.8, 0  ,  8,.01,-1],
			[5000,  32,16,1.5,-.5 ,  8,.05,-3],
			[5000,  32,16,1.5,-1  , 16,.007,-2],
			[10000, 64,12,1.2,-1.5, 16,.005,-3],
			[10000, 64,12,1.0,-2  , 32,.002,-3],
			[10000,256,12,1.2,-1.5, 32,.005,-3],
			[10000,256,12,1.0,-2  ,128,.002,-3],
			[25000,512, 6,1.0,-2  ,128,.001,-3.0],
			[25000,512, 6,1.0,-2  ,  0,.001,-3.0]
		]
	elif options.spt:
		stages=[
			[2**10, 32, 32, 1.8, -1,  1,.05, 9], # 0
			[2**10, 32, 32, 1.8, -1,  1,.05, 9], # 1
			[2**10, 32, 48, 1.8, -1,  2,.05,-3], # 2
			[2**12, 64, 48, 1.5, -2,  2,.04, 9], # 3
			[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 4
			[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 5
			[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 6
			[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 7
			[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 8
			[2**12, 64,   8, 1.5, -2,  8,.04,-2], # 9
			[2**12, 64, 48, 1.5, -2,  8,.02, 9], # 10
			[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 11
			[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 12
			[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 13
			[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 14
			[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 15
			[2**12, 64,   8, 1.5, -2,16,.02,-2], # 16
			[2**14,128,32, 1.2, -3,16,.01, 9], # 17
			[2**14,128,  8, 1.2, -3,16,.01,-2], # 18
			[2**14,128,  8, 1.2, -3,16,.01,-2], # 19
			[2**14,128,  8, 1.2, -3,16,.01,-2], # 20
			[2**14,128,  8, 1.2, -3,32,.01,-2], # 21
			[2**16,256,16, 1.2, -2,32,.01, 9], # 22
			[2**16,256,  8, 1.2, -2,32,.01,-2], # 23
			[2**16,256,  8, 1.2, -2,32,.01,-2], # 24
			[2**16,256,  8, 1.2, -2,32,.01,-2], # 25
			[2**16,256,  8, 1.2, -2,32,.01,-2], # 26
			[2**17,512,  8, 1.2, -2,32,.01, 9], # 27
			[2**17,512,  4, 1.2, -2,32,.01,-2], # 28
			[2**17,512,  4, 1.2, -2,32,.01,-2], # 29
			[2**17,512,  4, 1.2, -2,32,.01,-2], # 30
			[2**17,512,  4, 1.2, -2,32,.01,-2], # 31
			[2**17,512,16, 1.2, -2,  0,.01, 9]  # 32
		]
	elif options.quick:
		stages=[
			[512,   16,24,1.8,-3  , 1,.01,9],
			[512,   16,24,1.8, 0  , 4,.01,9],
			[512,   16,24,1.8,-1  , 1,.01,-3],
			[512,   16,24,1.8, 0  , 4,.01,-2],
			[1024,  32,24,1.5, 0  , 8,.005,-3],
			[4096,  64,24,1.2,-1.5,16,.003,-2],
			[8192, 256,24,1.0,-2  , 0,.003,-3],
		]
	else:
		stages=[
			[512,    16,32,1.8,-3  , 1,.01,9],   # 0
			[512,    16,32,1.8, 0  , 1,.01,9],   # 1
			[512,    16,32,1.8,-3  , 4,.01,-3],  # 2
			[1024,   32,32,1.5, 0  , 4,.005,9],  # 3
			[1024,   32, 8,1.5,-1  , 4,.005,-2], # 4
			[1024,   32, 8,1.5,-1  , 4,.005,-2], # 5
			[1024,   32, 8,1.5,-1  , 4,.005,-2], # 6
			[1024,   32, 8,1.5,-1  , 8,.005,-2], # 7
			[1024,   32,32,1.5,-1  , 8,.005,9],  # 8
			[1024,   32, 8,1.5,-1  , 8,.005,-2], # 9
			[1024,   32, 8,1.5,-1  , 8,.005,-2], # 10
			[1024,   32, 8,1.5,-1  , 8,.005,-2], # 11
			[1024,   32, 8,1.5,-1  ,16,.005,-2], # 12
			[4096,   64,32,1.2,-1.5,16,.003,9],  # 13
			[4096,   64, 8,1.2,-1.5,16,.003,-3], # 14
			[4096,   64, 8,1.2,-1.5,16,.003,-3], # 15
			[4096,   64, 8,1.2,-1.5,16,.003,-3], # 16
			[4096,   64, 8,1.2,-1.5,32,.003,-3], # 17
			[16384, 128,32,1.0,-2  ,32,.003,9],  # 18
			[16384, 128, 8,1.0,-2  ,32,.003,-3], # 19
			[16384, 128, 8,1.0,-2  ,32,.003,-3], # 20
			[16384, 128, 8,1.0,-2  ,32,.003,-3], # 21
			[16384, 128, 8,1.0,-2  ,32,.003,-3], # 22
			[16384, 256,32,1.0,-2  ,32,.003,9],  # 23
			[16384, 256,8,1.0,-2  ,32,.003,-3],  # 24
			[16384, 256,8,1.0,-2  ,32,.003,-3],  # 25
			[16384, 256,8,1.0,-2  ,32,.003,-3],  # 26
			[16384, 256,8,1.0,-2  ,32,.003,-3],  # 27
			[16384, 512,16,1.0,-2  ,32,.003, 9], # 28
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 29
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 30
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 31
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 32
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 33
			[65536, 512, 8,0.8,-2  ,32,.001,-3], # 34
			[65536, 512,32,0.8,-2  , 0,.001,9]  # 35
		]

	# limit sampling to (at most) the box size of the raw data
	# we do this by altering stages to help with jit compilation
	for i in range(len(stages)):
		stages[i][1]=min(stages[i][1],nxrawm2)

	batchsize=192
	ctf_refine=18
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
		if orts is None and options.fromscratch==False:
			print(f"""\nERROR: No orientations found and --fromscratch not provided. Check that the input has
				xform.projections in the header or rerun with --fromscratch""")
			sys.exit(1)
		elif options.fromscratch:
			tytx=np.zeros((stk.shape[0],3))
			orts=rng.random((stk.shape[0],3))-0.5
		else: tytx/=jnp.array((nxraw,nxraw,1)) # Don't divide the defocus
		if options.ctf>0:
			mindf = min(mindf, float(jnp.min(tytx[:, 2])))
			maxdf = max(maxdf, float(jnp.max(tytx[:, 2])))
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
			orts,tytx=stk.orientations
			if orts is None and options.fromscratch==False:
				print(f"""\nERROR: No orientations found and --fromscratch not provided. Check that the input has
					xform.projections in the header or rerun with --fromscratch""")
				sys.exit(1)
			elif options.fromscratch: # TODO: the orts at least here will need to be put in an Orientations object instead of just a numpy array
				tytx=np.zeros((stk.shape[0],3))
				orts=rng.random((stk.shape[0],3))-0.5
			else: tytx/= jnp.array((nxraw,nxraw, 1))
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
		try:
			ctf=EMData(args[0],0,True)["ctf"] # Some .lst files have ctf info in header (the ones that have been phase flipped--if the input does not it currently will not get the defocuses corrct)
		except:
			print("Warning: Did not find CTF info in header of input. Checking for micrograph level CTF info")
			js=js_open_dict(info_name(EMData(args[0],0,True)["ptcl_source_image"])) # Assuming SPR uses lst file ptcls_XX.lst created by spt refinement which we have to go back to the per micrograph I think
			ctf=js["ctf"][0]
			js.close()
		jctf = jaxCTF(ctf=ctf, dfmajor=ctf.defocus, dfminor=ctf.defocus)
		dfstep = jctf.defocus_step
		boxlen = apix*stages[-1][1]*sqrt(3) # stages[-1][1] is the largest downsampling for the particle
		df_buffer = (boxlen/20000) + dfstep # TODO: This is a sensible buffer for multislice tomography but isn't working for SPA--and doesn't work for when we refine defocus
		dfrange=(mindf - df_buffer, maxdf + df_buffer)
		# Create the ctf stack
		ctf_stack,dfstep = jctf.compute_2d_stack_complex(nxraw, "amplitude", dfrange, "defocus")

	if options.fromscratch:
		# Reseed orientations for global search at low resolution
		tstorts=[]
		for x in np.arange(-0.5,0.5,0.04):
			for y in np.arange(-0.5,0.5,0.04):
				for z in np.arange(-0.5,0.5,0.04):
					if hypot(x,y,z)<=0.5: tstorts.append((x,y,z))
		tst_orts=Orientations(np.array(tstorts))

	if options.verbose>1: print(f"\n{local_datetime()}: Refining")

	gaus=Gaussians()
	#Initialize Gaussians to random values with amplitudes over a narrow range
	rnd=rng.uniform(0.0,1.0,(options.initgauss*9//10,4))		# start with completely random Gaussian parameters
#	rnd=tf.random.uniform((options.initgauss,4))     # specify the number of Gaussians to start with here
	neg = rng.uniform(0.0, 1.0, (options.initgauss//10, 4))		# 10% of gaussians are negative WRT the background (zero)
	rnd+=(-.5,-.5,-.5,2.0)
	neg+=(-.5,-.5,-.5,-3.0)
	rnd = np.concatenate((rnd, neg))
	gaus._data=rnd/(1.5,1.5,1.5,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size

	times.append(time.time())
	ptcls=[]
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")
		ccache = caches[stage[1]]

#		nliststg=range(sn,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity

		if options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tshift_grad\tamp_grad\timshift\tgrad_scale")
		# TODO: Ok, this should really use one of the proper optimization algorithms available from the deep learning toolkits
		# this basic conjugate gradient gets the job done, but not very efficiently I suspect...
		optim = optax.adam(.003)		# parm is learning rate
#		optim = optax.lion(.003)		# tried, seems not quite as good as Adam in test, but maybe worth another try
#		optim = optax.lamb(.005)		# tried, slightly better than adam, worse than lion
#		optim = optax.fromage(.01)		# tried, not as good
		optim_state=optim.init(gaus._data)		# initialize with data
		for i in range(stage[2]):		# training epochs
			if nptcl>stage[0]*2: idx0=sn+i
			else: idx0=0
			nliststg=range(idx0,nptcl,max(1,nptcl//stage[0]+1))		# all of the particles to use in the current epoch in the current stage, sn+i provides stochasticity
			imshift=0.0
			for j in range(0,len(nliststg)-10,batchsize):	# compute the gradient step piecewise due to memory limitations, batchsize particles at a time. The "-10" prevents very small residual batches from being computed
				ptclsfds,orts,tytx=ccache.read(nliststg[j:j+batchsize])
				if len(orts)<5 :
					print("Abort tiny batch: ",len(nliststg),j,batchsize)
					continue
				# standard mode, optimize gaussian parms only
				if options.ctf==0:
					step0,qual0,shift0,sca0=gradient_step_optax(gaus,ptclsfds,orts,tytx,stage[3])
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
					step0,qual0,shift0,sca0=gradient_step_layered_ctf_optax(gaus,ptclsfds,orts,jax_downsample_2d(ctf_stack, ptclsfds.shape[1]),tytx,dfrange,dfstep,dsapix,stage[3])
					step0=jnp.nan_to_num(step0)
					if j==0:
						step,qual,shift,sca=step0,-qual0,shift0,sca0
					else:
						step+=step0
						qual-=qual0
						shift+=shift0
						sca+=sca
				elif options.ctf==1:
					step0,qual0,shift0,sca0=gradient_step_ctf_optax(gaus,ptclsfds,orts,jax_downsample_2d(ctf_stack, ptclsfds.shape[1]),tytx,dfrange,dfstep,stage[3])
					step0=jnp.nan_to_num(step0)
					if j==0:
						step,qual,shift,sca=step0,-qual0,shift0,sca0
					else:
						step+=step0
						qual-=qual0
						shift+=shift0
						sca+=sca
			norm=len(nliststg)//batchsize+1
			qual/=norm
			shift/=norm
			sca/=norm
			imshift/=norm

			if isnan(shift) or isnan(sca) :
				if i==0:
					print("ERROR: nan on gradient descent, saving crash images and exiting")
					ptclsfds.do_ift().write_images("crash_lastb_images.hdf",0)
					out=open("crash_lastb_ortdydx.txt","w")
					for io in range(len(orts)):
						out.write(f"{orts[io][0]:1.6f}\t{orts[io][1]*1000:1.6f}\t{orts[io][2]*1000:1.6f}\t{tytx[io][0]*1000:1.2f}\t{tytx[io][1]*1000:1.2f} (/1000)\n")
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

			print(f"{i}\t{qual:1.5f}\t{shift*1000:1.6f}\t\t{sca*1000:1.6f}\t{imshift*1000:1.6f}  # /1000")

			if qual>0.99: break

		if stage[7]<9:
			if options.verbose: print(f"Adjusting translational alignment of particles")
			for j in range(0,nptcl,1000):	# compute the alignments piecewise due to memory limitations, 500 particles at a time
				ptclsfds,orts,tytx=ccache.read(range(j,min(j+1000,nptcl)))
				oldtytx=tytx
				tytx=ccf_step_align(gaus,ptclsfds,orts,tytx)
				ccache.tytx[range(j,min(j+1000,nptcl)),:2]=tytx
				dif=(tytx-oldtytx[:,:2])**2
				print(f"{j}-{j+1000}: shift rmsd: {sqrt(float(jnp.mean(dif)))*nxraw:.2f}")

			if options.fromscratch:
				# reseed orientations of particles with low FRCs
				# we do this by finding the best orientation with fixed angular sampling and a fixed box size of 24
				nseeded=0
				frcs=ccache.frcs			# not ideal, stealing the actual list from the object, but good enough for now
				lowfrc=frcs[frcs<1.5]
				if len(lowfrc)>0:
					frcm=np.mean(lowfrc)
					frcsg=np.std(lowfrc)
					reseed_idx=np.where(frcs<frcm+frcsg*stage[7])[0]			# [0] required because of odd np.where return
					nseeded=len(reseed_idx)

					if options.verbose: print(f"Making {len(tst_orts)} projections for reseeding")
					seedprojsf=gaus.project_simple(tst_orts,24).do_fft()		# fixed box size

					ptcls,tmpo,tmpxy=caches[24].read(reseed_idx)				# read all of the particle images we need to seed with new orientations, each is tiny with the fixed size of 24x24

					if options.verbose: print(f"Optimize {nseeded} orientations")
					for i in range(len(ptcls)):
						if options.verbose>1: print(f"{i}/{len(ptcls)}")
						ofrcs=jax_frc(seedprojsf.jax,ptcls[i],-1)
						maxort=jnp.argmax(ofrcs)			# best orientation for this particle
						ccache.orts[i]=tst_orts[maxort]
						#ccache.tytx[ii]=(0,0)			# just keep the current center?
					print(f"{nseeded} orts reseeded ({frcm+frcsg*stage[7]} thr)   {local_datetime()}")

			ort_optim = optax.adam(.003)		# parm is learning rate
			ort_optim_state=ort_optim.init((ccache.orts, ccache.tytx))		# initialize with data

			if options.verbose: print(f"\tIterating orientations parms x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tort_grad\tcen_grad")
			fout=open(f"{options.path}/fscs.txt","w")
			for i in range(stage[2]):		# training epochs
				if nptcl>stage[0]*2: idx0=sn+i
				else: idx0=0
				nliststg=range(idx0,nptcl,max(1,nptcl//stage[0]+1))		# all of the particles to use in the current stage, sn start gives some stochasticity
				#grads=jnp.zeros((ccache.orts.shape[0], ccache.	orts.shape[1]+ccache.tytx.shape[1]))
				ort_grads = jnp.zeros(ccache.orts.shape)
				tytx_grads = jnp.zeros(ccache.tytx.shape)

				norm=len(nliststg)//batchsize+1
				for j in range(0,len(nliststg),batchsize):	# compute the gradient step piecewise due to memory limitations, batchsize particles at a time
					ptclsfds,orts,tytx=ccache.read(nliststg[j:j+batchsize])
					if options.ctf > 0 and sn>ctf_refine:
						ort_step,tytx_step,qual0,ortstd0,dydxstd0=gradient_step_ort_ctf_optax(gaus,ptclsfds,orts,jax_downsample_2d(ctf_stack, ptclsfds.shape[1]),tytx,dfrange,dfstep,stage[3])
					else:
						ort_step,tytx_step,qual0,ortstd0,dydxstd0=gradient_step_ort_optax(gaus,ptclsfds,orts,tytx,stage[3])

					if j==0:
						qual,ortstd,dydxstd=-qual0,ortstd0,dydxstd0
					else:
						qual-=qual0
						ortstd+=ortstd0
						dydxstd+=dydxstd0

					ort_grads = ort_grads.at[jnp.array(nliststg[j:j+batchsize])].add(ort_step)
					tytx_grads = tytx_grads.at[jnp.array(nliststg[j:j+batchsize])].add(tytx_step)

					#print(len(nliststg[j:j+500]),ortstep.shape,dydxstep.shape)
#					ccache.add_orts(nliststg[j:j+batchsize],ortstep,dydxstep)
#					ccache.set_frcs(nliststg[j:j+batchsize],frcs)
#					for ii,n in enumerate(nliststg[j:j+batchsize]): fout.write(f"{n}\t{frcs[ii]:1.4f}\t{orts[ii][0]:1.4f}\t{orts[ii][1]:1.4f}\t{orts[ii][2]:1.4f}\n")

				ort_grads=jnp.nan_to_num(ort_grads)
				tytx_grads=jnp.nan_to_num(tytx_grads)
				qual/=norm
				ortstd/=norm
				dydxstd/=norm

				if isnan(ortstd) or isnan(dydxstd) :
					if i==0:
						print("ERROR: nan on gradient descent, saving crash images and exiting")
						ptclsfds.do_ift().write_images("crash_lastb_images.hdf",0)
						out=open("crash_lastb_ortdydx.txt","w")
						for io in range(len(orts)):
							out.write(f"{orts[io][0]:1.6f}\t{orts[io][1]*1000:1.6f}\t{orts[io][2]*1000:1.6f}\t{tytx[io][0]*1000:1.2f}\t{tytx[io][1]*1000:1.2f} (/1000)\n")
						sys.exit(1)
					else:
						print("ERROR: encountered nan on gradient descent, skipping epoch. Image numbers saved to crash_img_S_E.lst")
						try: os.unlinK("crash_img.lst")
						except: pass
						out=LSXFile(f"crash_img_{sn}_{i}.lst")
						for ii in nliststg: out.write(-1,ii,args[0])
						out=None
						continue

				ort_update, ort_optim_state = ort_optim.update((ort_grads, tytx_grads), ort_optim_state, (ccache.orts, ccache.tytx))
				(ccache.orts, ccache.tytx) = optax.apply_updates((ccache.orts, ccache.tytx), ort_update)
				ccache.orts = np.array(ccache.orts)
				ccache.tytx = np.array(ccache.tytx)

				print(f"{i}: {qual:1.4f}\t{ortstd:1.4f}\t\t{dydxstd:1.4f}")
		else: print("Skipping orientation gradient this step")

		# end of epoch, save images and projections for comparison
		if options.verbose>3:
			dsapix=apix*nxraw/ptclsfds.shape[1]
			mx2d=orts.to_mx2d(swapxy=True)
			gausary=gaus.jax
			ny=ptclsfds.shape[1]
			projs=EMStack2D(gauss_project_simple_fn(gausary,mx2d,ny,tytx))
			if options.ctf>0:
				mx3d=orts.to_mx3d()
				ctfaryds=jax_downsample_2d(ctf_stack, ny)
				ctf_projs=EMStack2D(gauss_project_ctf_fn(gausary,mx2d,ctfaryds,dfrange[0],dfstep,ny,tytx))
				layered_ctf_projs=EMStack2D(gauss_project_layered_ctf_fn(gausary,mx3d,ctfaryds,dfrange[0],dfstep,dsapix,ny,tytx))
				ctf = EMAN2Ctf()
				ctf.from_dict({"defocus":1.0, "voltage":jctf.voltage, "bfactor":0., "cs":jctf.cs, "ampcont":jctf.ampcont, "apix":jctf.apix})
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
					ctf.defocus = tytx[i,2]
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
					a["ctf"] = ctf
					b["ctf"] = ctf
					c["ctf"] = ctf
					d["ctf"] = ctf
					a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4)
					b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+1)
					c.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+2)
					d.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*4+3)
				else:
					a.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2)
					b.write_image(f"debug_img_{projs.shape[1]}.hdf:8",i*2+1)

		# filter results and prepare for stage 2
		if stage[5]>0:			# no filter/replicate in the last stage
			g0=len(gaus)
#			if options.tomo: gaus.norm_filter(sig=stage[4], cyl_mask=1/outsz)		# gaussians outside the box may be important!
			gaus.norm_filter(sig=stage[4],rad_downweight=0.33)
			g1=len(gaus)
			# Replicate gaussians to produce a specified total number for each stage. Critical that these numbers
			# fall in a small set of specific N for JIT compilation
			gaus.replicate_abs(stage[5]*options.initgauss,stage[6])
			g2=len(gaus)
		else: g0=g1=g2=len(gaus)
		print(f"{local_datetime()}: Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians")
		times.append(time.time())

		# do this at the end of each stage in case of early termination
		# Particle orientations
		lsxout=LSXFile(f"{options.path}/ptcls_{sn:02d}.lst")
		for i in range(len(lsxin)):
			a,b,c=lsxin[i]
			lsxout[i]=(a,b,{"xform.projection":Transform({"type":"spinvec","v1":float(ccache.orts[i][0]),"v2":float(ccache.orts[i][1]),"v3":float(ccache.orts[i][2]),"tx":float(ccache.tytx[i][1]*nxraw),"ty":float(ccache.tytx[i][0]*nxraw)}),"frc":float(ccache.frcs[i])})
		lsxout=None

		# Gaussian locations
		np.savetxt(f"{options.path}/threed_{sn:02d}.txt",gaus.numpy,fmt="%0.4f",delimiter="\t")

		# show individual shifts at high verbosity
#		if options.verbose>2:
#			print("TYTX: ",(caches[stage[1]].tytx*nxraw).astype(np.int32))

	times.append(time.time())
	vol=gaus.volume(outsz,zmax).center_clip(outsz)
	vol=vol.emdata[0]
	if options.sym not in ("c1","C1","I","i"):
		if options.verbose>0 : print(f"Apply {options.sym} symmetry to map (not gaussians)")
		vol.process_inplace("xform.applysym",{"sym":options.sym})
	times.append(time.time())
	vol["apix_x"]=apix*nxraw/outsz
	vol["apix_y"]=apix*nxraw/outsz
	vol["apix_z"]=apix*nxraw/outsz
	if options.ptcl3d_id is not None and options.ptcl3d_id>=0 : vol["ptcl3d_id"]=options.ptcl3d_id
	vol.write_image(f"{options.path}/threed_{sn:02d}_unfilt.hdf:12",-1)
	if options.volfilthp>0: vol.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.volfilthp})
	if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
	vol.process_inplace("normalize.edgemean")
	times.append(time.time())
	vol.write_image(f"{options.path}/threed_{sn:02d}.hdf:12",-1)

	times=np.array(times)
	#times-=times[0]
	times=times[1:]-times[:-1]
	if options.verbose>1 : print(times.astype(np.int32))

	E3end(llo)

# @profile
def gradient_step_optax(gaus,ptclsfds,orts,tytx,weight=1.0):
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

	frcs,grad=gradvalfnl(gausary,mx,tytx,ptcls,weight)

	qual=frcs			# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational (gauss) std
	sca=grad[:,3].std()			# amplitude std

	return (grad,float(qual),float(shift),float(sca))

# @profile
def prj_frc_loss(gausary,mx2d,tytx,ptcls,weight):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations. Returns -frc since optax wants to minimize, not maximize"""

	ny=ptcls.shape[1]
	#pfn=jax.jit(gauss_project_simple_fn,static_argnames=["boxsize"])
	#prj=pfn(gausary,mx2d,ny,tytx)
	prj=gauss_project_simple_fn(gausary,mx2d,ny,tytx)
#	print(prj.shape,ptcls.shape,weight,frc_Z)
#	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,2,frc_Z)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,1,3) # last arg is frc_z which we are trying to remove

gradvalfnl=jax.jit(jax.value_and_grad(prj_frc_loss))

def align_2d(gaus,orts,tytx,ptclsfds):
	ny=ptclsfds.shape[1]
	#ptcls=ptclsfds.jax
	mx=orts.to_mx2d(swapxy=True)
	prj=jax_fft2d(gauss_project_simple_fn(gaus.jax,mx,ny,tytx))	# FFT of gaussian projection for each particle
	return ptclsfds.align_translate(prj)/ny			# ccf between each particle and its projection

def gradient_step_ctf_optax(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,weight=1.0):
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

	frcs,grad=gradvalfnl_ctf(gausary,mx,ctfaryds,dfrange[0],dfstep,tytx,ptcls,weight)

	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std

	return (grad,float(qual),float(shift),float(sca))

def prj_frc_loss_ctf(gausary,mx2d,ctfary,dfmin,dfstep,tytx,ptcls,weight):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	prj=gauss_project_ctf_fn(gausary,mx2d,ctfary,dfmin,dfstep,ny,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,1,3)  # last arg is frc_z which we are trying to remove

gradvalfnl_ctf=jax.jit(jax.value_and_grad(prj_frc_loss_ctf), static_argnames=["dfmin","dfstep"])

def gradient_step_layered_ctf_optax(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,dsapix,weight=1.0):
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

	frcs,grad=gradvalfnl_layered_ctf(gausary,mx,ctfaryds,dfrange[0],dfstep,dsapix,tytx,ptcls,weight)

	qual=frcs					# functions used in jax gradient can't return a list, so frcs is a single value now
	shift=grad[:,:3].std()		# translational std
	sca=grad[:,3].std()			# amplitude std

	return (grad,float(qual),float(shift),float(sca))

def prj_frc_layered_ctf_loss(gausary,mx3d,ctfary,dfmin,dfstep,apix,tytx,ptcls,weight):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations."""

	ny=ptcls.shape[1]
	prj=gauss_project_layered_ctf_fn(gausary,mx3d,ctfary,dfmin,dfstep,apix,ny,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,1,3)  # last arg is frc_z which we are trying to remove

gradvalfnl_layered_ctf=jax.jit(jax.value_and_grad(prj_frc_layered_ctf_loss), static_argnames=["dfmin","dfstep","apix"])

def gradient_step_ort_optax(gaus,ptclsfds,orts,tytx,weight=1.0):
	"""Computes one gradient step on the orientation coordinates given a set of particle FFTs at the approprate scale,
	computing FRC to axial Nyquist, with the specified linear weighting factor (def 1.0). Linear weight goes from 0-2. 1 is
	unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns ort_step, dytx_step, qual, ort_std, tytx_std, and the 1-tensor with the per particle integrated FRCs for potential
	quality control.
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]
	gausary=gaus.jax
	ortary=orts.jax
	ptcls=ptclsfds.jax

	frcs, [gradort, gradtytx] = gradval_ol(gausary,ortary,tytx,ptcls,weight)

	qual=frcs
	stdort=gradort.std()		# orientation spinvec std
	stdtytx=gradtytx.std()	# tytx std

	return (gradort, gradtytx,float(qual),float(stdort),float(stdtytx))

def prj_frc_ort_loss(gausary,ortary, tytx, ptcls, weight):
	"""Aggregates the functions we need to take the gradient through. Computes the frc array resulting from the comparison
	of the Gaussians in gaus to particles in their current orientation"""
	ny=ptcls.shape[1]
	mx2d=jax_to_mx2d(ortary, swapxy=True)
	prj=gauss_project_simple_fn(gausary,mx2d,ny,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,1,3) # last arg is frc_z which we are trying to remove

gradval_ol=jax.jit(jax.value_and_grad(prj_frc_ort_loss, argnums=(1,2)))

def gradient_step_ort_ctf_optax(gaus,ptclsfds,orts,ctfaryds,tytx,dfrange,dfstep,weight=1.0):
	"""Computes one gradient step on the orientation coordinates given a set of particle FFTs at the approprate scale,
	computing FRC to axial Nyquist, with the specified linear weighting factor (def 1.0). Linear weight goes from 0-2. 1 is
	unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns ort_step, dytx_step, qual, ort_std, tytx_std, and the 1-tensor with the per particle integrated FRCs for potential
	quality control.
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]
	gausary=gaus.jax
	ortary=orts.jax
	ptcls=ptclsfds.jax

	frcs, [gradort, gradtytx] = gradval_olc(gausary,ortary,ctfaryds,dfrange[0],dfstep,tytx,ptcls,weight)
#	print("gradtytx", gradtytx)

	qual=frcs
	stdort=gradort.std()		# orientation spinvec std
	stdtytx=gradtytx.std()	# tytx std

	return (gradort, gradtytx,float(qual),float(stdort),float(stdtytx))

def prj_frc_ort_ctf_loss(gausary,ortary,ctfaryds,dfmin,dfstep,tytx, ptcls, weight):
	"""Aggregates the functions we need to take the gradient through. Computes the frc array resulting from the comparison
	of the Gaussians in gaus to particles in their current orientation"""
	ny=ptcls.shape[1]
	mx2d=jax_to_mx2d(ortary, swapxy=True)
	prj=gauss_project_ctf_fn(gausary,mx2d,ctfaryds,dfmin,dfstep,ny,tytx)
	return -jax_frc_jit(jax_fft2d(prj),ptcls,weight,1,3) # last arg is frc_z which we are trying to remove

gradval_olc=jax.jit(jax.value_and_grad(prj_frc_ort_ctf_loss, argnums=(1,5)), static_argnames=["dfmin","dfstep"])


def ccf_step_align(gaus,ptclsfds,orts,tytx):
	"""Uses CCF to update all translational alignments in one step with CCF"""
	ny=ptclsfds.shape[1]
	mx = orts.to_mx2d(swapxy=True)
	# we are determining absolute shifts, so we get rid of the original shift
	projsf=EMStack2D(gauss_project_simple_fn(gaus.jax,mx,ny,jnp.zeros(tytx.shape))).do_fft()
	newtytx= ptclsfds.align_translate(projsf).astype(jnp.float32)/float(-ny)

	# if ny>64:
	# 	out=open("dbug_xy.txt","w")
	# 	out.write("# dx;dy\n")
	# 	for i in range(newtytx.shape[0]): out.write(f"{newtytx[i][1]}\t{newtytx[i][0]}\n")
 #
	# 	projsf=gaus.project_simple(orts,ny,newtytx).do_fft()
	# 	newtytx=tf.cast(ptclsfds.align_translate(projsf),tf.float32)/float(-ny)
 #
	# 	out=open("dbug_xy2.txt","w")
	# 	out.write("# dx;dy\n")
	# 	for i in range(newtytx.shape[0]): out.write(f"{newtytx[i][1]}\t{newtytx[i][0]}\n")
	# 	# ccf=projsf.calc_ccf(ptclsfds).do_ift()
	# 	# ccf.write_images("dbug_ccf.hdf",0)
	# 	# projs=projsf.do_ift()
	# 	# ptclsds=ptclsfds.do_ift()
	# 	# projs.write_images("dbug_prj.hdf")
	# 	# ptclsds.write_images("dbug_ptcl.hdf")
	# 	sys.exit(1)

	return newtytx


if __name__ == '__main__':
	main()

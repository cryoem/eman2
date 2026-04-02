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
from collections import defaultdict

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

	usage="""e3spa_refine_point.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--fromscratch", action="store_true",help="Ignore orientations from input file and refine from scratch")
	parser.add_argument("--volfiltlp", type=float, help="Lowpass filter to apply to output volume in A, 0 disables, default=5", default=5)
	parser.add_argument("--frc_weight", type=float, help="Testing only at present", default=-1)
	parser.add_argument("--apix", type=float, help="A/pix override for raw data", default=-1)
	# parser.add_argument("--preclip", type=int, help="Trim the input images to the specified (square) box size in pixels", default=-1) # Not supported by new caching mechanism
	parser.add_argument("--outbox",type=int,help="output boxsize, permitting over/undersampling (impacts A/pix)", default=-1)
	parser.add_argument("--initpoint",type=int,help="Points in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--keep", type=float, help="The fraction of images to use, based on quality scores (1.0 = use all)",default=1.0)
	parser.add_argument("--ptcl3d_id", type=str, help="only use 2-D particles with matching ptcl3d_id parameter (lst file/header, use + for range)",default=None)
	parser.add_argument("--class", dest="classid", type=int, help="only use 2-D particles with matching class parameter (lst file/header)",default=-1)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	# parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	# parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--fscdebug", type=str,help="Compute the FSC of the final map with a reference volume for debugging",default=None)
	parser.add_argument("--score",type=str,help="If set, will generate the specified .lst file with score set and plottable scores.txt",default=None)
	parser.add_argument("--cachepath",type=str,help="path for storing the cached images, ideally on a high speed drive. Default='.'",default=".")
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")


	(options, args) = parser.parse_args()
	# jax_set_device(dev=0,maxmem=options.gpuram)
	llo=E3init(sys.argv,options.ppid)
	rng=np.random.default_rng()

	os.putenv("EMAN3_CACHE_PATH",options.cachepath)
	# options.keep=[float(k) for k in options.keep.split(',')]
	# if len(options.keep)<3: options.keep=[options.keep[0]]*3
	# TODO: Make sure works since forced keep to be single number

	nptcl=EMUtil.get_image_count(args[0]) # Initial number of particles
	if options.path == None:
		options.path=num_path_new("r3dpoint_")
		print("Writing to ", options.path)

	if args[0][-4:]!=".lst" : error_exit("Only LST input files are supported. If working with an HDF with embedded parameters, create a .lst file and extract the parameters with e2proclst.py")
	# Store the input arguments in the refinement directory
	db = js_open_dict(f"{options.path}/0_refine_parms.json")
	db.update(vars(options))
	db["commandline"]=" ".join(sys.argv)
	db["timestamp"]=str(time.ctime())
	db.close()
	lsxin=LSXFile(args[0])

	# Particle selection based on various options
	if (options.keep==1.0) and (options.ptcl3d_id is None): selimg=set(range(nptcl))
	else:
		p2ds=[]						# construct list of (score,#) for 2-D particles at the same time
		for i,l in enumerate(lsxin):
			p2ds.append((l[2]["score"],i))
		p2ds.sort()
		selimg=set([i for s,i in p2ds[:int(len(p2ds)*options.keep)]])

	if options.classid>=0:
		selcls=set([i for i in range(len(lsxin)) if lsxin[i][2]["class"]==options.classid])
		selimg=selimg.intersection(selcls)

	selimg=np.array(list(selimg))   # need to convert set to list before going to array or we wind up with an array with a set in it
	selimg.sort()
	nptcl=len(selimg) # Update number of particles to only those selected
	if options.verbose>0: print(f"{nptcl}/{len(lsxin)} 2-D images selected")


	nxraw=EMData(args[0],0,True)["nx"]
#	if options.preclip>0: nxraw=options.preclip
	if options.apix>0: apix=options.apix
	else: apix=EMData(args[0],0,True)["apix_x"]
	zmax=0.5
	if options.outbox>0: outsz=options.outbox
	else: outsz=min(1024,nxraw)

	if options.verbose: print(f"Input data box size {nxraw}x{nxraw} at {apix} A/pix. Maximum downsampled size for refinement {nxraw}. Thickness limit +-{zmax}. {nptcl} input images")

	if options.fscdebug is not None:
		dbugvol=EMStack3D(EMData(options.fscdebug)).do_fft().jax
		rad_vol_int(dbugvol.shape[2])
	else: dbugvol=None

	if options.verbose: print(f"{nptcl} particles at {nxraw}^3")

	# definition of downsampling sequence for stages of refinement
	# 0) #ptcl, 1) downsample, 2) iter, 3) frc weight, 4) amp threshold, 5) replicate, 6) repl. spread, 7) frc loc threshold for ort refinement, (0 indicates Point stage, negative means ort refinement)
	# replication skipped in final stage
	if options.fromscratch:
		stages=[
			[200,   32,16,1.8, 0  ,  1,.05,9],
			[200,   32,16,1.8, 0  ,  1,.05,9],
			[200,   32,16,1.8,-1  ,  2,.01,-1],
			[5000,  32,16,1.8,-1  ,  2,.1, -1],
			[5000,  32,16,1.8, 0  ,  8,.01,-1],
			[5000,  32,16,1.5,-.5 ,  8,.05,-3],
			[5000,  32,16,1.5,-1  , 16,.007,-2],
			[10000, 64,12,1.2,-1.5, 16,.005,-3],
			[10000, 64,12,1.0,-2  , 32,.002,-3],
			[10000,256,12,1.2,-1.5, 32,.005,-3],
			[10000,256,12,1.0,-2  ,128,.002,-3],
			[25000,512, 6,1.0,-2  ,128,.001,-3.0],
			[25000,512, 6,1.0,-2  ,  0,.001,-3.0]
		]
	# elif options.spt: # TODO: Leaving for now until move to spt refine program, but will be removed
	# 	stages=[
	# 		[2**10, 32, 32, 1.8, -1,  1,.05, 9], # 0
	# 		[2**10, 32, 32, 1.8, -1,  1,.05, 9], # 1
	# 		[2**10, 32, 48, 1.8, -1,  2,.05,-3], # 2
	# 		[2**12, 64, 48, 1.5, -2,  2,.04, 9], # 3
	# 		[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 4
	# 		[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 5
	# 		[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 6
	# 		[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 7
	# 		[2**12, 64,   8, 1.5, -2,  2,.04,-2], # 8
	# 		[2**12, 64,   8, 1.5, -2,  8,.04,-2], # 9
	# 		[2**12, 64, 48, 1.5, -2,  8,.02, 9], # 10
	# 		[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 11
	# 		[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 12
	# 		[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 13
	# 		[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 14
	# 		[2**12, 64,   8, 1.5, -2,  8,.02,-2], # 15
	# 		[2**12, 64,   8, 1.5, -2,16,.02,-2], # 16
	# 		[2**14,128,32, 1.2, -3,16,.01, 9], # 17
	# 		[2**14,128,  8, 1.2, -3,16,.01,-2], # 18
	# 		[2**14,128,  8, 1.2, -3,16,.01,-2], # 19
	# 		[2**14,128,  8, 1.2, -3,16,.01,-2], # 20
	# 		[2**14,128,  8, 1.2, -3,32,.01,-2], # 21
	# 		[2**16,256,16, 1.2, -2,32,.01, 9], # 22
	# 		[2**16,256,  8, 1.2, -2,32,.01,-2], # 23
	# 		[2**16,256,  8, 1.2, -2,32,.01,-2], # 24
	# 		[2**16,256,  8, 1.2, -2,32,.01,-2], # 25
	# 		[2**16,256,  8, 1.2, -2,32,.01,-2], # 26
	# 		[2**17,512,  8, 1.2, -2,32,.01, 9], # 27
		# 	[2**17,512,  4, 1.2, -2,32,.01,-2], # 28
		# 	[2**17,512,  4, 1.2, -2,32,.01,-2], # 29
		# 	[2**17,512,  4, 1.2, -2,32,.01,-2], # 30
		# 	[2**17,512,  4, 1.2, -2,32,.01,-2], # 31
		# 	[2**17,512,16, 1.2, -2,  0,.01, 9]  # 32
		# ]
	else:
		stages=[
			[512,   16,32,1.8,-3  ,1,.01, 0], # 0: Point
			[512,   16,32,1.8, 0  ,4,.01, 0], # 1: Point
			[1024,  32,32,1.5, 0  ,4,.005,0], # 2: Point
			[1024,  32,32,1.5,-1  ,8,.005,0], # 3: Point
			[4096,  64,32,1.2,-1.5,16,.003,0], # 4: Point
			[16384, 256,32,1.0,-2 ,32,.003,0], # 5: Point
			[65536*2, 512,32,1.0,-2 ,32,.001,0], # 6: Point
			[65536*2, 512, 16,0.8,-2  ,0,.001,0], # 7: Points: don't filter
			[65536*2, 512, 16,0.8,-2  ,32,.001,-3], # 8: Orientations
			[65536*2, 512, 16,0.8,-2  ,0,.001,0], # 9: Points: don't filter
			[65536*2, 512, 16,0.8,-2  ,32,.001,-3], # 10: Orientations
			[65536*2, 512, 16,0.8,-2  ,0,.001,0], # 11: Points: don't filter
			[65536*2, 512, 16,0.8,-2  ,32,.001,-3], # 12: Orientations
			[65536*2, 512, 16,0.8,-2  ,0,.001,0], # 13: Points: don't filter
			[65536*2, 512, 16,0.8,-2  ,32,.001,-3], # 14: Orientations
			[65536*2, 512,16,1.0,-2 ,0,.001,0], # 15: Points
		]

	# limit sampling to (at most) the box size of the raw data
	# we do this by altering stages to help with jit compilation
	for i in range(len(stages)):
		stages[i][1]=min(stages[i][1],nxraw)
		if options.frc_weight>0: stages[i][3]=options.frc_weight

	# Setting up symmetry
	sym=parsesym(options.sym)
	sym_orts = Orientations()
	sym_orts.init_from_transforms(sym.get_syms())
	symmx = sym_orts.to_mx3d()

	# prior to jax 0.7.x there was a sharding problem in JAX causing a crash related to
	# each image in the batch getting turned into a separate argument when JIT compiling
	# in 0.7.x, larger batches can be used and have some advantages
	if int(jax.__version__.split(".")[1])>6 :
		batchsize=1024//len(sym_orts)
	else:
		batchsize=192

	times=[time.time()]
	# Cache initialization
	if options.verbose: print(f"{local_datetime()}: Caching particle data")
	downs=sorted(set([s[1] for s in stages]))

	# critical for later in the program, this initializes the radius images for all of the samplings we will use
	for d in downs:
		rad_img_int(d)
		if options.ctf>0:
			rad2_img(d)

	### Note that StackCache stores the entire .lst file. selimg must be handled when reading from the cache
	cache=StackCache(args[0])

	if options.fromscratch:
		# Reseed orientations for global search at low resolution
		tstorts=[]
		for x in np.arange(-0.5,0.5,0.04):
			for y in np.arange(-0.5,0.5,0.04):
				for z in np.arange(-0.5,0.5,0.04):
					if hypot(x,y,z)<=0.5: tstorts.append((x,y,z))
		tst_orts=Orientations(np.array(tstorts))

	if options.verbose>1: print(f"\n{local_datetime()}: Refining")

	# Initialize Points
	point=Points()
	#Initialize Points to random values with amplitudes over a narrow range
	rng = np.random.default_rng()
	rnd=rng.uniform(0.0,1.0,(options.initpoint*9//10,4))		# start with completely random Point parameters
	neg = rng.uniform(0.0, 1.0, (options.initpoint//10, 4))		# 10% of points are negative WRT the background (zero)
	rnd+=(-.5,-.5,-.5,2.0)
	neg+=(-.5,-.5,-.5,-3.0)
	rnd = np.concatenate((rnd, neg))
	point._data=rnd/(1.5,1.5,1.5,3.0)	# amplitudes set to ~1.0, positions random within 2/3 box size

	frchist=[]
	times.append(time.time())
	ptcls=[]
	qualities=np.zeros((nptcl,len(stages)))			# only used with --score
	weights=[None for i in range(len(stages))]		# saved, but not used at present
	thresholds=[None for i in range(len(stages))]	# saved, but not used at present
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")

		if stage[7]==0: # Refining Points
			if  options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]} to refine Points\n    FRC\t\tshift_grad\tamp_grad\timshift\tgrad_scale")
			optim = optax.adam(.003)		# parm is learning rate
			optim_state=optim.init(point._data)		# initialize with data
			for i in range(stage[2]):		# training epochs
				if nptcl>stage[0]*2: idx0=sn+i
				else: idx0=0
				nliststg=selimg[range(idx0,nptcl,max(1,nptcl//stage[0]+1))]		# all of the particles to use in the current epoch in the current stage, sn+i provides stochasticity
				imshift=0.0
				for j in range(0,len(nliststg)-10,batchsize):	# compute the gradient step piecewise due to memory limitations, batchsize particles at a time. The "-10" prevents very small residual batches from being computed
					ptclsfds=cache.read(stage[1],nliststg[j:j+batchsize])	# metadata stored in ptclsfds, which is a EMStack2D
				# Since we aren't modifying the metadata in this program, just make a JAX copy at the get-go
					meta=jnp.array(ptclsfds.metadata)		# 0:ty,1:tx,2:ortx,3:orty,4:ortz,5:defocus,6:phase,7:dfdiff,8:astigangle,9:score,10:class

					if len(ptclsfds)<5 :
						print("Abort tiny batch: ",len(nliststg),j,batchsize)
						continue

					# on the first epoch of each stage we look at the variance of the fsc curve to estimate weighting
					# only using a single batch for this right now. May be sufficient
					# TODO - Do we need to do this differently with other CTF modes? Variance estimate might not be impacted by the phase flipping, so may be ok...
					if i in (0,8) and j==0:
						frcs=prj_frcs(point.jax,ptclsfds,meta)
						try:
							thresh=1.25*np.std(frcs,0)/sqrt(batchsize)
							weight=1.0/np.array(thresh)		# this should make all of the standard deviations the same
							weight[0:2]=0			# low frequency cutoff
							weight[ptclsfds.shape[1]//2:]=0
							weight/=np.sum(weight)	# normalize to 1
							weight=jnp.array(weight*len(weight))	# the *len(weight) is dumb, but due to mean() being returned
						except:
							print(f"Weighting failed {sn},{i},{j}")
							weight=np.ones((len(frcs.shape[1])))
						thresholds[sn]=thresh
						weights[sn]=weight
						frchist.append((np.array(np.mean(frcs,0)),thresh,weight))


					if options.ctf==0:
						step0,qual0,shift0,sca0=point_gradient_step_optax(point,ptclsfds,meta,symmx,weight,thresh)
						# TODO: These nan_to_num shouldn't be necessary. Not sure what is causing nans. Could be there is an implicit sqrt somewhere we're taking the gradient of, or a roundoff error?
						step0=jnp.nan_to_num(step0)
						qual0=jnp.nan_to_num(qual0)
						shift0=jnp.nan_to_num(shift0)
						sca0=jnp.nan_to_num(sca0)
						if j==0:
							step,qual,shift,sca=step0,-qual0,shift0,sca0
						else:
							step+=step0
							qual-=qual0
							shift+=shift0
							sca+=sca0
					elif options.ctf==1:
						dsapix=ptclsfds.apix
						wavelength=12.2639/np.sqrt(ptclsfds.voltage*1000.0+0.97845*ptclsfds.voltage*ptclsfds.voltage)
						dfstep=2*apix*apix/(wavelength*10000)
						step0,qual0,shift0,sca0=point_gradient_step_ctf_optax(point,ptclsfds,meta,jnp.array([wavelength, ptclsfds.cs]),dfstep,dsapix,symmx,weight,thresh)
						step0=jnp.nan_to_num(step0)
						if j==0:
							step,qual,shift,sca=step0,-qual0,shift0,sca0
						else:
							step+=step0
							qual-=qual0
							shift+=shift0
							sca+=sca0
					elif options.ctf==2:
						dsapix=ptclsfds.apix
						wavelength=12.2639/np.sqrt(ptclsfds.voltage*1000.0+0.97845*ptclsfds.voltage*ptclsfds.voltage)
						dfstep=2*apix*apix/(wavelength*10000)
						step0,qual0,shift0,sca0=point_gradient_step_layered_ctf_optax(point,ptclsfds,meta,jnp.array([wavelength, ptclsfds.cs]),dfstep,dsapix,symmx,weight,thresh)
						step0=jnp.nan_to_num(step0)
						if j==0:
							step,qual,shift,sca=step0,-qual0,shift0,sca0
						else:
							step+=step0
							qual-=qual0
							shift+=shift0
							sca+=sca0

				norm=len(nliststg)//batchsize+1
				if norm==0: raise Exception("ERROR: norm zero. This shouldn't happen")
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

				update, optim_state = optim.update(step, optim_state, point._data)
				point._data = optax.apply_updates(point._data, update)

				if dbugvol is not None:
					nyd=dbugvol.shape[1]
					if options.sym not in ("c1","C1","I","i"):
						vol=point.volume(nyd,zmax)
						vol.emdata[0].process_inplace("xform.applysym",{"sym":options.sym})
						vol=vol.do_fft().jax
					else: vol=point.volume(nyd,zmax).do_fft().jax
					fsc=jax_fsc_jit(vol,dbugvol)
					out=open(f"fscm3d_{sn:02d}_{i:02d}.txt","w")
					for s in range(nyd//2): out.write(f"{s/nyd:1.4f}\t{float(fsc[0][s]):1.5f}\n")
					out=None
					dbfsc=float(fsc[0][:nyd//2].mean())
					print(f"{dbfsc:0.5f}\t",end="")

				print(f"{i}\t{qual:1.5f}\t{shift*1000:1.6f}\t\t{sca*1000:1.6f}\t{imshift*1000:1.6f}  # /1000")

				if qual>0.99: break

				# Particle quality assessment
				if options.score is not None and i==stage[2]-1:
					if options.verbose: print("Compute quality")
					for j in range(0,nptcl,batchsize):
						ptclsfds=cache.read(stage[1],nliststg[j:j+batchsize])	# metadata stored in ptclsfds, which is a EMStack2D
						meta=jnp.array(ptclsfds.metadata)		# 0:ty,1:tx,2:ortx,3:orty,4:ortz,5:defocus,6:phase,7:dfdiff,8:astigangle,9:score,10:class

						# need to recompute this here since we may not have hit this stage yet
						if j==0:
							frcs=prj_frcs(point.jax,ptclsfds,meta)
							try:
								thresh=1.25*np.std(frcs,0)/sqrt(batchsize)
								weight=1.0/np.array(thresh)		# this should make all of the standard deviations the same
								weight[0:2]=0			# low frequency cutoff
								weight[ptclsfds.shape[1]//2:]=0
								weight/=np.sum(weight)	# normalize to 1
								weight=jnp.array(weight*len(weight))	# the *len(weight) is dumb, but due to mean() being returned
							except:
								print(f"Weighting failed {sn},{i},{j}")
								weight=np.ones((len(frcs.shape[1])))

						# we measure per particle quality and save it
						quality=jnp.nan_to_num(sym_prj_frcs_loss(point.jax,symmx,ptclsfds.jax,meta,weight,thresh), nan=2.0)

						#print(quality.shape,quality)
						for ii,q in enumerate(np.array(quality)): qualities[ii+j][sn]=q

		if stage[7]<0:
			if options.verbose: print(f"Adjusting translational alignment of particles")
			for j in range(0,nptcl,1000):	# compute the alignments piecewise due to memory limitations, 1000 particles at a time
				ptclsfds=cache.read(stage[1],selimg[range(j,min(j+1000,nptcl))])
				orts, tytx=ptclsfds.orientations
				tytx=jnp.array(tytx)
				oldtytx=tytx
				tytx=ccf_step_align(point,ptclsfds,orts,tytx)
				ptclsfds.metadata[:,:2]=np.array(tytx)
				dif=(tytx-oldtytx)**2
				print(f"{j}-{j+1000}: shift rmsd: {sqrt(float(jnp.mean(dif)))*nxraw:.2f}")

##########################Needs altering to use meta--hence commented out
			# if options.fromscratch:
			# 	# reseed orientations of particles with low FRCs
			# 	# we do this by finding the best orientation with fixed angular sampling and a fixed box size of 24
			# 	nseeded=0
			# 	frcs=ptclsfds			# not ideal, stealing the actual list from the object, but good enough for now
			# 	lowfrc=frcs[frcs<1.5]
			# 	if len(lowfrc)>0:
			# 		frcm=np.mean(lowfrc)
			# 		frcsg=np.std(lowfrc)
			# 		reseed_idx=np.where(frcs<frcm+frcsg*stage[7])[0]			# [0] required because of odd np.where return
			# 		nseeded=len(reseed_idx)
   #
			# 		if options.verbose: print(f"Making {len(tst_orts)} projections for reseeding")
			# 		seedprojsf=point.project_simple(tst_orts,24).do_fft()		# fixed box size
   #
			# 		ptcls,tmpo,tmpxy,tmpa=caches[24].read(reseed_idx)				# read all of the particle images we need to seed with new orientations, each is tiny with the fixed size of 24x24
   #
			# 		if options.verbose: print(f"Optimize {nseeded} orientations")
			# 		for i in range(len(ptcls)):
			# 			if options.verbose>1: print(f"{i}/{len(ptcls)}")
			# 			ofrcs=jax_frc(seedprojsf.jax,ptcls[i],-1)
			# 			maxort=jnp.argmax(ofrcs)			# best orientation for this particle
			# 			ccache.orts[i]=tst_orts[maxort]
			# 			#ccache.tytx[ii]=(0,0)			# just keep the current center?
			# 		print(f"{nseeded} orts reseeded ({frcm+frcsg*stage[7]} thr)   {local_datetime()}")
			# Make jax copies of tytx and orts for gradient descent
			tytx=jnp.array(cache._meta[:,:2])
			orts=jnp.array(cache._meta[:,2:5])

			ort_optim = optax.adam(.001)
			ort_optim_state=ort_optim.init((orts, tytx))		# initialize with data

			if options.verbose: print(f"\tIterating orientations parms x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tort_grad\tcen_grad")
			fout=open(f"{options.path}/fscs.txt","w")
			for i in range(stage[2]):		# training epochs
				ort_grads = jnp.zeros(orts.shape)
				tytx_grads = jnp.zeros(tytx.shape)

				norm=nptcl//batchsize+1
				for j in range(0,nptcl,batchsize):
					ptclsfds=cache.read(stage[1],selimg[range(j,min(j+batchsize,nptcl))])
					meta=jnp.array(ptclsfds.metadata)# 0:ty,1:tx,2:ortx,3:orty,4:ortz,5:defocus,6:phase,7:dfdiff,8:astigangle,9:score,10:class
					dsapix = ptclsfds.apix

					# TODO: Same question--do we need to do this differently with other CTF modes?
					if i in (0,8) and j==0:
						frcs=prj_frcs(point.jax,ptclsfds,meta)
						try:
							thresh=1.25*np.std(frcs,0)/sqrt(batchsize)
							weight=1.0/np.array(thresh)		# this should make all of the standard deviations the same
							weight[0:2]=0			# low frequency cutoff
							weight[ptclsfds.shape[1]//2:]=0
							weight/=np.sum(weight)	# normalize to 1
							weight=jnp.array(weight*len(weight))	# the *len(weight) is dumb, but due to mean() being returned
						except:
							print(f"Weighting failed {sn},{i},{j}")
							weight=np.ones((len(frcs.shape[1])))
						frchist.append((np.array(np.mean(frcs,0)),thresh,weight))

					if options.ctf ==0:
						ort_step,tytx_step,qual0,ortstd0,dydxstd0=ort_gradient_step_optax(point,ptclsfds,meta,symmx,weight,thresh)
					elif options.ctf ==1:
						wavelength=12.2639/np.sqrt(ptclsfds.voltage*1000.0+0.97845*ptclsfds.voltage*ptclsfds.voltage)
						dfstep=2*apix*apix/(wavelength*10000)
						ort_step,tytx_step,qual0,ortstd0,dydxstd0=ort_gradient_step_ctf_optax(point,ptclsfds,meta,jnp.array([wavelength,ptclsfds.cs]),dfstep,dsapix,symmx,weight,thresh)
					elif options.ctf == 2:
						wavelength=12.2639/np.sqrt(ptclsfds.voltage*1000.0+0.97845*ptclsfds.voltage*ptclsfds.voltage)
						dfstep=2*apix*apix/(wavelength*10000)
						ort_step,tytx_step,qual0,ortstd0,dydxstd0=ort_gradient_step_layered_ctf_optax(point,ptclsfds,meta,jnp.array([wavelength,ptclsfds.cs]),dfstep,dsapix,symmx,weight,thresh)

					if j==0:
						qual,ortstd,dydxstd=-qual0,ortstd0,dydxstd0
					else:
						qual-=qual0
						ortstd+=ortstd0
						dydxstd+=dydxstd0

					ort_grads = ort_grads.at[jnp.arange(j, min(j+batchsize, nptcl))].add(ort_step)
					tytx_grads = tytx_grads.at[jnp.arange(j, min(j+batchsize, nptcl))].add(tytx_step)

				# TODO: again, nan_to_num shouldn't be necessary. What is causing it?'
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

				ort_update, ort_optim_state = ort_optim.update((ort_grads, tytx_grads), ort_optim_state, (orts, tytx))
				(orts, tytx) = optax.apply_updates((orts, tytx), ort_update)

				print(f"{i}: {qual:1.4f}\t{ortstd:1.4f}\t\t{dydxstd:1.4f}")

			# Save the changes we've made to the np array so it goes to all levels of downsampling
			cache._meta[:,:2]=np.array(tytx)
			cache._meta[:,2:5]=np.array(orts)

		# end of epoch, save images and projections for comparison
		if options.verbose>3:
			dsapix=ptclsfds.apix
			pointary=point.jax
			ny=ptclsfds.shape[1]
			orts, tytx=ptclsfds.orientations
			projs=EMStack2D(point_project_simple_sym_fn(pointary, orts.jax, ny, tytx, symmx))
			if options.ctf>0:
				dsapix=ptclsfds.apix
				ctf=ptclsfds.ctf
				wavelength=12.2639/np.sqrt(ptclsfds.voltage*1000.0+0.97845*ptclsfds.voltage*ptclsfds.voltage)
				dfstep=2*apix*apix/(wavelength*10000)
				ctf_projs=EMStack2D(point_project_ctf_sym_fn(pointary, orts.jax, jnp.array([wavelength,ptclsfds.cs]), dfstep, dsapix, ny, tytx, ctf, symmx))
				layered_ctf_projs=EMStack2D(point_project_layered_ctf_sym_fn(pointary,orts.jax,jnp.array([wavelength,ptclsfds.cs]),dfstep,dsapix,ny,tytx,ctf, symmx))
			ptclds=ptclsfds.do_ift()
			transforms=orts.transforms(tytx=tytx)
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

		# filter results and prepare for stage 2
		if stage[5]>0:			# no filter/replicate in the last stage
			g0=len(point)
			point.norm_filter(sig=stage[4],rad_downweight=0.33)
			g1=len(point)
			# Replicate points to produce a specified total number for each stage. Critical that these numbers
			# fall in a small set of specific N for JIT compilation
			point.replicate_abs(stage[5]*options.initpoint,stage[6])
			g2=len(point)
		else: g0=g1=g2=len(point)
		print(f"{local_datetime()}: Stage {sn} complete: {g0} -> {g1} -> {g2} points")
		times.append(time.time())

		# do this at the end of each stage in case of early termination
		# Particle orientations
		lsxout=LSXFile(f"{options.path}/ptcls_{sn:02d}.lst")
		if lsxin.parent is not None: lsxout.parent=lsxin.parent
		else: lsxout.parent=args[0]
		for i in range(len(lsxin)):
			a,b,c=lsxin[i]
			ctf=EMAN2Ctf()
			ctf.from_dict({"defocus":float(cache._meta[i][5]), "dfdiff":float(cache._meta[i][7]), "dfang":float(cache._meta[i][8]), "apix":float(apix), "voltage":int(ptclsfds.voltage), "cs":float(ptclsfds.cs)})
			ctf.set_phase(float(cache._meta[i][6]))
			lsxout[i]=(a,b,{"xform.projection":Transform({"type":"spinvec","v1":float(cache._meta[i][2]),"v2":float(cache._meta[i][3]),"v3":float(cache._meta[i][4]),"tx":float(cache._meta[i][1]*nxraw),"ty":float(cache._meta[i][0]*nxraw)}),"score":float(cache._meta[i][9]), "ctf":ctf})
		lsxout=None

		# Point locations
		np.savetxt(f"{options.path}/threed_{sn:02d}.txt",point.numpy,fmt="%0.4f",delimiter="\t")

		# show individual shifts at high verbosity
#		if options.verbose>2:
#			print("TYTX: ",(caches[stage[1]].tytx*nxraw).astype(np.int32))

	times.append(time.time())

	if options.score is not None:
		np.savetxt(options.score.replace(".lst","_score.txt"),qualities,fmt="%.6f")
		np.savetxt("scores_t.txt",qualities.transpose(),fmt="%.6f")

		lsxout=LSXFile(options.score)
		for i in range(nptcl):
			n,f,d=lsxin[selimg[i]]
			d["score"]=qualities[i][-1]
			lsxout[i]=n,f,d

	times.append(time.time())
	vol=point.volume(outsz,zmax).center_clip(outsz)
	vol=vol.emdata[0]
	if options.sym not in ("c1","C1","I","i"):
		if options.verbose>0 : print(f"Apply {options.sym} symmetry to map (not points)")
		vol.process_inplace("xform.applysym",{"sym":options.sym})
	times.append(time.time())
	vol["apix_x"]=apix*nxraw/outsz
	vol["apix_y"]=apix*nxraw/outsz
	vol["apix_z"]=apix*nxraw/outsz
	if options.ptcl3d_id is not None : vol["ptcl3d_id"]=options.ptcl3d_id
	vol.write_image(f"{options.path}/threed_{sn:02d}_unfilt.hdf:12",-1)
	if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
	vol.process_inplace("normalize.edgemean")
	times.append(time.time())
	vol.write_image(f"{options.path}/threed_{sn:02d}.hdf:12",-1)


	outf=open("frcstats.txt","w")
	for i in range(len(frchist[-1][0])):
		outf.write(f"{i}")
		for j in range(len(frchist)):
			try: outf.write(f"\t{frchist[j][0][i]:1.6f}\t{frchist[j][1][i]:1.6f}\t{frchist[j][2][i]:1.6f}")
			except:
				outf.write("\t0.0\t0.0\t0.0")
		outf.write("\n")

	# this is just to save some extra processing steps
	if options.fscdebug is not None:
		os.system(f'e2proc3d.py threed_{sn:02d}.hdf threed_{sn:02d}_fsc.txt --calcfsc {options.fscdebug}')

	times=np.array(times)
	times=times[1:]-times[:-1]
	if options.verbose>1 : print(times.astype(np.int32))

	E3end(llo)



if __name__ == '__main__':
	main()

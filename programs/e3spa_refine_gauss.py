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
import os

def main():

	usage="""e3spa_refine_gauss.py <projections>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--volfilt", type=float, help="Lowpass filter to apply to output volume, absolute, Nyquist=0.5", default=0.3)
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass, scaled with stage, default=500", default=500)
	parser.add_argument("--savesteps", action="store_true",help="Save the gaussian parameters for each refinement step, for debugging and demos")
	parser.add_argument("--fromscratch", action="store_true",help="Ignore orientations from input file and refine from scratch")
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	tf_set_device(dev=0,maxmem=options.gpuram)

	llo=E3init(sys.argv)
	rand = np.random.default_rng()

	if options.path == None:
		options.path=num_path_new("r3dgauss_")
		print("Writing to ",options.path)

	# store the input arguments in the refinement directory
	db = js_open_dict(f"{options.path}/0_refine_parms.json")
	db.update(vars(options))
	db["commandline"]=" ".join(sys.argv)
	db["timestamp"]=str(time.ctime())
	db.close()

	nptcl=EMUtil.get_image_count(args[0])
	nxraw=EMData(args[0],0,True)["nx"]
	nxrawm2=good_size_small(nxraw-2)
	apix=EMData(args[0],0,True)["apix_x"]

	if options.savesteps: 
		try: os.unlink("steps.hdf")
		except: pass

	if options.verbose: print(f"{nptcl} particles at {nxraw}^3")


	times=[time.time()]

	# Cache initialization
	if options.verbose: print("Caching particle data")
#	downs=sorted(set([s[1] for s in stages]))
	downs=sorted(set([min(i,nxrawm2) for i in (24,32,64,256,512)]))		# note that 24 is also used in reseeding
#	caches={down:StackCache(f"tmp_{os.getpid()}_{down}.cache",nptcl) for down in downs} 	# dictionary keyed by box size
	caches={down:StackCache(f"{options.path}/tmp_{down}.cache",nptcl) for down in downs} 	# dictionary keyed by box size
	fromscratch=options.fromscratch
	for i in range(0,nptcl,2500):
		if options.verbose>1: print(f"Caching {i}/{nptcl}")
		stk=EMStack2D(EMData.read_images(args[0],range(i,min(i+2500,nptcl))))
		orts,tytx=stk.orientations
		if orts is None or fromscratch:
			fromscratch=True
			tytx=np.zeros((stk.shape[0],2))
			orts=rand.random((stk.shape[0],3))-0.5
		else: tytx/=nxraw
		stkf=stk.do_fft()
		for down in downs:
			if down!=nxraw:
				stkfds=stkf.downsample(down)
				caches[down].write(stkfds,i,orts,tytx)
			else:
				caches[down].write(stkf,i,orts,tytx)

	# Reseed orientations for global search at low resolution
	tstorts=[]
	for x in np.arange(-0.5,0.5,0.04):
		for y in np.arange(-0.5,0.5,0.04):
			for z in np.arange(-0.5,0.5,0.04):
				if hypot(x,y,z)<=0.5: tstorts.append((x,y,z))
	tst_orts=Orientations(np.array(tstorts))

	# Forces all of the caches to share the same orientation information so we can update them simultaneously below (FRCs not jointly cached!)
	for down in downs[1:]:
		caches[down].orts=caches[downs[0]].orts
		caches[down].tytx=caches[downs[0]].tytx

	# definition of downsampling sequence for stages of refinement
	# 0)#ptcl, 1)downsample, 2)iter, 3)frc weight, 4)amp threshold, 5)replicate, 6)replicate spread, 7) gradient scale 8)frc loc threshold (9 disables gradient)
	# replication skipped in final stage. thresholds are mean+coef*std
	if fromscratch:
		print("Notice: refining from scratch without orientations")
		stages=[
			[200,   24,16,1.8, 0  ,2,.05, 2.0,9],
			[200,   24,16,1.8, 0  ,2,.05, 2.0,9],
			[200,   24,16,1.8,-1  ,2,.01, 2.0,-1],
			[5000,  24,16,1.8,-1  ,3,.1, 2.0,-1],
			[5000,  24,16,1.8, 0  ,1,.01, 2.0,-1],
			[5000,  32,16,1.5,-.5 ,2,.05,1.5,-3],
			[5000,  32,16,1.5,-1  ,3,.007,1.0,-2],
			[10000, 64,12,1.2,-1.5,3,.005,1.0,-3],
			[10000, 64,12,1.0,-2  ,3,.002,0.75,-3],
			[10000,256,12,1.2,-1.5,3,.005,1.0,-3],
			[10000,256,12,1.0,-2  ,3,.002,0.75,-3],
			[25000,512, 6,1.0,-2  ,1,.001,0.5,-3.0],
			[25000,512, 6,1.0,-2  ,1,.001,0.5,-3.0]
		]
	else:
		stages=[
			[1000,  24,16,1.8,-3  ,1,.01, 2.0, 9],
			[1000,  24,16,1.8, 0  ,2,.03, 1.5, 9],
			[1000,  24,16,1.8,-1  ,1,.01, 1.5, -3],
			[1000,  24,16,1.8, 0  ,2,.01, 1.5, -2],
			[2000,  32,16,1.5, 0  ,2,.02,1.5, -3],
			[2000,  32,16,1.5,-0.5,2,.01,1.25, -2],
			[5000,  64,24,1.2,-1  ,2,.005,1.0, -3],
			[10000,256,24,1.0,-1  ,2,.002,1.0,-3],
			[25000,512,12,1.0,-2  ,1,.001,1.0, -3]
		]

	for l in stages: l[1]=min(l[1],nxrawm2)		# make sure we aren't "upsampling"


	gaus=Gaussians()
	#Initialize Gaussians to random values with amplitudes over a narrow range
	rnd=tf.random.uniform((options.initgauss,4))     # specify the number of Gaussians to start with here
	rnd+=(-.5,-.5,-.5,10.0)
	gaus._data=rnd/(1.5,1.5,1.5,100.0)	# amplitudes set to ~1.0, positions random within 2/3 box size
	lsxin=LSXFile(args[0])

	times.append(time.time())
	ptcls=[]
	for sn,stage in enumerate(stages):
		if options.verbose: print(f"Stage {sn} - {local_datetime()}:")
		ccache=caches[stage[1]]

		nliststg=range(sn,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity
#		nliststg=range(0,nptcl,max(1,nptcl//stage[0]))		# all of the particles to use in the current stage, sn start gives some stochasticity
		norm=len(nliststg)//500+1

#	print(ptclsfds.shape,tytx.shape)
		
		if options.verbose: print(f"\tIterating Gaussian parms x{stage[2]} at size {stage[1]} with frc weight {stage[3]}\n    FRC\t\tshift_grad\tamp_grad")
		lqual=-1.0
		rstep=1.0
		for i in range(stage[2]):		# training epochs
			for j in range(0,len(nliststg),500):	# compute the gradient step piecewise due to memory limitations, 500 particles at a time
				ptclsfds,orts,tytx=ccache.read(nliststg[j:j+500])
				step0,qual0,shift0,sca0=gradient_step_gauss(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
				if j==0:
					step,qual,shift,sca=step0,qual0,shift0,sca0
				else:
					step+=step0
					qual+=qual0
					shift+=shift0
					sca+=sca0
			qual/=norm
			if qual<lqual: rstep/=2.0	# if we start falling or oscillating we reduce the step within the epoch
			step*=rstep/norm
			shift/=norm
			sca/=norm
			gaus.add_tensor(step)
			lqual=qual
			if options.savesteps: from_numpy(gaus.numpy).write_image(f"{options.path}/steps.hdf",-1)

			print(f"{i}: {qual:1.4f}\t{shift:1.4f}\t\t{sca:1.4f}")


		# reseed orientations of particles with low FRCs
		# we do this by finding the best orientation with fixed angular sampling and a fixed box size of 24
		nseeded=0
		if stage[8]<9:

			if options.verbose: print(f"Adjusting translational alignment of particles")
			for j in range(0,len(nliststg),500):	# compute the alignments piecewise due to memory limitations, 500 particles at a time
				ptclsfds,orts,tytx=ccache.read(nliststg[j:j+500])
				oldtytx=tytx
				tytx=ccf_step_align(gaus,ptclsfds,orts,tytx)
				ccache.tytx[nliststg[j:j+500],:2]=tytx
				dif=(tytx-oldtytx[:,:2])**2
				print(f"{j}-{j+500}: shift rmsd: {sqrt(tf.math.reduce_mean(dif)):.2f}")

			frcs=ccache.frcs			# not ideal, stealing the actual list from the object, but good enough for now
			lowfrc=frcs[frcs<1.5]
			if len(lowfrc)>0:
				frcm=np.mean(lowfrc)
				frcsg=np.std(lowfrc)
				reseed_idx=np.where(frcs<frcm+frcsg*stage[8])[0]			# [0] required because of odd np.where return
				nseeded=len(reseed_idx)

				if options.verbose: print(f"Making {len(tst_orts)} projections for reseeding")
				seedprojsf=gaus.project_simple(tst_orts,24).do_fft()		# fixed box size

				ptcls,tmpo,tmpxy=caches[24].read(reseed_idx)				# read all of the particle images we need to seed with new orientations, each is tiny with the fixed size of 24x24

				if options.verbose: print(f"Optimize {nseeded} orientations")
				for i in range(len(ptcls)):
					if options.verbose>1: print(f"{i}/{len(ptcls)}")
					ofrcs=tf_frc(seedprojsf.tensor,ptcls[i],-1)
					maxort=tf.argmax(ofrcs)			# best orientation for this particle
					ccache.orts[i]=tst_orts[maxort]
					#ccache.tytx[ii]=(0,0)			# just keep the current center?
				print(f"{nseeded} orts reseeded ({frcm+frcsg*stage[8]} thr)   {local_datetime()}")

			if options.verbose: print(f"\tIterating orientations parms x{stage[2]} with frc weight {stage[3]}\n    FRC\t\tort_grad\tcen_grad")
			fout=open(f"{options.path}/fscs.txt","w")
			for i in range(stage[2]):		# training epochs
				for j in range(0,len(nliststg),500):	# compute the gradient step piecewise due to memory limitations, 500 particles at a time
					ptclsfds,orts,tytx=ccache.read(nliststg[j:j+500])
					ortstep,dydxstep,qual0,ortstd0,dydxstd0,frcs=gradient_step_ort(gaus,ptclsfds,orts,tytx,stage[3],stage[7])
					if j==0:
						qual,ortstd,dydxstd=qual0,ortstd0,dydxstd0
					else:
						qual+=qual0
						ortstd+=ortstd0
						dydxstd+=dydxstd0

					#print(len(nliststg[j:j+500]),ortstep.shape,dydxstep.shape)
					ccache.add_orts(nliststg[j:j+500],ortstep,dydxstep)
					ccache.set_frcs(nliststg[j:j+500],frcs)
					for ii,n in enumerate(nliststg[j:j+500]): fout.write(f"{n}\t{frcs[ii]:1.4f}\t{orts[ii][0]:1.4f}\t{orts[ii][1]:1.4f}\t{orts[ii][2]:1.4f}\n")

				qual/=norm
				ortstd/=norm
				dydxstd/=norm

				print(f"{i}: {qual:1.4f}\t{ortstd:1.4f}\t\t{dydxstd:1.4f}")
		else: print("Skipping orientation gradient this step")

		# if options.savesteps:
		# 	vol=gaus.volume(nxraw)
		# 	vol.emdata[0].process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt})
		# 	vol.write_images(f"A_vol_opt_{sn}.hdf")

		# filter results and prepare for next stage
		g0=len(gaus)
		gaus.norm_filter(sig=stage[4],rad_downweight=0.33)	# remove gaussians below threshold
		g1=len(gaus)
		if stage[5]>0: gaus.replicate(stage[5],stage[6])	# make copies of gaussians with local perturbation
		g2=len(gaus)

		print(f"Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians  no orts reseeded   {local_datetime()}")




			# frcs=ccache.frcs			# not ideal, stealing the actual list from the object, but good enough for now
			# frcm=np.mean(frcs[frcs<1.5])
			# frcsg=np.std(frcs[frcs<1.5])
			# nseeded=0
			# for ii,f in enumerate(frcs):
			# 	if f<frcm+frcsg*stage[8]:
			# 		ccache.orts[ii]=rand.random((3,))-0.5
			# 		ccache.tytx[ii]=(0,0)
			# 		frcs[ii]=2.0
			# 		nseeded+=1
		times.append(time.time())
	
		# do this at the end of each stage in case of early termination

		# Particle orientations
		lsxout=LSXFile(f"{options.path}/ptcls_{sn:02d}.lst")
		for i in range(len(lsxin)):
			a,b,c=lsxin[i]
			lsxout[i]=(a,b,{"xform.projection":Transform({"type":"spinvec","v1":float(ccache.orts[i][0]),"v2":float(ccache.orts[i][1]),"v3":float(ccache.orts[i][2]),"tx":float(ccache.tytx[i][1]*nxraw),"ty":float(ccache.tytx[i][0]*nxraw)}),"frc":float(ccache.frcs[i])})
		lsxout=None

		# Gaussian locations
		out=open(f"{options.path}/threed_{sn:02d}.txt","w")
		for x,y,z,a in gaus.tensor: out.write(f"{x:1.5f}\t{y:1.5f}\t{z:1.5f}\t{a:1.3f}\n")

		# Filtered volume
		vol=gaus.volume(nxraw).emdata[0]
		vol["apix_x"]=apix
		vol["apix_y"]=apix
		vol["apix_z"]=apix
		vol.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.volfilt*min(stage[1],nxraw)/nxraw})
		vol.write_image(f"{options.path}/threed_{sn:02d}.hdf:12")


	times=np.array(times)
	times-=times[0]
	if options.verbose>1 : print(times)

	E3end(llo)

#@tf.function
def gradient_step_gauss(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0):
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

def gradient_step_ort(gaus,ptclsfds,orts,tytx,weight=1.0,relstep=1.0):
	"""Computes one gradient step on the orientation coordinates given a set of particle FFTs at the appropriate scale,
	computing FRC to axial Nyquist, with specified linear weighting factor (def 1.0). Linear weight goes from
	0-2. 1 is unweighted, >1 upweights low resolution, <1 upweights high resolution.
	returns ort_step, dydx_step, qual, ort_std, tytx_std and the 1-tensor with the per particle integrated FRCs for potential quality control
	step - one gradient step to be applied with (gaus.add_tensor)
	qual - mean frc
	shift - std of xyz shift gradient
	scale - std of amplitude gradient"""
	ny=ptclsfds.shape[1]

	with tf.GradientTape() as gt:
		gt.watch(orts.tensor)
		gt.watch(tytx)
		projs=gaus.project_simple(orts,ny,tytx=tytx)
		projsf=projs.do_fft()
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,3)	# specifying ny/2 radius explicitly so weight functions

	gradort,gradtytx=gt.gradient(frcs,(orts.tensor,tytx))
	qual=tf.math.reduce_mean(frcs)			# this is the average over all projections, not the average over frequency
	stdort=tf.math.reduce_std(gradort)		# orientation spinvec std
	stdtytx=tf.math.reduce_std(gradtytx)	# tytx std

	return (gradort*relstep/(stdort*1000),gradtytx*relstep/(stdtytx*1000),float(qual),float(stdort),float(stdtytx),frcs)

def ccf_step_align(gaus,ptclsfds,orts,tytx):
	"""Uses CCF to update all translational alignments in one step with CCF"""
	ny=ptclsfds.shape[1]
	projsf=gaus.project_simple(orts,ny,tytx=tytx).do_fft()
	newtytx=tf.cast(ptclsfds.align_translate(projsf),tf.float32)/float(-ny)
	return newtytx

if __name__ == '__main__':
	main()

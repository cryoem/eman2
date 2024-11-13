#!/usr/bin/env python
#
# Author: Anya Porter  10/10/2024
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

	usage="""e3tomogram_gauss.py <tiltseries>


	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--volout", type=str,help="Volume output file", default="threed.hdf")
	parser.add_argument("--gaussout", type=str,help="Gaussian list output file",default=None)
	parser.add_argument("--volfiltlp", type=float, help="Lowpass filter to apply to output volume in A, 0 disables, default=40", default=40)
	parser.add_argument("--volfilthp", type=float, help="Highpass filter to apply to output volume in A, 0 disables, default=2500", default=2500)
	parser.add_argument("--apix", type=float, help="A/pix override for raw data", default=-1)
	parser.add_argument("--thickness", type=float, help="For tomographic data specify the Z thickness in A to limit the reconstruction domain", default=-1)
	parser.add_argument("--tilesize",type=int,help="Controls how large the tiles are, default=1024", default=1024)
	parser.add_argument("--preclip",type=int,help="Trim the input images to the specified (square) box size in pixels", default=-1)
	parser.add_argument("--bin", type=int, help="Binning level for output file (will still use full resolution data to reconstruct", default=-1) 
	parser.add_argument("--initgauss",type=int,help="Gaussians in the first pass for each tile, scaled with stage, default=1000", default=1000)
	parser.add_argument("--savesteps", action="store_true",help="Save the gaussian parameters for each refinement step, for debugging and demos")
	parser.add_argument("--savetiles", action="store_true",help="Save the tiles as an hdf file outside of tmp. Currently only used for debugging but possibly could be used when multiple runs with same tiles")
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	tf_set_device(dev=0,maxmem=options.gpuram)

	llo=E3init(sys.argv)

	# Things to fix:
		# Figure out why the FRC was going negative for some tile--mostly at first few iterations so may be fine
		# FRC not currently remaining high for each tile--convergence issue and need more iterations or something else?
		# Edit number of Gaussians to be more reasonable: Upper limit: 1M Gaussians for 1kx1kx1k at 40A resolution (which is the default filtering)


	# Getting tiles and saving to temp files
	ntilts = EMUtil.get_image_count(args[0])
	nxraw, nyraw = EMData(args[0],0,True)["nx"], EMData(args[0],0,True)["ny"]
	if options.preclip>0: nxraw=nyraw=options.preclip
	if options.apix>0: apix=options.apix
	else: apix=EMData(args[0],0,True)["apix_x"]
	if options.thickness>0: zmax=options.thickness/(apix*nxraw*2.0)		# instead of +- 0.5 Z range, +- zmax range
	else: zmax=0.5

	# Define stages
	# 0) #ptcl, 1) downsample, 2) iter, 3) frc weight, 4) amp threshold, 5) replicate, 6) repl. spread, 7) step coef
	# replication skipped in final stage
	# Dr Ludtke's original from --tomo option of e3make3d_gauss.py
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
	# This treats tiles as fixed and downsamples the tiles, another thought would be to change the tiling for each stage

	step=options.tilesize//2

	outx = (nxraw//options.tilesize)*options.tilesize
	outy = (nyraw//options.tilesize)*options.tilesize
	xtiles = int(outx/step) + 1
	ytiles = int(outy/step) - 1 # Remove the tiles going off the edge along the non-tilt axis
	nptcls=int(ntilts*xtiles*ytiles)
	nxstep=int(outx/step/2)
	nystep=int(outy/step/2 - 1)

	if options.savesteps:
		try: os.unlink("steps.hdf")
		except: pass

	# Get information from info file
	js=js_open_dict(info_name(args[0]))
	try:
		ttparams = np.array(js["tlt_params"])
	except:
		traceback.print_exc()
		print(f"""\nERROR: tlt_params missing in {info_name(args[0])}. This will happen if a tomogram has not yet been constructed in EMAN. Support for this
			program being the first volume constructed should come later. For now reconstruct a tomogram (large bin should be fine) and then try again.""")
		sys.exit(1)
	try:
		cs=float(js["cs"])
	except:
		print(f"""\nWarning: Could not get Cs from {info_name(args[0])}""") #TODO: Change warning, give option to overwrite or something
	try:
		volt = float(js["voltage"])
	except:
		print(f"""\nWarning: Could not get Voltage from {info_name(args[0])}""") #TODO: Change warning, give option to overwrite or something
	try:
		defocus = np.array(js["defocus"])
	except:
		print(f"""\nWarning: Could not get defocus values from {info_name(args[0])}. Have you done CTF correction?""") # TODO: Change warning?
	try:
		phase = np.array(js["phase"])
		if min(phase)==max(phase): phase = [phase[0]]
	except:
		print(f"""\nWarning: Could not get phases from {info_name(args[0])}.""") #TODO: Change warning?
	js.close()

	# Loading tilts
	tilt_img = EMData(args[0], 0)
	if tilt_img["nz"]>1:
		tilt_imgs = [tilt_img.get_clip(Region(0,0,i, tilt_img["nx"], tilt_img["ny"],1)).copy() for i in range(img["nz"])]
	else:
		tilt_imgs = EMData.read_images(args[0])
	tilt_img=None

	for tilt in tilt_imgs:
		if options.preclip > 0: tilt = tilt.get_clip(Region(tilt["nx"]//2-nxraw//2, tilt["ny"]//2-nxraw//2, nxraw, nxraw), fill=0)
		tilt.process_inplace("threshold.clampminmax.nsigma",{"nsigma":10}) # removing x-ray pixels
		tilt.process_inplace("normalize.edgemean") # Normalizing
	ctf=EMAN2Ctf() # TODO: I checked and this works with import EMAN3 but will it eventually update to EMAN3Ctf or are we moving away from ctf object?
	ctf.from_dict({"defocus":1.0, "voltage":volt, "bfactor":0., "cs":cs, "ampcont":0, "apix":apix})
	ctf.set_phase(phase[0]*np.pi/180.)

	# Creating and caching tiles
	#TODO: Provide option to use saved tiling file from --savetiles so don't have to spend the time to extract them all?
	times= [time.time()]
	if options.verbose: print("Tiling and caching tilt series")
	downs=sorted(set([s[1] for s in stages]))
	caches={down:StackCache(f"tmp_{os.getpid()}_{down}.cache",nptcls) for down in downs}
	mindf=float('inf')
	maxdf=0
	for stepx in range(nxstep, -nxstep-1,-1):
		for stepy in range(nystep, -nystep-1,-1):
			full_tiles = []
			for i in range(ntilts):
				pos = [stepx*step, stepy*step, 0]
				pxf = get_xf_pos(ttparams[i],pos)
				tx = tilt_imgs[i]["nx"]//2 + pxf[0]
				ty = tilt_imgs[i]["ny"]//2 + pxf[1]
				m=tilt_imgs[i].get_clip(Region(int(tx) - step, int(ty) - step, options.tilesize, options.tilesize), fill=0) # Step used instead of recalculating tilesize//2
				m.mult(-1) # Make sure to invert contrast
				xform = Transform({"type":"xyz","ytilt":ttparams[i][3],"xtilt":ttparams[i][4],"ztilt":ttparams[i][2], "tx":tx-int(tx), "ty":ty-int(ty)}) # I skipped the dxfs part from e2spt_extract for now...was that important?
				m["xform.projection"]=xform
				rot=Transform({"type":"xyz","xtilt":float(ttparams[i][4]), "ytilt":float(ttparams[i][3])}) #TODO: Figure out why made separate transform instead of reusing xform
				p1=rot.transform(pos)
				pz=p1[2]*apix/10000. # Convert distance from center into defocus change due to tilt
				tilted_defocus = defocus[i]-pz
				ctf.defocus=tilted_defocus
				if len(phase) > 1:
					ctf.set_phase(phase[i]*np.pi/180.) # Avoid resetting phase if they are all the same
				m["ctf"]=ctf
				if tilted_defocus < mindf: mindf=tilted_defocus
				if tilted_defocus > maxdf: maxdf=tilted_defocus
				if options.ctf==0:
					fft1=m.do_fft()
					flipim=fft1.copy()
					ctf.compute_2d_complex(flipim, Ctf.CtfType.CTF_SIGN) 
					# TODO: See if this should be CTF_SIGN or if I should be doing full ctf correction here too instead of just phase flipping
					fft1.mult(flipim)
					m=fft1.do_ift()
				full_tiles.append(m)
			if options.savetiles:
				EMData.write_images("debug_tiling.hdf", full_tiles, ntilts*(-stepx+nxstep)*(ytiles)+ntilts*(-stepy+nystep))
				#Appears to tile correctly--starts in bottom left and goes up in columns moving right ^^^^ ->
			if options.verbose>1:
				print(f" Caching {ntilts*(-stepx+nxstep)*(ytiles)+ntilts*(-stepy+nystep)}/{nptcls}",end="\r",flush=True)
				sys.stdout.flush()
			stk=EMStack2D(full_tiles)
#			if options.preclip>0 : stk=stk.center_clip(options.preclip)
			#TODO: Should there be an option to clip here too that isn't tilesize?
			orts,tytx=stk.orientations
			try: tytx/= (nxraw, nyraw, 1) # TODO: Did I break this this by accepting non-square tiltseries?
						# I think probably not but maybe with tiling I did
			except: pass    # The try/except can probably be removed now--is just from when xform.projections was not set
			for im in stk.emdata: im.process_inplace("normalize.edgemean")
			stkf=stk.do_fft()
			for down in downs:
				stkfds=stkf.downsample(min(down,options.tilesize))
				caches[down].write(stkfds,ntilts*(-stepx+nxstep)*(ytiles)+ntilts*(-stepy+nystep),orts,tytx)
	tilt_imgs=None

	# Forces all of the caches to share the same orientation information so we can update them simultaneously below (FRCs not jointly cached!)
	for down in downs[1:]:
		caches[down].orts=caches[downs[0]].orts
		caches[down].tytx=caches[downs[0]].tytx


	if options.ctf > 0:
		ampcont=np.sin(phase[0]*np.pi/100)*100
		#TODO: When rewriting the create_CTF_stack in jax would it make sense to just pass phase instead of ampcont?
		dfstep = apix*apix/100
		boxlen = apix*options.tilesize*sqrt(3)
		df_buffer = (boxlen/2)/(dfstep*10000) + dfstep
		ctf_stack,dfstep=create_ctf_stack((mindf-df_buffer,maxdf+df_buffer),volt,cs,ampcont,options.tilesize,apix)

	if options.verbose>1: print("")

	# Initializing Gaussians to random values with amplitudes over a narrow range
	rnd = tf.random.uniform((options.initgauss*xtiles*ytiles, 4)) # all between 0-1
	rnd+= (-0.5, -0.5, -0.5, 2.0) # Coord between -.5 and .5, amp between 2 and 3
	rnd *= (xtiles//2, ytiles//2, 1., 1.) # Each tile has -.5 to .5 (with some offset n*0.5)
	rnd/= (0.9,0.9,1.0/zmax, 3.0) # Spread out points
	all_gaus = Gaussians()
	all_gaus._data = rnd
	cur_gaus= Gaussians()
	times.append(time.time())
	for sn,stage in enumerate(stages):
		all_gaus.coerce_numpy()
		for tnx in range(xtiles):
			xrange = (-0.5*(xtiles//2-tnx)-1, -0.5*(xtiles//2-tnx)+1)
			for tny in range(ytiles):
				yrange = (-0.5*(ytiles//2-tny)-1, -0.5*(ytiles//2-tny)+1)
				if options.verbose: print(f"Stage {sn}, Tile {tnx*ytiles+tny} - {local_datetime()}:")
				if options.verbose: print(f"\tIterating x{stage[2]} with frc weight {stage[3]}\n        FRC\t\tshift_grad\tamp_grad\timshift\tgrad_scale")
				cur_gaus._data = select_tile_gauss(all_gaus, xrange, yrange)
				cur_gaus.coerce_tensor()
				print(f"split types: {cur_gaus._data.dtype}, {all_gaus._data.dtype}")
				print(f"{len(cur_gaus)} Gaussians in tile {tnx*ytiles+tny}")
				print(f"{len(all_gaus)} Gaussians total")
				nliststg=range((tnx*ytiles+tny)*ntilts, (tnx*ytiles+tny)*ntilts+ntilts)
				imshift=0.0
				lqual=-1.0
				rstep=1.0
				for i in range(stage[2]):
					for j in range(0, len(nliststg), 500):
						ptclsfds, orts, tytx = caches[stage[1]].read(nliststg[j:j+500])
						if options.ctf==0:
							step0,qual0,shift0,sca0=gradient_step(cur_gaus,ptclsfds,orts,tytx,stage[3],stage[7])
							if j==0:
								step,qual,shift,sca=step0,qual0,shift0,sca0
							else:
								step+=step0
								qual+=qual0
								shift+=shift0
								sca+=sca0
						elif options.ctf==2:
							dsapix=apix*nxraw/ptclsfds.shape[1]
							step0,qual0,shift0,sca0=gradient_step_layered_ctf(cur_gaus,ptclsfds,orts,ctf_stack.downsample(ptclsfds.shape[1]),tytx,(mindf-df_buffer,maxdf+df_buffer),dfstep,dsapix,stage[3],stage[7])
							if j==0:
								step,qual,shift,sca=step0,qual0,shift0,sca0
							else:
								step+=step0
								qual+=qual0
								shift+=shift0
								sca+=sca
						elif options.ctf==1:
							step0,qual0,shift0,sca0=gradient_step_ctf(cur_gaus,ptclsfds,orts,ctf_stack.downsample(ptclsfds.shape[1]),tytx,(mindf-df_buffer,maxdf+df_buffer),dfstep,stage[3],stage[7])
							if j==0:
								step,qual,shift,sca=step0,qual0,shift0,sca0
							else:
								step+=step0
								qual+=qual0
								shift+=shift0
								sca+=sca
						# optimize gaussians and image shifts
						else:
							step0,stept0,qual0,shift0,sca0,imshift0=gradient_step_tytx(cur_gaus,ptclsfds,orts,tytx,stage[3],stage[7])
							if j==0:
								step,stept,qual,shift,sca,imshift=step0,stept0,qual0,shift0,sca0,imshift0
								caches[stage[1]].add_orts(nliststg[j:j+500],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
							else:
								step+=step0
								caches[stage[1]].add_orts(nliststg[j:j+500],None,stept0*rstep)	# we can immediately add the current 500 since it is per-particle
								qual+=qual0
								shift+=shift0
								sca+=sca0
								imshift+=imshift0
					# End of looping over j for memory
					norm=len(nliststg)//500+1
					qual/=norm
					if qual<lqual: rstep/=2.0	# if we start falling or oscillating we reduce the step within the epoch
					step*=rstep/norm
					shift/=norm
					sca/=norm
					imshift/=norm
					cur_gaus.add_tensor(step)
					lqual=qual
					if options.savesteps: from_numpy(cur_gaus.numpy).write_image("steps.hdf",-1)
					# Add debug_images to check geometry
					print(f"{i}: {qual:1.5f}\t{shift:1.5f}\t\t{sca:1.5f}\t{imshift:1.5f}\t{rstep:1.5f}")
					if qual>0.99: break
				# End of looping over epochs
				all_gaus._data = update_tomogram_gauss(cur_gaus, all_gaus, xrange, yrange)
		# End of looping over tiles
		# filter results and prepare for stage 2
		g0=len(all_gaus)
		all_gaus.norm_filter(sig=stage[4])		# gaussians outside the box may be important!
		g1=len(all_gaus)
		if stage[5]>0: all_gaus.replicate(stage[5],stage[6])
		g2=len(all_gaus)

		print(f"Stage {sn} complete: {g0} -> {g1} -> {g2} gaussians  {local_datetime()}")
		times.append(time.time())

		# do this at the end of each stage in case of early termination
		if options.gaussout is not None and g2 != 0:
			np.savetxt(options.gaussout,all_gaus.numpy,fmt="%0.4f",delimiter="\t")
			# out=open(options.gaussout,"w")
			# for x,y,z,a in gaus.tensor: out.write(f"{x:1.5f}\t{y:1.5f}\t{z:1.5f}\t{a:1.3f}\n")

		# show individual shifts at high verbosity
#		if options.verbose>2:
#			print("TYTX: ",(caches[stage[1]].tytx*nxraw).astype(np.int32))

	# End of looping over stages
	all_gaus._data=select_tile_gauss(all_gaus, (-1*(xtiles//2/2),xtiles//2/2), (-1*(ytiles//2/2), ytiles//2/2)) # Select only Gaussians in box for final volume (this was implicit with tensorflow not error at out of bounds index)
	times.append(time.time())
	vol=all_gaus.volume_tiled(outx,outy,options.tilesize,xtiles,ytiles,zmax)
	times.append(time.time())
	vol["apix_x"]=apix*nxraw/outx
	vol["apix_y"]=apix*nyraw/outy
	vol["apix_z"]=apix*nxraw/outx
	vol.write_image(options.volout.replace(".hdf","_unfilt.hdf"),0)
	if options.volfilthp>0: vol.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.volfilthp})
	if options.volfiltlp>0: vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.volfiltlp})
	if options.bin>0: 
		vol.process_inplace("math.meanshrink", {"n":options.bin}) # TODO: Ask if maybe this processing should go before the filters
		vol["apix_x"] = vol["apix_x"]*options.bin
		vol["apix_y"] = vol["apix_y"]*options.bin
		vol["apix_z"] = vol["apix_z"]*options.bin
	# Right now the unfilt will be at full resolution, maybe that's good maybe it should be downsampled too
	times.append(time.time())
	vol.write_image(options.volout,0)


	times=np.array(times)
	#times-=times[0]
	times=times[1:]-times[:-1]
	if options.verbose>1 : print(times.astype(np.int32))

	E3end(llo)


def get_xf_pos(tpm, pk):
	"""Function taken from e2tomogram.py (a version also exists in e2spt_extract.py), gets a 2D position on a tilt given a 3D location"""
	xf0 = Transform({"type":"xyz","xtilt":tpm[4],"ytilt":tpm[3],"ztilt":tpm[2],"tx":tpm[0], "ty":tpm[1]})
	p1 = xf0.transform([pk[0], pk[1], pk[2]]) # why doesn't it just say xf0.transform(pk)?
	return [p1[0], p1[1]] # Same question, why not return p1? is it an indexable class thats not a list?

def select_tile_gauss(all_gaus, xrange, yrange):
	"""Finds the subset of the Gaussians in all_gauss within xrange and yrange. 
	Returns a numpy array which are the Gaussians within xrange and yrange, translated to be centered around (0,0).
	Input:
		all_gaus: A Gaussian object with the all Gaussians in the tiltseries
		xrange: A (min, max) tuple of the x values we want to select
		yrange: A (min, max) tuple of the y values we want to select
	Output:
		A nx4 numpy array that should be set as the gaus._data for the tile's refinement as it includes all the Gaussians that satisfy xrange and yrange, translated to be centered around 0,0
	"""
	all_gaus.coerce_numpy()
	in_tile_gaus = all_gaus._data[(all_gaus._data[:,0] > xrange[0]) & (all_gaus._data[:,0] <= xrange[1]) & (all_gaus._data[:,1] > yrange[0]) & (all_gaus._data[:,1] <= yrange[1])]
	in_tile_gaus -= ((xrange[0]+xrange[1])/2, (yrange[0]+yrange[1])/2, 0., 0.)
	return in_tile_gaus
	# Overall tiltseries may not be square, don't put the calculation of where the tile is in this function give it a range (min,max) and 
	# Calculate that range in main code.
#	col=tn//tiles_per_side
#	row=tn-tiles_per_side*col
#	all_gaus.coerce_numpy()
#	in_tile = np.logical_and(np.logical_and(all_gaus._data[:,0]> -0.5*(tiles_per_side//2-col)-1, all_gaus._data[:,0] <= -0.5*(tiles_per_side//2-col)+1),
#				np.logical_and(all_gaus._data[:,1]> -0.5*(tiles_per_side//2-row)-1, all_gaus._data[:,1] <= -0.5*(tiles_per_side//2-row)+1))
#	gaus_in_tile = all_gaus[in_tile]
	# Just use 'and' not np.logical_and
#	Boolean index method using tensorflow
#	in_tile = tf.logical_and(tf.logical_and(all_gaus._data[:,0]> -0.5*(tiles_per_side//2-col)-1, all_gaus._data[:,0] <= -0.5*(tiles_per_side//2-col)+1),
#				tf.logical_and(all_gaus._data[:,1]> -0.5*(tiles_per_side//2-row)-1, all_gaus._data[:,1] <= -0.5*(tiles_per_side//2-row)+1))
#	gaus_in_tile += (0.5*(tiles_per_side//2-col),0.5*(tiles_per_side//2-row),0.,0.)
#	gaus_out_tile = all_gaus[np.logical_not(in_tile)]
#	return gaus_in_tile, gaus_out_tile

def update_tomogram_gauss(in_tile_gaus, all_gaus, xrange, yrange):
	""" Updates all_tile_gaus with the Gaussians from in_tile_gaus post-refinement. Only the Gaussians within the tile should be updated, not the ones outside as a buffer for the tilt
	Input:
		in_tile_gaus: A Gaussian object with the refined Gaussians from the tile specified by xrange and yrange
		all_gaus: A Gaussian object with the all Gaussians in the tiltseries
		xrange: A (min, max) tuple of the x values we want to select
		yrange: A (min, max) tuple of the y values we want to select
	Output:
		Returns: The new ._data for all_gaus"""
	in_tile_gaus.coerce_numpy()
	add_back = in_tile_gaus._data[(in_tile_gaus._data[:,0] > -0.5) & (in_tile_gaus._data[:,0] <= 0.5) & (in_tile_gaus._data[:,1] > -0.5) & (in_tile_gaus._data[:,1] <= 0.5)]
	add_back += ((xrange[0]+xrange[1])/2, (yrange[0]+yrange[1])/2, 0.,0.)
	add_to = all_gaus._data[(all_gaus._data[:,0] <= xrange[0]+0.5) | (all_gaus._data[:,0] > xrange[1]-0.5) | (all_gaus._data[:,1] <= yrange[0]+0.5) | (all_gaus._data[:,1] > yrange[1]-0.5)]
	print(add_back.shape, add_to.shape)
	return np.concatenate((add_to,add_back))
#	col = tn//tiles_per_side
#	row = tn-tiles_per_side*col
#	in_tile_gaus._data -= (0.5*(tiles_per_side//2-col), 0.5*(tiles_per_side//2-row), 0., 0.)
#	in_tile_gaus.coerce_numpy()
#	return np.concatenate((in_tile_gaus._data, out_tile_gaus._data))


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
#	dist_tuning=ny
	batch_size=16

	with tf.GradientTape() as gt:
		gt.watch(gaus.tensor)
		projs=gaus.project_simple(orts,ny,tytx=tytx)
		projsf=projs.do_fft()
		frcs=tf_frc(projsf.tensor,ptclsfds.tensor,ny//2,weight,2)	# specifying ny/2 radius explicitly so weight functions
#		dists=gaus.calc_distance(dist_tuning, batch_size)

#	grad=gt.gradient([frcs,dists],gaus._data)
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


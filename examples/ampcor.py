#!/usr/bin/env python

from EMAN2 import *
import numpy as np
import os
import sys
from scipy.optimize import minimize
from scipy.optimize import fmin
from joblib import Parallel, delayed 
from EMAN2_utils import *

def main():
	parser = EMArgumentParser()
	parser.add_argument("--tomogram", required=True, help="Path to tomogram.")
	parser.add_argument("--raw_tiltseries", required=True, help="Path to tiltseries.")	
	parser.add_argument("--filterto",default=0.125,type=float,help="cutoff_abs for lowpass filter applied to tilts and reconstruction projections. Useful to exclude high-resolution information we can't currently correct.")
	parser.add_argument("--niters", default=3,type=int,help="Number of iterations to run.")
	parser.add_argument("--threads", default=8,type=int,help="Number of threads to run in parallel on a single computer.")
	parser.add_argument("--outpath", default="ampcor", help="Name of output directory. Numbers will be appended to ensure each instance is written to a unique folder. Default is ampcor_XX.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()

	niters = options.niters

	info_file = info_name(options.raw_tiltseries)
	if os.path.isfile(info_file):
		jd = js_open_dict(info_file)
		try: 
			tiltparams = np.asarray(jd['tlt_params'])
			jd.close()
		except:
			jd.close()
			print("Could not locate tlt_params key in {}. Try running this through e2tomogram.py before proceeding.".format(info_file))
			sys.exit(1)
	else:
		print("Could not locate metadata file ({}) for {}. Exiting.".format(info_file,options.raw_tiltseries))
		sys.exit(1)

	if options.verbose: print("Loading input data...")

	if options.verbose: print("Reading in tomogram: {}".format(options.raw_tiltseries))
	r_phase=EMData(options.tomogram)
	rnx = r_phase["nx"]
	rny = r_phase["ny"]
	rnz = r_phase["nz"]
	apix = r_phase["apix_x"]

	## Not sure if this is correct, but it seems like it might help
	## amplitude contrast should result in img["minimum"]>=0, 
	## but phase effects could make that more complicated
	# if r_phase["minimum"]<0: 
	# 	r_phase -= r_phase["minimum"]

	if options.verbose: print("Reading in raw tiltseries: {}".format(options.raw_tiltseries))
	binfac = 1
	for f in [2,4,8,16]:
		if "_bin{}".format(f) in options.tomogram:
			binfac = f
	if binfac == 1: 
		tilts = EMData.read_images(options.raw_tiltseries)
	else:
		# not sure if this is read/write limited, but the fft resample process takes quite a bit of time on large images.
		ntlt = EMUtil.get_image_count(options.raw_tiltseries)
		tilts = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(get_binned_tilt) (i,binfac,options) for i in range(ntlt))

	outdir = make_path(options.outpath)

	# initial parameters for optimization. Likely not the best values.
	incid = []
	for t in tilts:
		inc = t+t["minimum"]
		incid.append(inc["maximum"])
	incid = np.max(incid) # incident intensity (approximated)

	alpha = 1.0
	beta = 0.1 # ...just something to start with

	x0 = [alpha/beta,incid] # initial guess for optimization

	if options.verbose: print("Generating initial projections of reconstructed tomogram")

	# Generate initial projections of reconstructed tomogram
	projs = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(make_projection) (tp,r_phase,binfac) for tp in tiltparams)
	for i,p in enumerate(projs):
		p.write_image("{}/projs_00.hdf".format(outdir),i)
	
	for itr in range(niters):

		thisi = str(itr).zfill(2)
		nexti = str(itr+1).zfill(2)
		
		if options.verbose: print("Generating amplitude contrast representation")
		
		# Get differences between reconstruction projections and original tilt series
		for i in range(len(tiltparams)):
			p = EMData("{}/projs_{}.hdf".format(outdir,thisi),i)
			dmask = p.process("threshold.notzero")
			p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})

			m = tilts[i].copy()
			reg = Region((m["nx"]-rnx)/2,(m["ny"]-2-rny)/2-1,rny,rny)
			m = m.get_clip(reg)

			d = p - m*dmask # not sure about this...
			d.process_inplace("math.squared")
			d.process_inplace("math.sqrt")
			d.process_inplace("math.log")
			d.process_inplace("threshold.belowtozero",{"minval":0.0})

			d.write_image("{}/logdiffs_{}.hdf".format(outdir,nexti),i)
		
		# Reconstruct amplitude contrast representation
		cmd = "e2make3dpar.py --input {inp} --output {out} --keep {keep} --apix {apix} --pad {pad} --sym {sym} --mode {mode} --altedgemask --threads {ncores} --outsize {nx} --clipz {nz}".format(
			inp="{}/logdiffs_{}.hdf".format(outdir,nexti), 
			out="{}/ramplitude_{}.hdf".format(outdir,nexti), 
			keep=1.0, apix=apix, pad=-1, sym="c1", 
			mode="gauss_2", ncores=options.threads, nx=rnx, nz=rnz)
		run(cmd,options)
		r_amp = EMData("{}/ramplitude_{}.hdf".format(outdir,nexti))

		if options.verbose: print("Optimizing parameters relating phase and amplitude contrast representations")

		# Optimize parameters relating phase and amplitude contrast representations		
		x_opt = fmin(cost, x0, disp=True,args=(r_phase,r_amp))
		print(x_opt)
		
		# construct phi volume from phase and amplitude contrast representations
		new_phi = (r_phase + r_amp * x_opt[0] /np.log(x_opt[1])) / 2.
		new_phi.write_image("{}/phivol_{}.hdf".format(outdir,nexti),0)
		
		if options.verbose: print("Generating corrected tiltseries")

		# generate corrected tiltseries
		projs = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(make_projection) (tp,new_phi,binfac) for tp in tiltparams)
		for i,p in enumerate(projs):
			p.write_image("{}/projs_{}.hdf".format(outdir,nexti),i)
		
		if options.verbose: print("Reconstructing corrected tiltseries")

		# Reconstruct new phase contrast representation
		cmd = "e2make3dpar.py --input {inp} --output {out} --keep {keep} --apix {apix} --pad {pad} --sym {sym} --altedgemask --mode {mode} --threads {ncores} --outsize {nx} --clipz {nz}".format(
			inp="{}/projs_{}.hdf".format(outdir,nexti),
			out="{}/threed_{}.hdf".format(outdir,nexti),
			keep=1.0,apix=apix,pad=-1,sym='c1',
			mode='gauss_2',ncores=options.threads,nx=rnx,nz=rnz)
		run(cmd,options)
		
		# define new phase contrast representation and tiltseries
		r_phase = EMData("{}/threed_{}.hdf".format(outdir,nexti))
		
		# probably DONT do this...
		#tilts = EMData.read_images("{}/projs_{}.hdf".format(outdir,nexti)) 

	print("DONE")

# Optimize parameters by minimizing ||R_phase+(a/b)*R_amp/log(I0)||^2
def cost(x,r_phase,r_amp):
	if x[1] <= 0.0 : return np.inf
	energy = r_phase - x[0] * ( r_amp / np.log(x[1]) )
	return energy["square_sum"]

def get_binned_tilt(i,binfac,options):
	t = EMData(options.raw_tiltseries,i)
	t.process_inplace("math.fft.resample",{"n":binfac})
	t.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	return t

def make_projection(param,threed,binfac):
	xf=Transform({"type":"xyz","tx":-param[0]/binfac,"ty":-param[1]/binfac,"ytilt":param[3],"xtilt":param[4]})
	p = threed.project("standard",xf)
	p.rotate(0,0,param[2])
	return p#-p["minimum"]

def run(cmd,options):
	if options.verbose: print("{}: {}".format(time.ctime(time.time()),cmd))
	ret=launch_childprocess(cmd)
	if ret != 0: print("Error running: {}".format(cmd))
	return

if __name__ == "__main__":
	main()
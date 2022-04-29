#!/usr/bin/env python
# Muyuan Chen 2020-05
from EMAN2 import *
import numpy as np

def main():
	
	usage=""" New single particle refinement routine. Still under construction. For simple usage,
	e2spa_refine.py --ptcl <particle list file> --ref <reference map> --res <inital resoution>
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcl", type=str,help="Input particle stack.", default="")
	parser.add_argument("--ref", type=str,help="Reference map. This will be scaled/clipped to match the particles automatically.", default="")
	parser.add_argument("--path", type=str,help="Path for refinement output files. default is r3d_xx", default=None)
	parser.add_argument("--parallel", type=str,help="Run in parallel, specify type:<option>=<value>:<option>=<value>.", default="thread:1")
	parser.add_argument("--sym", type=str,help="sym", default="c1")
	parser.add_argument("--res", type=float,help="The resolution that reference map is lowpass filtered to (with phase randomization) at the begining of the refinement. ", default=10)
	parser.add_argument("--keep", type=float,help="Fraction of best particles to keep in each iteration.", default=.9)
	parser.add_argument("--niter", type=int,help="Number of iterations. Default is 10.", default=10)
	parser.add_argument("--startiter", type=int,help="Start from a specified iteration in an existing refinement ", default=0)
	parser.add_argument("--setsf", type=str,help="Text file containing structure factor for map sharpening. Can be produced during CTF estimation, or from an existing high resolution map.", default=None)
	parser.add_argument("--tophat", type=str, default="local",help="Options for filtering maps. Run 'e2help.py tophat' for more information. Default=local.")
	parser.add_argument("--threads", type=int,help="Threads to use during postprocessing of 3d volumes", default=4)
	parser.add_argument("--mask", default=None, type=str,help="Specify a mask file for each iteration of refinement. Otherwise will generate mask automatically.")
	parser.add_argument("--compressbits", type=int,help="Bits to keep when writing images. 4 generally safe for raw data. 0-> true lossless (floating point). Default 6", default=6)
	parser.add_argument("--localsize",type=float,default=-1,help="Override the automatic local region size (in A) used for local resolution calculation and filtration.")
	parser.add_argument("--m3dthread",action="store_true", default=False ,help="do make3d in threading mode with shared memory. safer for large boxes")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	res=options.res
	sym=options.sym
		
	tophat=amask=localsize=sf=""
	
	npt=EMUtil.get_image_count(options.ptcl)
	if options.path==None: options.path=num_path_new("r3d_")
	
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_spa_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	if options.m3dthread:
		m3dpar=f" --threads {options.threads}"
	else:
		m3dpar=f" --parallel {options.parallel}"
		
	if options.startiter==0:
		
		r=1/res
		opt="--process filter.lowpass.gauss:cutoff_freq={:.4f} --process filter.lowpass.randomphase:cutoff_freq={:.4f}".format(r,r)
		
		er=EMData(options.ref,0,True)
		ep=EMData(options.ptcl,0,True)
		if abs(1-ep["apix_x"]/er["apix_x"])>0.01 or ep["ny"]!=er["ny"]:
			print("Reference-particle apix or box size mismatch. will scale/clip reference to match particles")
				
			rs=er["apix_x"]/ep["apix_x"]
			if rs>1.:
				opt+=" --clip {} --scale {} --process mask.soft:outer_radius=-1".format(ep["nx"], rs)
			else:
				opt+=" --scale {} --clip {} --process mask.soft:outer_radius=-1".format(rs, ep["nx"])
		
		#### prepare even/odd split
		lst=load_lst_params(options.ptcl)
		for i, l in enumerate(lst):
			l["class"]=i%2
		save_lst_params(lst,"{}/ptcls_00.lst".format(options.path))
		
		for eo in ["even","odd"]:
			run("e2proc3d.py {} {}/threed_00_{}.hdf {}".format(options.ref, options.path, eo, opt))
	
	
	for i in range(options.startiter, options.startiter+options.niter):
		
		run("e2spa_align.py --ptclin {pt}/ptcls_{i0:02d}.lst --ptclout {pt}/ptcls_{i1:02d}.lst --ref {pt}/threed_{i0:02d}.hdf --parallel {par} --sym {s} --maxres {rs:.2f} --goldcontinue --verbose {verbose}".format(pt=options.path, i0=i, i1=i+1, rs=res, s=sym, par=options.parallel, verbose=options.verbose))
			
		for eo in ["even","odd"]:
			run("e2spa_make3d.py --input {pt}/ptcls_{i1:02d}.lst --output {pt}/threed_{i1:02d}_{eo}.hdf --keep {kp} --sym {s} {par} --clsid {eo}".format(pt=options.path, i1=i+1, eo=eo, s=sym, par=m3dpar, kp=options.keep))

		if i==0:
			res/=2
		
		if i>0:
			tophat="--tophat {}".format(options.tophat)	
				
		if options.mask:
			amask="--mask {}".format(options.mask)

		if options.localsize>0:
			localsize="--localsize {}".format(options.localsize)
			
		if options.setsf:
			sf="--setsf {}".format(options.setsf)

		run("e2refine_postprocess.py --even {pt}/threed_{i1:02d}_even.hdf --sym {s} {sf} --restarget {rs:.1f} --threads {th} {top} {asize} {amask}".format(pt=options.path, i1=i+1, s=sym, sf=sf, rs=res*.8, th=options.threads, top=tophat,  amask=amask, asize=localsize))

		fsc=np.loadtxt("{}/fsc_masked_{:02d}.txt".format(options.path, i+1))
		fi=fsc[:,1]<0.2
		try: res=1./fsc[fi, 0][0]
		except:
			print("resolution approaching Nyquist !?")
			res=1.0/(2.0*er["apix_x"])
		res*=.8

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

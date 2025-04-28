#!/usr/bin/env python
# Muyuan Chen 2020-05
from EMAN2 import *
import numpy as np

def main():
	
	usage=""" New single particle refinement routine. Still under construction. For simple usage,
	e2spa_refine.py --ptcl <particle list file> --ref <reference map> --res <inital resoution>
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcl", type=str,help="Input particle stack.", default="", guitype='filebox', browser='EMSetsTable(withmodal=True,multiselect=False)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--ref", type=str,help="Reference map. This will be scaled/clipped to match the particles automatically.", default="", guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', filecheck=False, row=3, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--path", type=str,help="Path for refinement output files. default is r3d_xx", default=None)
	parser.add_argument("--parallel", type=str,help="Run in parallel, specify type:<option>=<value>:<option>=<value>.", default="thread:4", guitype='strbox', row=30, col=0, rowspan=1, colspan=2, mode="refinement[thread:4]")
	parser.add_argument("--sym", type=str,help="sym", default="c1", guitype='strbox', row=10, col=1, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--res", type=float,help="The resolution that reference map is lowpass filtered to (with phase randomization) at the begining of the refinement. ", default=16, guitype='floatbox', row=10, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--minrespx", type=int,default=-1, help="skip the first x pixels in fourier space", guitype='intbox', row=11, col=0, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--keep", type=float,help="Fraction of best particles to keep in each iteration.", default=.9, guitype='floatbox', row=12, col=1, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--niter", type=int,help="Number of iterations. Default is 8.", default=8, guitype='intbox', row=10, col=0, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--startiter", type=int,help="Start from a specified iteration in an existing refinement ", default=0)
	parser.add_argument("--setsf", type=str,help="Text file containing structure factor for map sharpening. Can be produced during CTF estimation, or from an existing high resolution map.", default=None)
	parser.add_argument("--tophat", type=str, default="global",help="Options for filtering maps. Run 'e2help.py tophat' for more information. Default=global (local is often a better choice)", guitype='strbox', row=12, col=0, rowspan=1, colspan=1, mode="refinement['global']")
	parser.add_argument("--gaussrecon",type=int,default=0,help="Use new e3make3d_gauss for 3-D reconstruction with N starting gaussians (1000 typ)")
	parser.add_argument("--threads", type=int,help="Threads to use during postprocessing of 3d volumes", default=4, guitype='intbox', row=30, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--mask", default=None, type=str,help="Specify a mask file for each iteration of refinement. Otherwise will generate mask automatically.", guitype='filebox', row=29, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--compressbits", type=int,help="Bits to keep when writing images. 4 generally safe for raw data. 0-> true lossless (floating point). Default 8", default=8)
	parser.add_argument("--localsize",type=float,default=-1,help="Override the automatic local region size (in A) used for local resolution calculation and filtration.")
	parser.add_argument("--m3dthread",action="store_true", default=False ,help="do make3d in threading mode with shared memory. safer for large boxes", guitype='boolbox', row=11, col=1, rowspan=1, colspan=1, mode="refinement[True]")
	parser.add_argument("--curve",action="store_true", default=False ,help="curve mode for filaments")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--localrefine", type=int, default=-1 ,help="local refinement. larger value correspond to smaller local region")

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
		
	er=EMData(options.ref,0,True)
	if options.startiter==0:
		
		r=1/res
		opt="--process filter.lowpass.gauss:cutoff_freq={:.4f} --process filter.lowpass.randomphase:cutoff_freq={:.4f}".format(r,r)
		
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
			
		run("e2proc3d.py {}/threed_00_even.hdf {}/threed_00.hdf --addfile {}/threed_00_odd.hdf --mult .5".format(options.path,options.path,options.path,))
	
	for i in range(options.startiter, options.startiter+options.niter):
		
		etc=""
		if options.curve: etc+=" --curve"
		
		if options.minrespx>0: mrp=f"--minrespx {options.minrespx}"
		else: mrp=""

		run("e2spa_align.py --ptclin {pt}/ptcls_{i0:02d}.lst --ptclout {pt}/ptcls_{i1:02d}.lst --ref {pt}/threed_{i0:02d}.hdf --parallel {par} --sym {s} --maxres {rs:.2f} --goldcontinue --verbose {verbose} --localrefine {lc} {etc} {mrp}".format(pt=options.path, i0=i, i1=i+1, rs=res, s=sym, par=options.parallel, verbose=options.verbose, lc=options.localrefine, etc=etc,mrp=mrp))
			
		if i==0:
			res/=1.5

		for ieo,eo in enumerate(["even","odd"]):
			if options.gaussrecon>0 :
				if options.parallel[:6]=="thread" and options.parallel.count(":")==2: 
					cache=options.parallel.split(":")[2]
				else: cache="." 
				run(f"e3make3d_gauss.py {options.path}/ptcls_{i+1:02d}.lst --volout {options.path}/threed_{i+1:02d}_{eo}.hdf:12 --gaussout {options.path}/threed_{i+1:02d}_{eo}.txt --sym {sym} --volfiltlp={res*0.75:.2f} --class {ieo} --cachepath {cache} --initgauss {options.gaussrecon}")
			else:
				run("e2spa_make3d.py --input {pt}/ptcls_{i1:02d}.lst --output {pt}/threed_{i1:02d}_{eo}.hdf --keep {kp} --sym {s} {par} --clsid {eo}".format(pt=options.path, i1=i+1, eo=eo, s=sym, par=m3dpar, kp=options.keep))
			run("e2proc3d.py {pt}/threed_{i1:02d}_{eo}.hdf {pt}/threed_raw_{eo}.hdf".format(pt=options.path, i1=i+1, eo=eo))

		
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
	

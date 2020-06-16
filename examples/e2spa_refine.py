#!/usr/bin/env python
# Muyuan Chen 2020-05
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcl", type=str,help="particle input", default="")
	parser.add_argument("--ref", type=str,help="reference", default="")
	parser.add_argument("--path", type=str,help="path", default="r3d_00")
	parser.add_argument("--emd", type=str,help="emd for comparison", default="")
	parser.add_argument("--parallel", type=str,help="", default="thread:1")
	parser.add_argument("--sym", type=str,help="sym", default="c1")
	parser.add_argument("--res", type=float,help="initial resolution", default=10)
	parser.add_argument("--keep", type=float,help="keep", default=.8)
	parser.add_argument("--startiter", type=int,help="iter", default=0)
	parser.add_argument("--niter", type=int,help="iter", default=10)
	parser.add_argument("--scipy",action="store_true",help="test scipy optimizer.",default=False)
	parser.add_argument("--seed",action="store_true",help=".",default=False)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	res=options.res
	sym=options.sym
	etcrefine=""
	if options.scipy:
		etcrefine+=" --scipytest"
	else:
		etcrefine+=" --fromscratch"
		
	if options.seed:
		etcrefine+=" --seedmap"
		
	etcpostp=""
	if options.emd:
		etcpostp+=" --emd {}".format(options.emd)
		
	npt=EMUtil.get_image_count(options.ptcl)
	
	if options.startiter==0:
		if not os.path.isdir(options.path):
			os.mkdir(options.path)
		
		r=1/res
		for i,eo in enumerate(["even","odd"]):
			run("e2proclst.py {} --create {}/ptcls_00_{}.lst --range {},{},2".format(options.ptcl, options.path, eo, i, npt))
			run("e2proc3d.py {} {}/threed_00_{}.hdf --process filter.lowpass.gauss:cutoff_freq={:.4f} --process filter.lowpass.randomphase:cutoff_freq={:.4f}".format(options.ref, options.path, eo, r,r))
	
	
	for i in range(options.startiter, options.startiter+options.niter):
		
		for eo in ["even","odd"]:
			run("e2spt_tiltrefine_oneiter.py --ptclin {pt}/ptcls_{i0:02d}_{eo}.lst --ptclout {pt}/ptcls_{i1:02d}_{eo}.lst --ref {pt}/threed_{i0:02d}_{eo}.hdf --threedout {pt}/threed_{i1:02d}_{eo}.hdf --threads 24 --parallel {par} --keep {kp} --sym {s} --maxres {rs:.2f} --minres -1 --padby 2 --nkeep 1 {etc}".format(pt=options.path, i0=i, i1=i+1, rs=res, eo=eo, s=sym, etc=etcrefine, par=options.parallel, kp=options.keep))

			
		run("e2refine_postprocess.py --even {pt}/threed_{i1:02d}_even.hdf --sym {s} --setsf strucfac.txt  {etc}".format(pt=options.path, i1=i+1, etc=etcpostp, s=sym))
		
		
		if res>0:
			fsc=np.loadtxt("{}/fsc_masked_{:02d}.txt".format(options.path, i+1))
			fi=fsc[:,1]<0.2
			res=1./fsc[fi, 0][0]
			if res<7:
				res=-1
				etcpostp+=" --tophat local"

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
#!/usr/bin/env python
# Muyuan Chen 2020-05
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcl", type=str,help="particle input", default="")
	parser.add_argument("--ncls", type=int,help="number of classes", default=2)
	parser.add_argument("--path", type=str,help="path", default="r3dcls_00")
	parser.add_argument("--parallel", type=str,help="", default="thread:1")
	parser.add_argument("--sym", type=str,help="sym", default="c1")
	parser.add_argument("--mask", type=str,help="mask file", default=None)
	parser.add_argument("--maxres", type=float,help="max resolution", default=10)
	parser.add_argument("--minres", type=float,help="min resolution", default=100)
	parser.add_argument("--niter", type=int,help="iter", default=10)
	parser.add_argument("--setsf", type=str,help="setsf", default=None)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	sym=options.sym
	path=options.path
	ncls=options.ncls
	if not os.path.isdir(path):
		os.mkdir(path)
		
	if options.setsf:
		setsf=" --setsf {}".format(options.setsf)
	else:
		setsf=""
	
	pinput="{}/ptcls_input.lst".format(path)
	if "even" in options.ptcl:
		print("Merging even/odd particle list")
		run("e2proclst.py {} {} --create {} --mergeeo".format(options.ptcl, options.ptcl.replace("even", "odd"), pinput))
	else:
		run("e2proclst.py {} {} --create {}".format(options.ptcl,options.ptcl, pinput))
	
	npt=EMUtil.get_image_count(pinput)
	cls=np.random.randint(0, ncls, npt)
	np.savetxt("{}/class_00.txt".format(path), cls)
	
	
	e=EMData(pinput, 0, True)
	apix=e["apix_x"]
	box=e["ny"]
	pad=good_size(box*1.5)
	etc=""
	if options.mask:
		etc+=" --mask {}".format(options.mask)
	
	for itr in range(options.niter+1):
		classify_list(pinput, cls, "{}/ptcls_{:02d}".format(path, itr))
		
		c0=cls.copy()
		scr=np.zeros((npt, ncls))
		refs=[]
		for ic in range(ncls):
			lname="{}/ptcls_{:02d}_{:02d}.lst".format(path, itr, ic)
			threed="{}/threed_{:02d}_{:02d}.hdf".format(path, itr, ic)
			refs.append(threed)
			
			run("e2spa_make3d.py --input {inp} --output {out} --keep 1 --sym {sm} --parallel {par}".format(inp=lname, out=threed, sm=sym, par=options.parallel))
			
			run("e2proc3d.py {} {} {} --process filter.lowpass.gauss:cutoff_freq={:.4f} --process normalize.edgemean".format(threed, threed, setsf, 1./options.maxres))
			
		
		sfile="{}/score_{:02d}.txt".format(path, itr)
		run("e2spa_classify.py {rf} --ptclin {inp} --output {out} --maxres {rsx:.1f} --minres {rsn:.1f} --parallel {par} {etc}".format(rf=' '.join(refs), inp=pinput, out=sfile, rsx=options.maxres, rsn=options.minres, par=options.parallel, etc=etc))
				
		scr=np.loadtxt(sfile)
		cls=np.argmin(scr, axis=1)
		np.savetxt("{}/class_{:02d}.txt".format(path, itr), cls)
		
		for ic in range(ncls):
			print("class {}:  {} particles".format(ic, np.sum(cls==ic)))
			
		
		
		print("iter {}, {:.1f}% particles change class".format(itr, 100*np.mean(c0!=cls)))
		
	
	ps=classify_list("r3dcls_02/ptcls_input.lst", cls, "r3dcls_02/ptcls_final", True)
	thd=[p.replace("ptcls_final","threed_final")[:-3]+"hdf" for p in ps]
	for pt,td, in zip(ps, thd):
		run("e2spa_make3d.py --input {inp} --output {out} --keep 1 --sym {sm} --parallel {par}".format(inp=pt, out=td, sm=sym, par=options.parallel))
		
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
def classify_list(lstname, cls, outprefix, spliteo=False):
	lst=LSXFile(lstname)
	ncls=np.max(cls)+1
	outnames=["{}_{:02d}.lst".format(outprefix, i) for i in range(ncls)]
	if spliteo:
		om=[]
		for o in outnames:
			om.append(o[:-4]+"_even.lst")
			om.append(o[:-4]+"_odd.lst")
		outnames=om	
		
	for o in outnames:
		if os.path.isfile(o):
			os.remove(o)

	lout=[LSXFile(o, False) for o in outnames]
	for i in range(lst.n):
		l=lst.read(i)
		ii=int(cls[i])
		
		if spliteo:
			ii=ii*2+i%2
			
		lout[ii].write(-1, l[0], l[1], l[2])

	lout=None
	return outnames
		
	
if __name__ == '__main__':
	main()
	
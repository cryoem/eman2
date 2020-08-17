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

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	sym=options.sym
	path=options.path
	ncls=options.ncls
	if not os.path.isdir(path):
		os.mkdir(path)
	
	pinput="{}/ptcls_input.lst".format(path)
	if "even" in options.ptcl:
		print("Merging even/odd particle list")
		run("e2proclst.py {} {} --create {} --mergeeo".format(options.ptcl, options.ptcl.replace("even", "odd"), pinput))
	
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
			
			run("e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pd} --outsize {bx} --apix {apx} --mode trilinear --keep 1 --sym {sm} --parallel {par}".format(inp=lname, out=threed, bx=box, pd=pad, apx=apix, sm=sym, par=options.parallel))
			
			run("e2proc3d.py {} {} --process filter.lowpass.gauss:cutoff_freq={:.4f} --process normalize.edgemean".format(threed, threed, 1./options.maxres))
			
		
		sfile="{}/score_{:02d}.txt".format(path, itr)
		run("e2spa_classify.py {rf} --ptclin {inp} --output {out} --maxres {rsx:.1f} --minres {rsn:.1f} --parallel {par} {etc}".format(rf=' '.join(refs), inp=pinput, out=sfile, rsx=options.maxres, rsn=options.minres, par=options.parallel, etc=etc))
				
		scr=np.loadtxt(sfile)
		cls=np.argmin(scr, axis=1)
		np.savetxt("{}/class_{:02d}.txt".format(path, itr), cls)
		
		for ic in range(ncls):
			print("class {}:  {} particles".format(ic, np.sum(cls==ic)))
			
		
		
		print("iter {}, {:.1f}% particles change class".format(itr, 100*np.mean(c0!=cls)))
		
	
		
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
def classify_list(lstname, cls, outprefix):
	lst=LSXFile(lstname)
	ncls=np.max(cls)+1
	outnames=["{}_{:02d}.lst".format(outprefix, i) for i in range(ncls)]
	for o in outnames:
		if os.path.isfile(o):
			os.remove(o)
	lout=[LSXFile(o, False) for o in outnames]
	for i in range(lst.n):
		l=lst.read(i)
		lout[int(cls[i])].write(-1, l[0], l[1], l[2])
		
	lout=None
		
	
if __name__ == '__main__':
	main()
	
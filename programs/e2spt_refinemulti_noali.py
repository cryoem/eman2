#!/usr/bin/env python
# Muyuan Chen 2020-03
from EMAN2 import *
import numpy as np
from shutil import copy2

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default="sptcls_00")
	parser.add_argument("--ncls", type=int,help="number of class", default=2)
	parser.add_argument("--niter", type=int,help="number of iterations", default=10)
	parser.add_argument("--maxres", type=float,help="max resolution in A", default=15)
	parser.add_argument("--parm", type=str,help="particle_parms_xx.json file from spt_refine", default="")
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:1")
	parser.add_argument("--mask", type=str,help="mask file", default="")
	parser.add_argument("--sym", type=str,help="", default="c1")
	parser.add_argument("--applysym", type=str,help="", default="c1")	
	parser.add_argument("--symcopy", action="store_true", default=False ,help="Copy each particle for each asymetrical unit. need a maskclass to focus on one unit. do not work with applysym")
	parser.add_argument("--randsym", action="store_true", default=False ,help="randomize symmetry unit assignment at the begining. seems to be helpful for high symmetry without obvious asymetrical features..")


	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	path=options.path
	ncls=options.ncls
	
	if not os.path.isdir(path):
		os.mkdir(path)
		
	if options.randsym and options.sym!="c1":
		js=dict(js_open_dict(options.parm)).copy()
		x0=Transform()
		nsym=x0.get_nsym(options.sym)
		out={}
		for k in js.keys(): 
			d=dict(js[k]).copy()
			xf=Transform(d["xform.align3d"]) 
			x=x0.get_sym(options.sym, np.random.randint(nsym)) 
			d["xform.align3d"]=x*xf 
			out[k]=d
			
		jo=js_open_dict("{}/particle_parms_00.json".format(path))
		jo.update(out)
		jo.close()
		
	else:
		copy2(options.parm, "{}/particle_parms_00.json".format(path))
	
	opt=""
	if options.mask!="":
		opt+=" --maskclass {}".format(options.mask)
	if options.symcopy:
		opt+=" --symcopy"

	for itr in range(options.niter+1):
		if itr==0:
			cmd="e2spt_average_multi.py --path {pp} --iter 0 --sym {sm} --applysym {asm} --noali --simthr 1 --simthr2 1 --threads 24 --parallel {par} --randnclass {c} --nolstout {opt}".format(pp=path, c=ncls, par=options.parallel, sm=options.sym, asm=options.applysym, opt=opt)
			run(cmd)
		else:
			tds=["{}/threed_{:02d}_{:02d}.hdf".format(path, itr-1, c) for c in range(ncls)]
			tds=' '.join(tds)
			cmd="e2spt_average_multi.py {mp} --path {pp} --iter {it} --sym {sm} --applysym {asm} --noali --simthr 1 --simthr2 1 --threads 24 --parallel {par} --maxres {res:.1f}  --nolstout {opt}".format(mp=tds, pp=path, it=itr, par=options.parallel, res=options.maxres, opt=opt, sm=options.sym, asm=options.applysym)
			run(cmd)
			
			mt0=np.loadtxt("{}/avg_multi_{:02d}.txt".format(path, itr-1))
			mt1=np.loadtxt("{}/avg_multi_{:02d}.txt".format(path, itr))
			chg=np.mean(np.any(mt0[:,2:]!=mt1[:,2:], axis=1))*100
			print("iter {}, {:.1f}% particles change class".format(itr, chg))
		
		for i in range(ncls):
			cmd="e2proc3d.py {pp}/threed_{it:02d}_{cl:02d}.hdf {pp}/threed_{it:02d}_{cl:02d}.hdf --setsf sf.txt --process filter.lowpass.gauss:cutoff_freq={res:.3f} --process normalize.edgemean --process mask.soft:outer_radius=-10:width=10".format(it=itr, cl=i, pp=path, res=1./options.maxres)
			run(cmd)
			
	itr=options.niter
	mt=np.loadtxt("{}/avg_multi_{:02d}.txt".format(path, itr))
	js=js_open_dict("{}/particle_parms_00.json".format(path))
	s=Transform()
	for i in range(ncls):
		jout={}
		for k in list(js.keys()):
			src,ii=eval(k)
			ic=np.where(mt[:,0]==ii)[0][0]
			if mt[ic,2]==i:
				jout[k]=js[k]
				xf=js[k]["xform.align3d"]
				x=s.get_sym(options.sym, int(mt[ic,3]))
				jout[k]["xform.align3d"]=x*xf
				
		fm="{}/particle_parms_{:02d}_cls_{:02d}.json".format(path, itr, i)
		jo=js_open_dict(fm)
		jo.update(jout)
		jo.close()
		print("writing to {}".format(fm))

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

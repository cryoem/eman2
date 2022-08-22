#!/usr/bin/env python
# Muyuan Chen 2020-03
from EMAN2 import *
import numpy as np
from sklearn import cluster,mixture
from sklearn.decomposition import PCA

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--pts", type=str,help="point input", default="")
	parser.add_argument("--pcaout", type=str,help="pca output", default="")
	parser.add_argument("--ptclsin", type=str,help="ptcl input", default="")
	parser.add_argument("--ptclsout", type=str,help="ptcl out suffix", default="")
	parser.add_argument("--pad", type=int,help="pad for make3d", default=-1)
	parser.add_argument("--ncls", type=int,help="number of classes", default=3)
	parser.add_argument("--nbasis", type=int,help="PCA dimensionality", default=2)
	parser.add_argument("--width", type=float,help="width of the vector. 1 covers all points. default 0.95", default=.95)
	parser.add_argument("--setsf", type=str,help="setsf", default="")
	parser.add_argument("--mode", type=str,help="classify/regress", default="classify")
	parser.add_argument("--axis", type=str,help="axis for regress. one number for a line, and two numbers separated by comma to draw circles.", default='0')
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--nptcl", type=int,help="number of particles per class in regress mode", default=2000)
	parser.add_argument("--threads", default=12,type=int,help="Number of threads to run in parallel on a single computer. This is the only parallelism supported by e2make3dpar")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	pts=np.loadtxt(options.pts)
	
	pca=PCA(options.nbasis)
	p2=pca.fit_transform(pts[:,1:])
	if options.pcaout:
		np.savetxt(options.pcaout, p2)

	if options.mode=="classify":
		clust=cluster.KMeans(options.ncls)
		lbs=clust.fit_predict(p2)
		lbunq=np.unique(lbs)
		
	else:
		axis=[int(i) for i in options.axis.split(',')]
		print('regress along axis', axis)
		if len(axis)==1:
			p=p2[:,axis[0]]
			rg=np.arange(options.ncls)
			rg=rg/np.max(rg)-.5
			mx=2*np.sort(abs(p))[int(len(p)*options.width)]
			rg=rg*mx+np.mean(p)
			print(rg)
			
		else:
			p=np.linalg.norm(p2[:, axis], axis=1)
			mx=np.sort(abs(p))[int(len(p)*options.width)]
			t=np.arange(options.ncls)/options.ncls
			t=t*np.pi*2
			rg=np.array([np.cos(t), np.sin(t)]).T
			rg*=mx
			print(rg)
			
		
		
	onames=[]
	fname=options.ptclsin
	lstinp=fname.endswith(".lst")
	
	lin=load_lst_params(fname)
		
	lout=[]
	for j in range(options.ncls):
		
		
		if options.mode=="classify":
			l=lbunq[j]
			ii=(lbs==l)
			print(j, np.sum(ii))
		else:
			if len(axis)==1:
				d=abs(p2[:,axis[0]]-rg[j])
			else:
				d=np.linalg.norm(p2[:,axis]-rg[j], axis=1)
				
			ii=np.argsort(d)[:options.nptcl]
			print(j, rg[j], d[ii[-1]])
		
		idx=pts[ii,0].astype(int)
		lo=[lin[i].copy() for i in idx]
		for l in lo:
			l["class"]=j
		
		lout.extend(lo)
	
	save_lst_params(lout, options.ptclsout)
	
	e=EMData(fname, 0, True)
	if options.pad<1: options.pad=good_size(e["nx"]*1.25)
	if options.setsf:
		options.setsf=" --setsf "+options.setsf
	
	name3d=options.ptclsout[:-4]+".hdf"
	if os.path.isfile(name3d):
		os.remove(name3d)
		
	for j in range(options.ncls):
		t="tmp_classify.hdf"
		cmd="e2spa_make3d.py --input {} --output {} --pad {} --keep 1 --threads {} {} --sym {} --clsid {}".format(options.ptclsout, t, options.pad, options.threads, options.setsf, options.sym, j)
		launch_childprocess(cmd)
		e=EMData(t)
		e.write_image(name3d,-1)
	
	os.remove(t)
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

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
	parser.add_argument("--setsf", type=str,help="setsf", default="")
	parser.add_argument("--mode", type=str,help="classify/regress", default="classify")
	parser.add_argument("--axis", type=str,help="axis for regress. one number for a line, and two numbers separated by comma to draw circles.", default='0')
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--nptcl", type=int,help="number of particles per class in regress mode", default=2000)
	parser.add_argument("--threads", default=12,type=int,help="Number of threads to run in parallel on a single computer. This is the only parallelism supported by e2make3dpar")
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
			mx=2*np.sort(abs(p))[int(len(p)*.9)]
			rg=rg*mx+np.mean(p)
			print(rg)
			
		else:
			p=np.linalg.norm(p2[:, axis], axis=1)
			mx=np.sort(abs(p))[int(len(p)*.9)]
			t=np.arange(options.ncls)/options.ncls
			t=t*np.pi*2
			rg=np.array([np.cos(t), np.sin(t)]).T
			rg*=mx
			print(rg)
			
		
		
	onames=[]
	fname=options.ptclsin
	lstinp=fname.endswith(".lst")
	
	if lstinp:
		lin=LSXFile(fname)
		
	for j in range(options.ncls):
		
		onm="{}_{:02d}.lst".format(options.ptclsout, j)
		
		if options.mode=="classify":
			l=lbunq[j]
			ii=(lbs==l)
			print(onm, np.sum(ii))
		else:
			if len(axis)==1:
				d=abs(p2[:,axis[0]]-rg[j])
				ii=np.argsort(d)[:options.nptcl]
				print(onm, rg[j], d[ii[-1]])
			else:
				d=np.linalg.norm(p2[:,axis]-rg[j], axis=1)
				ii=np.argsort(d)[:options.nptcl]
				print(onm, rg[j], d[ii[-1]])
		
		idx=pts[ii,0].astype(int)
		
		if os.path.isfile(onm):
			os.remove(onm)
		lout=LSXFile(onm, False)
		for i in idx:
			if lstinp:
				l=lin.read(i)
				lout.write(-1, l[0], l[1],l[2])
			else:
				lout.write(-1, i, fname)
			
		lout=None
		onames.append(onm)
	
	e=EMData(fname, 0, True)
	if options.pad<1: options.pad=good_size(e["nx"]*1.25)
	if options.setsf:
		options.setsf=" --setsf "+options.setsf
	
	for o in onames:
		t=o[:-3]+"hdf"
		print(o,t)
		cmd="e2spa_make3d.py --input {} --output {} --pad {} --keep 1 --threads {} {} --sym {}".format(o,t, options.pad, options.threads, options.setsf, options.sym)
		launch_childprocess(cmd)
	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

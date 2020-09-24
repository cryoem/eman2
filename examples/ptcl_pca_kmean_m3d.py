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
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	pts=np.loadtxt(options.pts)
	
	pca=PCA(options.nbasis)
	p2=pca.fit_transform(pts[:,1:])
	if options.pcaout:
		np.savetxt(options.pcaout, p2)

	clust=cluster.KMeans(options.ncls)

	lbs=clust.fit_predict(p2)
	onames=[]
	fname=options.ptclsin
	lstinp=fname.endswith(".lst")
	
	if lstinp:
		lin=LSXFile("gmm_03/ptcls_set2.lst")
		
	for j,l in enumerate(np.unique(lbs)):
		ii=(lbs==l)
		idx=pts[ii,0].astype(int)
		
		onm="{}_{:02d}.lst".format(options.ptclsout, j)
		print(len(idx),onm)
		if os.path.isfile(onm):
			os.remove(onm)
		lout=LSXFile(onm, False)
		for i in idx:
			if lstinp:
				l=lin.read(i)
				lout.write(-1, l[0], l[1])
			else:
				lout.write(-1, i, fname)
			
		lout=None
		onames.append(onm)
	
	e=EMData(fname, 0, True)
	if options.pad<1: options.pad=good_size(e["nx"]*1.25)
	if options.setsf:
		options.setsf=" --setsf "+options.setsf
	
	for o in onames:
		t=o.replace("ptcl", "threed")[:-3]+"hdf"
		print(o,t)
		cmd="e2make3dpar.py --input {} --output {} --pad {} --mode trilinear --no_wt --keep 1 --threads 12 {}".format(o,t, options.pad, options.setsf)
		launch_childprocess(cmd)
	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

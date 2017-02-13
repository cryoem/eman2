#!/usr/bin/env python
# Muyuan Chen 2017-02

import numpy as np
from EMAN2 import *
import sklearn.manifold as manifold
import sklearn.decomposition as decomp

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--ncomponents", type=int,help="number of output components", default=2)
	parser.add_argument("--perplexity", type=int,help="perplexity for TSNE", default=300)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.path==None:
		folder=[f for f in os.listdir('.') if f.startswith("refine_") and os.path.isdir(f)]
		options.path=np.sort(folder)[-1]
	options.path=options.path.strip("\\/")
	print "Working on refinement {}...".format(options.path)
	
	
	pjsimx=[f for f in os.listdir(options.path) if f.startswith("proj_simmx_") and f.endswith("even.hdf")]
	pjsimx=os.path.join(options.path, np.sort(pjsimx)[-1])
	print "Using {}...".format(pjsimx)
	e=EMData(pjsimx)
	simx=e.numpy().copy()
	
	print "Doing PCA..."
	dc=decomp.PCA(n_components=options.ncomponents)
	dcout=dc.fit_transform(simx)
	write_out(dcout,"pca_projs_{}.txt".format(options.path))
	
	
	print "Doing TSNE..."
	mani=manifold.TSNE(n_components=options.ncomponents,verbose=3, perplexity=options.perplexity)
	mnout=mani.fit_transform(simx) 
	write_out(mnout,"tsne_projs_{}.txt".format(options.path))

	E2end(logid)
		
def write_out(data,fname=None,projname="refine_final/projections_03_even.hdf"):
	n=len(data[0])
	if fname:
		f=open(fname,'w')
		
	for di,d in enumerate(data):
		s=""
		for i in range(n):
			s+="{}\t".format(d[i])
		s+="# {};{}\n".format(di,projname)
		if fname:
			f.write(s)
		else:
			print s[:-1]
		
	if fname:
		f.close()
		return
	else:
		return s



	
if __name__ == '__main__':
	main()
	
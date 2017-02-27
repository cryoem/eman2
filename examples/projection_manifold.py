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
	parser.add_argument("--dopca", type=int,help="perform PCA for comparison", default=1)
	parser.add_argument("--perplexity", type=int,help="perplexity for TSNE", default=300)
	parser.add_argument("--mode", type=int,help="choose from: 0: TSNE, 1: Isomap, 2: LocallyLinearEmbedding, 3:SpectralEmbedding, default is TSNE", default=0)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.path==None:
		folder=[f for f in os.listdir('.') if f.startswith("refine_") and os.path.isdir(f)]
		options.path=np.sort(folder)[-1]
	options.path=options.path.strip("\\/")
	print "Working on refinement {}...".format(options.path)
	
	
	pjsimx=[f for f in os.listdir(options.path) if f.startswith("proj_simmx_") and f.endswith("even.hdf")]
	pjsimx=os.path.join(options.path, np.sort(pjsimx)[-1])
	projname=pjsimx.replace("proj_simmx_", "projections_")
	print "Using {}...".format(pjsimx)
	e=EMData(pjsimx)
	simx=e.numpy().copy()
	
	if options.dopca:
		print "Doing PCA..."
		dc=decomp.PCA(n_components=options.ncomponents)
		dcout=dc.fit_transform(simx)
		write_out(dcout,"pca_projs_{}.txt".format(options.path), projname)
	
	mode_dict={0: "TSNE", 1: "Isomap", 2: "LocallyLinearEmbedding", 3:"SpectralEmbedding"}
	mode=mode_dict[options.mode]
	if options.mode==0:
		print "Doing TSNE..."
		mani=manifold.TSNE(n_components=options.ncomponents,verbose=3, perplexity=options.perplexity)
	else:
		print "Doing {}...".format(mode)
		mani=eval("manifold.{}(n_components={:d})".format(mode, options.ncomponents))
		
	mnout=mani.fit_transform(simx) 
	write_out(mnout,"{}_projs_{}.txt".format(mode,options.path), projname)

	E2end(logid)
		
def write_out(data,fname=None,projname=None):
	n=len(data[0])
	if fname:
		f=open(fname,'w')
		
	for di,d in enumerate(data):
		s=""
		for i in range(n):
			s+="{}\t".format(d[i])
		
		if projname:
			s+="# {};{}\n".format(di,projname)
		else:
			s+="\n"
		if fname:
			f.write(s)
		else:
			print s[:-1]
		
	if fname:
		f.close()
		print "Output written to {}.".format(fname)
		return
	else:
		return s



	
if __name__ == '__main__':
	main()
	
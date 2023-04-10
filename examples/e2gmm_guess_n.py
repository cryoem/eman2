#!/usr/bin/env python
# Muyuan Chen 2023-03
from EMAN2 import *
from EMAN2_utils import *
from scipy.spatial import KDTree
import numpy as np
from sklearn.cluster import KMeans

def main():

	usage="""
	Guess the number of Gaussian needed to represent a given volume. 
	e2gmm_guess_n.py threed_xx.hdf --thr 3 --maxres 3
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--thr", type=float,help="threshold", default=-1)
	parser.add_argument("--maxres", type=float,help="resolution", default=3)
	parser.add_argument("--minres", type=float,help="resolution", default=100)
	parser.add_argument("--startn", type=int,help="", default=1000)
	parser.add_argument("--maxn", type=int,help="", default=20000)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	fname=args[0]
	if options.thr<0:
		tmp="tmp_input.hdf"
		run(f"e2proc3d.py {fname} {tmp} --process normalize.edgemean")
		fname=tmp
		options.thr=4
		
	e=EMData(fname)
	img=e.numpy().copy()
	pts=np.array(np.where(img>options.thr)).T
	
	val=img[pts[:,0], pts[:,1], pts[:,2]]
	#pts[:,1:]=e["nx"]-pts[:,1:]
	pts=pts[:,::-1]
	pts=pts*e["apix_x"]
	pts=pts[np.argsort(-val)]
	
	tree=KDTree(pts)
	d=tree.query(pts, k=2)[0][:,1]
	dm=np.sort(d)[len(d)//2]
	print(np.sort(d))
	print(len(pts), np.mean(d), np.median(d))
	#print("{:<5} {:.4f}".format(len(pts), dm))
	res=max(options.maxres/2., 1.3) ## shouldn't have atoms that close even for atomic model
	
	tokeep=np.ones(len(pts), dtype=bool)
	for i in range(len(pts)):
		if tokeep[i]:
			k=tree.query_ball_point(pts[i], res)
			tokeep[k]=False
			tokeep[i]=True
			
	print(np.sum(tokeep))
	
	options.maxn=min(options.maxn, len(pts))
	nrng=np.arange(options.startn, options.maxn+1, 1000)
	
	print("N:    dist")
	for n in nrng:
		#launch_childprocess(f"e2segment3d.py {fname} --pdb threed_seg.pdb --process=segment.kmeans:nseg={n}:thr={options.thr}:minsegsep=1")
		#p=pdb2numpy("threed_seg.pdb")
		p0=pts.copy()
		np.random.shuffle(p0)
		p0=p0[:n]
		#p0+=np.random.randn(len(p0),3)*.1
		km=KMeans(n,init=p0, n_init=1,max_iter=100)
		km.fit(pts)
		p=km.cluster_centers_
		
		tree=KDTree(p)
		d=tree.query(p, k=2)[0][:,1]
		dm=np.sort(d)[len(d)//3]
		print(np.sort(d))
		print("{:<5} {:.4f}".format(n, dm))
		numpy2pdb(p, "threed_seg.pdb")
		if dm<res:
			print(f"stop. using N={n}")
			break
			
	run(f"e2project3d.py {fname} --outfile tmp_projection.hdf --orientgen=eman:delta=4 --parallel=thread:16")
	run(f"e2gmm_refine_new.py --ptclsin tmp_projection.hdf --model threed_seg.pdb --maxres {options.maxres} --minres {options.minres} --modelout model_gmm.txt --niter 50 --trainmodel --evalmodel tmp_model_projections.hdf --learnrate 1e-5")
	run(f"e2spa_make3d.py --input tmp_model_projections.hdf --output tmp_avg.hdf --thread 32")
	run(f"e2proc3d.py tmp_avg.hdf model_avg.hdf --process mask.soft:outer_radius=-16 --matchto {fname}")
	run(f"e2proc3d.py model_avg.hdf model_fsc.txt --calcfsc {fname}")

	print("final pdb model: threed_seg.pdb")
	print("final GMM in text file: model_gmm.txt")
	print("final GMM in density map: model_avg.hdf")
	print("map-model FSC: model_fsc.txt")
	E2end(logid)
	
	
if __name__ == '__main__':
	main()

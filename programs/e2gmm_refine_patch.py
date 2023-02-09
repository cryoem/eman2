#!/usr/bin/env python
# Muyuan Chen 2022-12
from EMAN2 import *
import numpy as np
from sklearn.cluster import KMeans

	
def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--niter", type=int,help="iteration", default=3)
	parser.add_argument("--startres", type=float,help="starting resolution", default=3)
	parser.add_argument("--npatch", type=int,help="number of patches", default=8)
	parser.add_argument("--maxres", type=float,help="max resolution", default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	fname=args[0]
	oldpath=os.path.dirname(fname)
	
	c=[i for i in fname if i.isdigit()]
	olditer=int(''.join(c[-2:]))
	print(f"refinement path {oldpath}, iteration {olditer}.")
	p00=np.loadtxt(f"{oldpath}/model_{olditer:02d}_even.txt")
	
	
	if options.path==None: 
		path=options.path=num_path_new("gmm_")
	
	p=p00
	pn=options.npatch
	km=KMeans(pn,max_iter=30)
	km.fit(p[:,:3])
	pc=km.cluster_centers_
	lb=km.labels_
	path=options.path
	
	msk0=EMData(f"{oldpath}/mask.hdf")
	msk0.write_image(f"{path}/mask_00.hdf")
	msk=msk0.copy()
	
	etc=""
	if options.maxres>0:
		etc+=f" --maxres {options.maxres}"
	
	for ci in range(pn):
		c=pc[ci].copy()
		r=np.linalg.norm(p[lb==ci,:3]-c, axis=1)
		r=np.max(r)
		
		msk.to_one()
		c=(c*msk["nx"]).tolist()
		msk.process_inplace("mask.soft",{"dx":c[0], "dy":-c[1], "dz":-c[2], "outer_radius":r*msk["nx"],"width":5})
		msk.write_image(f"{path}/mask_patch_{ci:02d}.hdf")
		#msk=msk*.75+.25
		#msk.mult(msk0)
		#msk.write_image(f"{path}/mask_patch_{ci:02d}_soft.hdf")
		
	
	for ci in range(pn):
		for eo in ["even", "odd"]:
			run(f"e2proc3d.py {oldpath}/threed_{olditer:02d}_{eo}.hdf {path}/threed_{ci*10:02d}_{eo}.hdf")
			run(f"e2proclst.py {oldpath}/ptcls_{olditer:02d}_{eo}.lst --create {path}/ptcls_{ci*10:02d}_{eo}.lst")
			p=np.loadtxt(f"{oldpath}/model_{olditer:02d}_{eo}.txt")
			np.savetxt(f"{path}/model_{ci*10-1:02d}_{eo}.txt", p)
			
		run(f"gmm_refine_iter.py {oldpath}/threed_{olditer:02d}.hdf --startres {options.startres} --initpts {oldpath}/model_{olditer:02d}.txt --mask {path}/mask_patch_{ci:02d}.hdf --masksigma --path {path} --niter {options.niter} --maskpp {path}/mask_00.hdf --startiter {ci*10+1} {etc}")
		
		for eo in ["even", "odd"]:
			run(f"e2proc3d.py {path}/threed_raw_{eo}.hdf {path}/threed_patch_{ci:02d}_raw_{eo}.hdf")
	
	run(f"gmm_merge_patch.py {path}")
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
		

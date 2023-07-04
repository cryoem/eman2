#!/usr/bin/env python
# Muyuan Chen 2023-06
from EMAN2 import *
import numpy as np
from sklearn.decomposition import PCA

def main():
	
	usage="""
	Compile movement trajectories of part of proteins from two or more alignment parameters of the same set of particles. 
	
	e2spa_trajfromrefine.py gmm_01/ptcls_*3_even.lst --baseali gmm_01/ptcls_00_even.lst 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--baseali", type=str,help="first particle alignment file", default=None)
	parser.add_argument("--nframe", type=int,help="number of frames in the trajectory", default=5)
	parser.add_argument("--nbasis", type=int,help="number of pca basis. default is 2", default=2)
	parser.add_argument("--res", type=float,help="Filter final maps to resolution. default is 5", default=5)
	parser.add_argument("--threads", type=int,help="threads", default=32)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	
	lst0=load_lst_params(options.baseali)

	params=np.zeros((len(lst0),0))
	for fm1 in args:
		lst1=load_lst_params(fm1)
		print(f"Loading {len(lst1)} particles from {fm1}...")
		dxfs=[]
		for i in range(len(lst0)):
			xf0=lst0[i]["xform.projection"]
			xf1=lst1[i]["xform.projection"]    
			dxf=xf0.inverse()*xf1
			dxf=dxf.get_params("xyz")
			dxfs.append([dxf["xtilt"],dxf["ytilt"],dxf["ztilt"],dxf["tx"],dxf["ty"],dxf["tz"]])

		dxfs=np.array(dxfs)
		params=np.hstack([params, dxfs])
		
	print("Total parameters: ", params.shape)
	options.path=os.path.dirname(options.baseali)
	
	pmfile=os.path.join(options.path,"refinestats.txt")
	np.savetxt(pmfile, params)
	print("Refine stats saved to {}")

	print("Generating eigen-trajectory...")
	
	mean=np.mean(params,0)
	std=np.std(params, 0)
	p=(params-mean)/std

	neig=max(2,options.nbasis)
	pca=PCA(neig)
	pca.fit(p)
	pfit=pca.transform(p)
	
	egfile=os.path.join(options.path,"eigval.txt")
	np.savetxt(egfile, pfit)
	print("Eigen values saved to {}".format(egfile))
	
	
	print("Making averages along trajectories...")
	ncls=options.nframe
	
	for ie in range(options.nbasis): 
		srt=np.argsort(np.argsort(pfit[:,ie]))
		srt=srt/np.max(srt+1)*ncls
		srt=np.floor(srt).astype(int)

		lout=[]
		for i,l in enumerate(lst0):
			q=l.copy()
			q["class"]=srt[i]
			lout.append(q)
		save_lst_params(lout, f"{options.path}/ptcls_mov_{ie:02d}.lst") 
		
		tfile=f"{options.path}/threed_mov_{ie:02d}.hdf"
		if os.path.isfile(tfile): os.remove(tfile)
		for i in range(ncls):
			threedtmp=f"{options.path}/threed_mov_{ie:02d}_{i:02d}.hdf"
			run(f"e2spa_make3d.py --input {options.path}/ptcls_mov_{ie:02d}.lst --output {threedtmp} --parallel thread:{options.threads} --keep 1 --sym c1 --clsid {i}")
			
			e=EMData(threedtmp)
			e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./options.res})
			e.process_inplace("normalize.edgemean")
			
			e.write_image(tfile, -1)
			if os.path.isfile(threedtmp): os.remove(threedtmp)

		print("trajectories written to {}".format(tfile))
	E2end(logid)
		
	
if __name__ == '__main__':
	main()
	

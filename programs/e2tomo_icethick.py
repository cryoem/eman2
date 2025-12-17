#!/usr/bin/env python
# Muyuan Chen 2025-10
from EMAN2 import *
import numpy as np

def main():
	
	usage="""
	Measure ice thickness of tomograms. Check result in e2tomo_eval.py. 
	e2tomo_icethick.py  
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ntile", type=int, help="number of tiles along one axis. default is 5", default=5)
	parser.add_argument("--tilesize", type=int, help="tile size. default is 256", default=256)
	parser.add_argument("--skipexist", action="store_true",help="Skip existing tomograms.", default=False)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	ntile=options.ntile
	tsz=options.tilesize
	keep=.9
	
	for tm in sorted(os.listdir("tomograms")):
		fname="tomograms/"+tm
		if options.skipexist:
			info=js_open_dict(info_name(fname))
			if info.has_key("ice_thick"):
				continue
			
		try:
			tomo=EMData(fname)
		except:
			print(f"cannot load {fname}")
			continue
	
		tomo.process_inplace("normalize.edgemean")
		apix=tomo["apix_x"]
		tomo=tomo.numpy().copy()
		
		ind=np.indices((ntile, ntile)).T.reshape(-1,2)
		shp=np.array(tomo.shape[1:])-tsz
		ind=ind*shp/(ntile-1)
		ind=ind.astype(int)
		
		stds=[]
		for ii in ind:
			t=tomo[:, ii[0]:ii[0]+tsz, ii[1]:ii[1]+tsz]
			std=np.std(t, axis=(1,2))
			std-=np.min(std)
			std/=np.max(std)
			stds.append(std)

		stds=np.array(stds)
		ss=np.std(stds, axis=1)
		ss=np.argsort(-ss)[:int(len(stds)*keep)]
		std_sel=stds[ss]
		thick=np.mean(np.sum(std_sel>.5, axis=1))*apix/10.
		print(base_name(tm), f"thickness: {thick:.1f}")
		
		
		info=js_open_dict(info_name(fname))
		info["ice_thick"]=float(thick)
		info.close()


	E2end(logid)


if __name__ == '__main__':
	main()
	

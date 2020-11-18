#!/usr/bin/env python
# Muyuan Chen 2017-10
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
from sklearn.decomposition import PCA
def main():
	
	usage="""prog --path <spt_xx folder> --iter <iteration number> [options]
	This program reads a spt refinement run with the --refine option, looks for the difference between the current alignment and the initial one, and builds a motion trajectory from it. 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--iter", type=int,help="iteration number", default=2)
	parser.add_argument("--nframe", type=int,help="number of frames in the trajectory", default=5)
	parser.add_argument("--parminit", type=str,help="compare to another parm file instead of xform in header", default=None)
	parser.add_argument("--replace3d", type=str,help="replace the 3D map used for trajectory", default=None)
	parser.add_argument("--mask", type=str,help="specify mask to use. will use mask from refinement by default", default=None)
	parser.add_argument("--nstd", type=float,help="build trajectories from -n x std to n x std of eigenvalues. default is 2", default=2)
	parser.add_argument("--nptcl", type=int,help="number of particle per average. default is 500", default=500)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#js=js_open_dict(os.path.join(options.path, "0_spt_params.json"))
	#if js.has_key("refine") and js["refine"]:
		#print("Reading from {}".format(options.path))
	#else:
		#print("This program only works on SPT projects using the --refine option")
		#return
	
	#js.close()
	
	
	js0=js_open_dict(os.path.join(options.path, "particle_parms_{:02d}.json".format(options.iter)))
	dic=dict(js0)
	js0.close()
	
	if options.parminit:
		dic00=dict(js_open_dict(options.parminit)).copy()
		
	params=[]
	pts=[]
	keys=sorted(dic.keys())
	for i, ky in enumerate(keys):
		pm0=dic[ky]
		xf0=Transform(pm0['xform.align3d'])
		src, ii =eval(ky)
		
		if options.parminit:
			xf1=Transform(dic00[ky]["xform.align3d"])
		else:
			e=EMData(src, ii, True)
			xf1=e["xform.align3d"]

		dx=np.linalg.norm(xf0.get_trans()-xf1.get_trans())

		dt=xf0*xf1.inverse()
		rot=dt.get_params("spin")["omega"]
		s=pm0["score"]
		params.append([dx, rot, s])
		
		dt=dt.get_params("xyz")
		pts.append([dt["xtilt"], dt["ytilt"], dt["ztilt"], dt["tx"], dt["ty"], dt["tz"]])
	
	pts=np.array(pts)
	pts[:,:3]=np.sin(pts[:,:3]/180*np.pi)
	
	params=np.array(params)
	pmfile=os.path.join(options.path,"refinestats_{:02d}.txt".format(options.iter))
	np.savetxt(pmfile, params)
	print("Refine stats saved to {}. The three columns are: ".format(pmfile))
	print("  Translation - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,0]), np.std(params[:,0])))
	print("  Rotation    - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,1]), np.std(params[:,1])))
	print("  Score       - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,2]), np.std(params[:,2])))
	

	print("Generating eigen-trajectory...")
	mean=np.mean(pts,0)
	std=np.std(pts, 0)
	p=(pts-mean)/std

	neig=2
	pca=PCA(neig)
	pfit=pca.fit_transform(p)
	egfile=os.path.join(options.path,"eigval_{:02d}.txt".format(options.iter))
	np.savetxt(egfile, pfit)
	print("Eigen values saved to {}".format(egfile))
	
	if options.replace3d:
		ref=EMData(options.replace3d)
	else:
		ref=EMData(os.path.join(options.path,"threed_{:02d}.hdf".format(options.iter)))
		
	if options.mask:
		msk=EMData(options.mask)
	else:
		msk=EMData(os.path.join(options.path,"mask.hdf"))
	#n=options.nframe
	for ie in range(neig): 
		v= pca.components_[ie]*std+mean
		v[:3]=np.arcsin(v[:3])*180/np.pi
		print("Trajectory {}:".format(ie))
		print("  Rotation    - x: {:.2f}  y: {:.2f}  z: {:.2f}".format(v[0],v[1], v[2]))
		print("  Translation - x: {:.2f}  y: {:.2f}  z: {:.2f}".format(v[3],v[4], v[5]))
		
			
		pv=pfit[:,ie]
		psrt=np.sort(pv[abs(pv-np.mean(pv))<np.std(pv)*options.nstd])
		rg=np.arange(options.nframe)*(len(psrt)-1)/(options.nframe-1)
		rg=psrt[rg.astype(int)]
		tfile=os.path.join(options.path,"traj_it{:02d}eg{:02d}.hdf".format(options.iter, ie))
		if os.path.isfile(tfile):
			os.remove(tfile)
		
		for i,r in enumerate(rg):
			e=ref.copy()
			t=v*r
			xf=Transform({"type":"xyz","xtilt":t[0],"ytilt":t[1],"ztilt":t[2],
					"tx":t[3],"ty":t[4],"tz":t[5]})
			e.transform(xf)
			e.write_image(tfile, i)
    
		print("Trajectory file written to {}".format(tfile))
	
	print("Making averages along trajectories...")
	jstmp="{}/particle_parms_99.json".format(options.path)
	for ie in range(neig): 
		tfile=os.path.join(options.path,"class_it{:02d}eg{:02d}.hdf".format(options.iter, ie))
		if os.path.isfile(tfile): os.remove(tfile)
		
		pv=pfit[:,ie].copy()
		pv=(pv-np.mean(pv))/np.std(pv)
		rg=np.arange(options.nframe)/options.nframe
		rg-=np.mean(rg)
		rg=rg/np.max(rg)*options.nstd
		print(rg)
		for ii, r in enumerate(rg):
			d=abs(pv-r)
			idx=np.argsort(d)[:options.nptcl]
			dcout={keys[i]:dic00[keys[i]] for i in idx}
			if os.path.isfile(jstmp): os.remove(jstmp)
			jout=js_open_dict(jstmp)
			jout.update(dcout)
			jout.close()
			
			run("e2spt_average.py --path {} --iter 99 --skippostp --threads 4 --keep 1".format(options.path))
			e=EMData("{}/threed_99.hdf".format(options.path))
			e.process_inplace("filter.matchto",{"to":ref})
			e.mult(msk)
			e.write_image(tfile, -1)
			

	E2end(logid)
		
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	
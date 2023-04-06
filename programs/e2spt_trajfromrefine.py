#!/usr/bin/env python
# Muyuan Chen 2019-05
from EMAN2 import *
import numpy as np
from sklearn.decomposition import PCA,FastICA

def main():
	
	usage="""
	Compile movement trajectories of part of proteins from two alignment parameters of the same set of particles. For example, starting from the global refinement of particles in spt_00, run another iterative local refinement in spt_01 focusing on one part of the structure using e2spt_refine_new.py with --localrefine. Then to compare the difference between angle assignment of the particles in two refinement, run
	
	e2spt_trajfromrefine.py --ali3dold spt_00/aliptcls3d_03.lst --ali3dnew spt_01/aliptcls3d_03.lst --ali2d spt_00/aliptcls2d_03.lst --path spt_00
	
	This will generate 3D movies showing the movement trajectories between spt_00 and spt_01. --ali2d defines the reference frame of the movement, and --path specifies the directory to write output.	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ali3dold", type=str,help="first particle alignment file", default=None)
	parser.add_argument("--ali3dnew", type=str,help="second particle alignment file", default=None)
	parser.add_argument("--path", type=str,help="path to write output", default=None)
	parser.add_argument("--ali2d", type=str,help="2d particle alignment file for reconstruction.", default=None)
	parser.add_argument("--nframe", type=int,help="number of frames in the trajectory", default=5)
	parser.add_argument("--maxshift", type=float,help="ignore particles with drift/rotation (pixel/degree) larger than this. default 7", default=7)
	parser.add_argument("--nstd", type=float,help="build trajectories from -n x std to n x std of eigenvalues. default is 2", default=2)
	parser.add_argument("--nbasis", type=int,help="number of pca basis. default is 2", default=2)
	parser.add_argument("--nptcl", type=int,help="number of particle per average. default is 500", default=500)
	parser.add_argument("--threads", type=int,help="threads", default=12)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	params0=load_lst_params(options.ali3dold)
	params1=load_lst_params(options.ali3dnew)
	e=EMData(options.ali3dnew)
	outsz=e["nx"]
	e=EMData(options.ali2d)
	pad=e["nx"]
	
	params=[]
	pts=[]
	for i in range(len(params0)):
		pm0=params0[i]
		pm1=params1[i]
		xf0=Transform(pm0['xform.align3d'])
		xf1=Transform(pm1['xform.align3d'])
		
		dx=np.linalg.norm(xf0.get_trans()-xf1.get_trans())

		dt=xf0*xf1.inverse()
		rot=dt.get_params("spin")["omega"]
		if "score" in pm0:
			s=pm0["score"]
		else:
			s=-1
		params.append([dx, rot, s])
		
		dt=dt.get_params("xyz")
		pts.append([dt["xtilt"], dt["ytilt"], dt["ztilt"], dt["tx"], dt["ty"], dt["tz"]])
	
	pts=np.array(pts)
	pts[:,:3]=np.sin(pts[:,:3]/180*np.pi)
	
	params=np.array(params)
	goodi=np.linalg.norm(params[:,:2], axis=1)<options.maxshift
	#print(goodi)
	
	pmfile=os.path.join(options.path,"refinestats.txt")
	np.savetxt(pmfile, params[goodi])
	print("Refine stats saved to {}. The three columns are: ".format(pmfile))
	print("  Translation - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,0]), np.std(params[:,0])))
	print("  Rotation    - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,1]), np.std(params[:,1])))
	print("  Score       - mean: {:.2f}  std: {:.2f}".format(np.mean(params[:,2]), np.std(params[:,2])))
	

	print("Generating eigen-trajectory...")
	mean=np.mean(pts,0)
	std=np.std(pts, 0)
	p=(pts-mean)/std

	neig=max(2,options.nbasis)
	pca=PCA(neig)
	pca.fit(p[goodi])
	pfit=pca.transform(p)
	pfit=(pfit-np.mean(pfit[goodi], axis=0))/np.std(pfit[goodi], axis=0)
	pfit[~goodi]=1e5
	egfile=os.path.join(options.path,"eigval.txt")
	np.savetxt(egfile, pfit[goodi])
	print("Eigen values saved to {}".format(egfile))
	
	
	print("Making averages along trajectories...")
	info3d=load_lst_params("{}/particle_info_3d.lst".format(options.path))
	jstmp="{}/tmp_ptcl_list.txt".format(options.path)
	threedtmp="{}/threed_tmp.hdf".format(options.path)
	for ie in range(options.nbasis): 
		print("Eigenvector {}".format(ie))
		tfile=os.path.join(options.path,"threed_eig_{:02d}.hdf".format(ie))
		if os.path.isfile(tfile): os.remove(tfile)
		
		pv=pfit[:,ie].copy()
		#pv=(pv-np.mean(pv))/np.std(pv)
		rg=np.arange(options.nframe)/options.nframe
		rg-=np.mean(rg)
		rg=rg/np.max(rg)*options.nstd
		print(rg)
		for ii, r in enumerate(rg):
			d=abs(pv-r)
			idx=np.argsort(d)[:options.nptcl]
			idx2d=[info3d[i]["idx2d"] for i in idx]
			idx2d=sum(idx2d,[])
			np.savetxt(jstmp, idx2d)
			
			run(f"e2spa_make3d.py --input {options.ali2d} --output {threedtmp} --keep 1 --parallel thread:{options.threads} --outsize {outsz} --pad {pad} --listsel {jstmp}")
			e=EMData(threedtmp)
			#e.process_inplace("filter.matchto",{"to":ref})
			#e.mult(msk)
			e.write_image(tfile, -1)
			if os.path.isfile(threedtmp): os.remove(threedtmp)

		print("trajectories written to {}".format(tfile))
	E2end(logid)
		
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	

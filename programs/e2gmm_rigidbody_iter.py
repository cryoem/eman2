#!/usr/bin/env python
# Muyuan Chen 2024-02
from EMAN2 import *
import numpy as np
from sklearn.cluster import KMeans

	
def main():
	
	usage="""
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path for refinement. default is the next gmm_xx", default=None)
	parser.add_argument("--niter", type=int,help="number of iteration. default is 5", default=5)
	parser.add_argument("--npatch", type=int,help="number of patches. default is 8", default=8)
	parser.add_argument("--maxres", type=float,help="max resolution to consider in refinement. i.e. the alignment will not use information beyond this even if FSC goes further.", default=4.)
	parser.add_argument("--minres", type=float,help="minimum resolution for comparision.", default=25.)
	parser.add_argument("--startres", type=float,help="start resolution.", default=-1)
	parser.add_argument("--expandsym", type=str,help="symmetry expansion. i.e. start from an input refinement with the given symmetry and perform the new refinement with c1 by making copies of each particle at symmetrical equivalent positions. ", default=None)
	parser.add_argument("--masktight", action="store_true", default=False ,help="Use tight patch mask instead of spherical ones. seems safe.")
	parser.add_argument("--masks", type=str,help="masks for refinement. multiple files separated by ','. Replace npatch if specified.", default=None)
	parser.add_argument("--applysym", type=str,help="apply symmetry to reconstruction.", default="c1")
	parser.add_argument("--startiter", type=int,help="starting iteration. default is 1", default=1)
	parser.add_argument("--l2_weight", type=float, help="weighting factor for L2. default is 1",default=1.)
	parser.add_argument("--chunksize", type=int,help="Number of particles in each e2gmm_batch process. Increase will make the alignment slightly faster, but also increases CPU memory use. ", default=25000)
	parser.add_argument("--storexf", action="store_true", default=False ,help="store the transform after refinement and start from the transform in the last iteration each time.")
	parser.add_argument("--parallel", type=str,help="parallel options for 3d reconstruction", default="thread:32")
	parser.add_argument("--initpts", type=str,help="gmm input to replace model from initial folder", default=None)

	# parser.add_argument("--recover",action="store_true", default=False ,help="continue previous crashed refinement")


	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	fname=args[0]
	oldpath=os.path.dirname(fname)
	options.thread=int(options.parallel.split(':')[1])
	
	c=[i for i in fname if i.isdigit()]
	olditer=int(''.join(c[-2:]))
	print(f"refinement path {oldpath}, iteration {olditer}.")
	
	if options.masks!=None:
		maskfile=options.masks.split(',')
		pn=options.npatch=len(maskfile)
		
	#### prepare iteration 00
	if options.path==None: 
		path=options.path=num_path_new("gmm_")
		iiter=0
		
		for ieo, eo in enumerate(["even", "odd"]):
			
			it0=max(0, olditer-1)
			
			if os.path.isfile(f"{oldpath}/ptcls_{olditer:02d}_{eo}.lst"):
				run(f"e2proclst.py {oldpath}/ptcls_{olditer:02d}_{eo}.lst --create {path}/ptcls_00_{eo}.lst")
			else:
				run(f"e2proclst.py {oldpath}/ptcls_{olditer:02d}.lst --create {path}/ptcls_00_{eo}.lst --getclass {ieo}")
				
			if options.expandsym:
				run(f"e2proclst.py {path}/ptcls_00_{eo}.lst --sym {options.expandsym}")
			
			run(f"e2proc3d.py {oldpath}/threed_{olditer:02d}_{eo}.hdf {path}/threed_{iiter:02d}_{eo}.hdf")
			run(f"e2proc3d.py {oldpath}/threed_raw_{eo}.hdf {path}/threed_{iiter:02d}_raw_{eo}.hdf")
			
			
			if options.initpts==None:
				pts=np.loadtxt(f"{oldpath}/model_{it0:02d}_{eo}.txt")
				np.savetxt(f"{path}/model_00_{eo}.txt", pts)
			else:
				run(f"e2project3d.py {path}/threed_{iiter:02d}_{eo}.hdf --outfile {path}/tmp_projection.hdf --orientgen=eman:delta=4 --parallel={options.parallel}")
				run(f"e2gmm_refine_new.py --ptclsin {path}/tmp_projection.hdf --model {options.initpts} --maxres {options.maxres} --minres {options.minres} --modelout {path}/model_00_{eo}.txt --niter 40 --trainmodel --learnrate 1e-5")
			
			
			# run(f"e2spa_make3d.py --input {path}/ptcls_{iiter:02d}_{eo}.lst --output {path}/threed_{iiter:02d}_{eo}.hdf --parallel thread:{options.thread} --keep 1 --sym c1")
			# e=EMData(f"{path}/threed_{iiter:02d}_{eo}.hdf")
			# e.write_image(f"{path}/threed_{iiter:02d}_raw_{eo}.hdf")
			
		run(f"e2proc3d.py {oldpath}/mask.hdf {path}/mask_00.hdf")
		# run(f"e2refine_postprocess.py --even {path}/threed_{iiter:02d}_even.hdf --res 5 --tophat localwiener --sym c1 --thread {options.thread} --setsf sf.txt --align --mask {path}/mask_00.hdf")
		
		
		#### make masks
		if options.masks==None:
			pn=options.npatch
			
			p00=np.loadtxt(f"{oldpath}/model_{olditer-1:02d}_even.txt")
			km=KMeans(pn,max_iter=30)
			km.fit(p00[:,:3])
			pc=km.cluster_centers_
			lb=km.labels_
		
			msk=EMData(f"{path}/mask_00.hdf")
				
			if options.masktight:
				for ci in range(pn):
					ps=p00[lb==ci,:].copy()
					print(ps.shape)
					np.savetxt(f"{path}/model_tmp.txt", ps)
					run(f"e2gmm_refine_new.py --ptclsin {oldpath}/projections_even.hdf  --model {path}/model_tmp.txt --maxres -1 --modelout {path}/model_tmp.txt --niter 0 --trainmodel --evalmodel {path}/model_projs.hdf")
					run(f"e2spa_make3d.py --input {path}/model_projs.hdf --output {path}/model_avg.hdf --thread {options.thread}")
					e=EMData(f"{path}/model_avg.hdf")
					e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
					e.process_inplace("normalize.edgemean")
					e.process_inplace("mask.auto3d.thresh",{"nshells":3,"nshellsgauss":5,"threshold1":3,"threshold2":1,"return_mask":True})
					e.write_image(f"{path}/mask_patch_{ci:02d}.hdf")
				
			else:
				for ci in range(pn):
					c=pc[ci].copy()
					r=np.linalg.norm(p00[lb==ci,:3]-c, axis=1)
					r=np.max(r)
					
					msk.to_one()
					c=(c*msk["nx"]).tolist()
					msk.process_inplace("mask.soft",{"dx":c[0], "dy":-c[1], "dz":-c[2], "outer_radius":r*msk["nx"],"width":5})
					msk.write_image(f"{path}/mask_patch_{ci:02d}.hdf")
					
		else:
			
			for ci,m in enumerate(maskfile):
				e=EMData(m)
				e.write_image(f"{path}/mask_patch_{ci:02d}.hdf")
				
	else:
		path=options.path
		
		
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_gmm_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()	

	if options.storexf:
		for eo in ["even", "odd"]:
			nptcl=EMUtil.get_image_count(f"{path}/ptcls_00_{eo}.lst")
			last_xf=np.zeros((nptcl, options.npatch*6))
			last_xf=from_numpy(last_xf)
			last_xf.write_image(f"{path}/xf_{eo}.hdf")

	#### starting refinement
	maskfile=[f"{path}/mask_patch_{ci:02d}.hdf" for ci in range(options.npatch)]
	masks=','.join(maskfile)
	if options.startres<0:
		resmult=[1.5,1.5,1.3,1.1,1.0]+[1.0]*options.niter
	else:
		resmult=np.arange(options.niter)/(options.niter-1)
		resmult=resmult[::-1]*(options.startres-options.maxres)+options.maxres
		print(resmult)
		resmult/=options.maxres
		resmult=np.append(resmult[0], resmult)
		resmult=np.append(resmult, resmult[-1])
	
	for iiter in range(options.startiter, options.niter+1):
		res=options.maxres*resmult[iiter]
		
		for eo in ["even", "odd"]:
			
			etc=""
			if options.storexf: etc+=f" --xf_file {path}/xf_{eo}.hdf"
			
			run(f"e2gmm_batch.py 'e2gmm_rigidbody.py --ptclsin {path}/ptcls_00_{eo}.lst --ptclsout {path}/ptcls_{iiter:02d}_{eo}.lst --model {path}/model_{iiter-1:02d}_{eo}.txt --mask {masks} --maxres {res} --l2_weight {options.l2_weight} --minres {options.minres} {etc}' --batch {options.chunksize} --niter 0")
		
			
			for pid in range(options.npatch):
				run(f"e2spa_make3d.py --input {path}/ptcls_{iiter:02d}_p{pid:02d}_{eo}.lst --output {path}/threed_{iiter:02d}_p{pid:02d}_{eo}.hdf --parallel {options.parallel} --keep 1 --sym c1")
		
			
			#### merge patches here
			avg=EMData(f"{path}/threed_00_raw_{eo}.hdf")
			avg.to_zero()
			wt=avg.copy()
			
			for pid,mfile in enumerate(maskfile):
				m=EMData(mfile)
				e=EMData(f"{path}/threed_{iiter:02d}_p{pid:02d}_{eo}.hdf")
				e.mult(m)
				avg.add(e)
				wt.add(m)

			#### copy background noise from the first map
			wt0=1-wt
			wt0.process_inplace("threshold.belowtozero")
			e=EMData(f"{path}/threed_00_raw_{eo}.hdf")
			e.mult(wt0)
			avg.add(e)
			wt.add(wt0)

			wt.process_inplace("math.reciprocal", {"zero_to":0})
			avg.mult(wt)
			if options.applysym!="c1":
				avg.process_inplace("xform.applysym",{"sym":options.applysym})
			avg.write_image(f"{path}/threed_{iiter:02d}_{eo}.hdf")
			avg.write_image(f"{path}/threed_raw_{eo}.hdf")
			
			
		run(f"e2refine_postprocess.py --even {path}/threed_{iiter:02d}_even.hdf --res 5 --tophat localwiener --sym c1 --thread {options.thread} --setsf sf.txt --align --mask {path}/mask_00.hdf")
			
		rr=options.maxres*resmult[iiter+1]
		for eo in ["even","odd"]:
			run(f"e2project3d.py {path}/threed_{iiter:02d}_{eo}.hdf --outfile {path}/projections_{eo}.hdf --orientgen=eman:delta=4 --parallel={options.parallel}")
			
			run(f"e2gmm_refine_new.py --ptclsin {path}/projections_{eo}.hdf --model {path}/model_{iiter-1:02d}_even.txt --maxres {rr} --modelout {path}/model_{iiter:02d}_{eo}.txt --niter 40 --trainmodel --learnrate 1e-6")
			
			
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
		

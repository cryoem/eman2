#!/usr/bin/env python
# Muyuan Chen 2022-12
from EMAN2 import *
import numpy as np
from sklearn.cluster import KMeans

def main():
	
	usage="""Refine heterogeneious domain undergoing large scale continuous motion. Start from hetergenerity analysis using e2gmm_refine_new, then convert the conformation to orientation within a given mask, and run a few rounds of orientation refinement focusing on the mask. 
	
	e2gmm_heter_refine.py gmm_00/threed_05.hdf --mask mask.hdf --maxres 3.5
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path for refinement. default is the next gmm_xx", default=None)
	parser.add_argument("--mask", type=str,help="mask that defines the region of focusing.", default=None)
	parser.add_argument("--ngauss", type=int,help="number of gauss used in initial heterogeneity analysis.", default=8000)
	parser.add_argument("--niter", type=int,help="iterations of focused alignment afterwards.", default=5)
	parser.add_argument("--batchsz", type=int,help="batch size.", default=16)
	parser.add_argument("--expandsym", type=str,help="symmetry. the program does not apply symmetry so always specify symmetry here and the final structure will be in c1", default=None)
	parser.add_argument("--maxres", type=float,help="maximum resolution for the heterogeneity analysis. This will also be the starting resolution for the following focused alignment.", default=7)
	parser.add_argument("--minres", type=float,help="min resolution for the heterogeneity analysis.", default=50)	
	parser.add_argument("--rigidbody", action="store_true", default=False ,help="rigid body movement mode. still testing")
	parser.add_argument("--remodel", action="store_true", default=False ,help="rebuild model after heterogeneity analysis for later focused refinement")

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	
	fname=args[0]
	er=EMData(fname,0, True)
	oldpath=os.path.dirname(fname)
	
	c=[i for i in fname if i.isdigit()]
	itr=int(''.join(c[-2:]))
	print(f"refinement path {oldpath}, iteration {itr}.")
	
	if options.path==None: 
		path=options.path=num_path_new("gmm_")
		
		fsc=np.loadtxt(f"{oldpath}/fsc_masked_{itr:02d}.txt")
		np.savetxt(f"{path}/fsc_masked_00.txt", fsc)
		fsc=np.loadtxt(f"{oldpath}/fsc_maskedtight_{itr:02d}.txt")
		np.savetxt(f"{path}/fsc_maskedtight_00.txt", fsc)
		run(f"e2proc3d.py {oldpath}/threed_{itr:02d}.hdf {path}/threed_00.hdf")
		
		for ieo, eo in enumerate(["even", "odd"]):
			run(f"e2proclst.py {oldpath}/ptcls_{itr:02d}_{eo}.lst --create {path}/ptcls_00_{eo}.lst")
			if options.expandsym:
				run(f"e2proclst.py {path}/ptcls_00_{eo}.lst --sym {options.expandsym}")
				
			run(f"e2proc3d.py {oldpath}/threed_{itr:02d}_{eo}.hdf {path}/threed_00_{eo}.hdf")
			
			it0=max(0, itr-1)
			pts=np.loadtxt(f"{oldpath}/model_{it0:02d}_{eo}.txt")
			np.savetxt(f"{path}/model_00_{eo}.txt", pts)
			
			run(f"e2project3d.py {path}/threed_00_{eo}.hdf --outfile {path}/projections_{eo}.hdf --orientgen=eman:delta=4 --parallel=thread:12")
				
	
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_gmm_heter_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()	
	
	it0=0
	path=options.path
	res=options.maxres
	for ieo, eo in enumerate(["even", "odd"]):
		
		pts=np.loadtxt(f"{path}/model_{it0:02d}_{eo}.txt")
		rawnpt=len(pts)
		npt=options.ngauss
		if npt>0:
			np.random.shuffle(pts)
			pts=pts[:npt]
		np.savetxt(f"{path}/model_{it0:02d}_{eo}.txt", pts)
		run(f"e2gmm_refine_new.py --ptclsin {path}/projections_{eo}.hdf --model {path}/model_{it0:02d}_{eo}.txt --maxres {res} --modelout {path}/model_{it0:02d}_{eo}.txt --niter 40 --trainmodel --learnrate 1e-6")
		
		pn=16
		km=KMeans(pn,max_iter=30)
		km.fit(pts[:,:3])
		pc=km.cluster_centers_
		
		emmask=True
		try:
			msk=EMData(options.mask)
		except:
			emmask=False
			
		if emmask:
			## read selected Gaussian from mask file
			m=msk.numpy().copy()
			p=pts[:,:3].copy()
			p=p[:,::-1]
			p[:,:2]*=-1
			p=(p+.5)*msk["nx"]

			o=np.round(p).astype(int)
			v=m[o[:,0], o[:,1], o[:,2]]
			imsk=v.astype(float)
		else:
			i=np.loadtxt(options.mask).astype(int).flatten()
			imsk=np.zeros(len(pts), dtype=float)
			imsk[i]=1
		
		if options.rigidbody:
			mode="--rigidbody"
		else:
			pm=pts[imsk>.1]
			pn=16
			km=KMeans(pn,max_iter=30)
			km.fit(pm[:,:3])
			pc2=km.cluster_centers_
			
			pcx=np.vstack([pc, pc2])
			pp=np.hstack([pcx, np.zeros((len(pcx),1))+np.mean(pts[:,3]), np.zeros((len(pcx),1))+np.mean(pts[:,4])])
			np.savetxt(f"{path}/model_{it0:02d}_{eo}_anchor.txt", pp)
			mode=f"--anchor {path}/model_{it0:02d}_{eo}_anchor.txt"
		
		
		lst=load_lst_params(f"{path}/ptcls_{it0:02d}_{eo}.lst")
		np.random.shuffle(lst)
		save_lst_params(lst[:20000], f"{path}/ptcls_{it0:02d}_{eo}_sample.lst")

		etc=""
		etc+=f" --selgauss {options.mask}"
		etc+=f" --batchsz {options.batchsz}"
		etc3=""# --setsf sf_lp.txt"
		
		## pretrain from lower res to ensure convergence
		res0=max(res,7)
		run(f"e2gmm_refine_new.py --model {path}/model_{it0:02d}_{eo}.txt {mode} --conv --midout {path}/mid_00_{eo}.txt --maxres {res0} --minres {options.minres} --learnrate 1e-5 --niter 20 --ptclsin {path}/ptcls_{it0:02d}_{eo}_sample.lst --heter --encoderout {path}/enc_{eo}.h5 --decoderout {path}/dec_{eo}.h5 --pas 100 {etc} --maxgradres {res}")
		# run(f"e2gmm_refine_new.py --model {path}/model_{it0:02d}_{eo}.txt --rigidbody --midout {path}/mid_00_{eo}.txt --maxres {res0} --minres {options.minres} --learnrate 1e-5 --niter 20 --ptclsin {path}/ptcls_{it0:02d}_{eo}_sample.lst --heter --encoderout {path}/enc_{eo}.h5 --decoderout {path}/dec_{eo}.h5 --pas 100 {etc} --maxgradres {res}")
		
		run(f"e2gmm_eval.py --pts {path}/mid_00_{eo}.txt --pcaout {path}/mid_pca.txt --ptclsin {path}/ptcls_{it0:02d}_{eo}_sample.lst --ptclsout {path}/ptcls_{eo}_cls_00.lst --mode regress --ncls 4 --nptcl 8000 --axis 0")
		
		run(f"e2proc3d.py {path}/ptcls_{eo}_cls_00.hdf {path}/ptcls_{eo}_cls_00.hdf --process filter.lowpass.gauss:cutoff_freq={1./res} --process normalize.edgemean")
		
		run(f'e2gmm_batch.py "e2gmm_refine_new.py --model {path}/model_{it0:02d}_{eo}.txt --midout {path}/midall_00_{eo}.txt {mode} --maxres {res} --minres {options.minres} --learnrate 1e-5 --niter 10 --ptclsin {path}/ptcls_{it0:02d}_{eo}.lst --heter --encoderout {path}/enc_{eo}.h5 --decoderout {path}/dec_{eo}.h5 --pas 100 {etc}" --load --niter 1 --batch 20000')
		
		run(f"e2gmm_eval.py --pts {path}/midall_00_{eo}.txt --pcaout {path}/mid_pca.txt --ptclsin {path}/ptcls_{it0:02d}_{eo}.lst --ptclsout {path}/ptcls_{eo}_cls_00.lst --mode regress --ncls 4 --nptcl 5000 --axis 0")
		
		run(f"e2proc3d.py {path}/ptcls_{eo}_cls_00.hdf {path}/ptcls_{eo}_cls_00.hdf --process filter.lowpass.gauss:cutoff_freq={1./res} --process normalize.edgemean")
		
	run(f"e2gmm_heter_ali.py --path {path} --mask {options.mask}")
	run(f"e2refine_postprocess.py --even {path}/threed_01_even.hdf --res {res} --tophat localwiener --thread 32 --setsf sf.txt --align")
	
	if options.remodel:
		run(f"e2segment3d.py {path}/threed_01.hdf --pdb {path}/model_01.pdb --process=segment.kmeans:nseg={rawnpt}:thr=3")
		inp=f"{path}/model_01.pdb"
	else:
		inp=f"{path}/model_00.txt"
	
	run(f"e2gmm_refine_iter.py {path}/threed_01.hdf --startres {res} --initpts {inp} --mask {options.mask} --masksigma --path {path} --startiter 2 --niter {options.niter}")
	
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	



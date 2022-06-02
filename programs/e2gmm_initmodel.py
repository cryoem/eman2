#!/usr/bin/env python
# Muyuan Chen 2022-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default="gmm_00")
	parser.add_argument("--setsf", type=str,help="set structure factor", default=None)
	parser.add_argument("--batchsize", type=int,help="batch size", default=256)
	parser.add_argument("--npt", type=int,help="number of gaussian in final model", default=256)
	parser.add_argument("--minamp", type=float,help="remove gaussian with low amplitude", default=.8)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	path=options.path
	if options.setsf:
		sf=f" --setsf {options.setsf}"
	else:
		sf=""
	npt=8
	
	run(f"e2spa_sgd.py {args[0]} --niter 0 --path {path} --batch {options.batchsize}")

	run(f"e2proc3d.py {path}/output.hdf {path}/output.hdf --process filter.lowpass.gauss:cutoff_abs=.02 --process normalize.edgemean")

	run(f"e2segment3d.py {path}/output.hdf --pdb {path}/model_init.pdb --process=segment.kmeans:ampweight=1:nseg={npt}:verbose=1:thr=3")
	
	run(f"e2gmm_refine_new.py --projs {path}/ptcls.lst --npt 8 --model {path}/model_init.pdb --modelout {path}/model_gmm.txt --niter 40 --evalmodel {path}/model_projs.hdf --maxres 40")
	
	run(f"e2spa_make3d.py --input {path}/model_projs.hdf --output {path}/threed_model.hdf --threads 32 {sf}")

	res=40
	for i in range(40):
		run(f"e2proclst.py {path}/ptcls_input.lst --create {path}/tmp.lst")
		run(f"e2proclst.py {path}/tmp.lst --shuffle")
		run(f"e2proclst.py {path}/tmp.lst --create {path}/ptcls_sample.lst --range {options.batchsize}")
		
		run(f"e2gmm_refine_new.py --ptclsin {path}/ptcls_sample.lst --model {path}/model_gmm.txt --align --maxres {res} --batchsz 128 --ptclsout {path}/ptcls.lst --niter 40 --fromscratch")
		run(f"e2spa_make3d.py --input {path}/ptcls.lst --output {path}/threed_ptcl.hdf --threads 32 {sf}")
		
		p=np.loadtxt(f"{path}/model_gmm.txt")
		#np.random.shuffle(p)
		p=p[len(p)//8:]
		p[:,3]/=np.max(p[:,3])
		p=p[p[:,3]>options.minamp,:]
		p[:,:3]-=np.mean(p[:,:3], axis=0)
		#p[:,:3]+=np.random.randn(len(p), 3)*1e-2
		np.savetxt(f"{path}/model_gmm.txt", p)
		
		run(f"e2gmm_refine_new.py --projs {path}/ptcls.lst --npt {npt} --model {path}/model_gmm.txt --modelout {path}/model_gmm.txt --niter 20 --evalmodel {path}/model_projs.hdf --maxres {res}")
		run(f"e2spa_make3d.py --input {path}/model_projs.hdf --output {path}/threed_model.hdf --threads 32 {sf}")
		
		run(f"e2proc3d.py {path}/threed_ptcl.hdf {path}/output_all.hdf --append")
		run(f"e2proc3d.py {path}/threed_model.hdf {path}/threed_model_all.hdf --append")
		
		
		if npt<options.npt and (i+1)%5==0:
			npt*=2
		if i>20:
			res=20
			
	E2end(logid)
	
def run(cmd):
	print(cmd)
	os.system(cmd)
	
if __name__ == '__main__':
	main()
	

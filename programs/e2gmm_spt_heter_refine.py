#!/usr/bin/env python
# Muyuan Chen 2025-03
import numpy as np
from EMAN2 import *

def main():
	
	usage="""e2gmm_spt_heter_refine.py spt_xx/threed_xx.hdf
	It should find all relavent files to do the heterogeneity analysis. 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="output path", default=None)
	#parser.add_argument("--ptcls", type=str,help="particles input", default=None)
	parser.add_argument("--model", type=str,help="model input", default=None)
	parser.add_argument("--mask", type=str,help="masks for refinement.", default=None)
	parser.add_argument("--niter", type=int, help="number of iteration",default=20)
	parser.add_argument("--ntilt", type=int, help="number of 2D tilt per 3D particle. If set, only use the information from the central N tilts.",default=-1)
	
	parser.add_argument("--randsym", type=str,help="assign particle to random symmetry unit", default=None)
	parser.add_argument("--n_anchor", type=int,help="number of anchor points. default 32", default=32)
	# parser.add_argument("--xfin_starti", type=int,help="starting index for tranform input", default=0)

	parser.add_argument("--maxres", type=float, help="resolution",default=20.)
	parser.add_argument("--minres", type=float, help="resolution",default=200.)
	parser.add_argument("--learnrate", type=float, help="learning rate",default=1e-5)
	# parser.add_argument("--angle_range", type=str,help="search angle range", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. use 3 digit 0/1 input. default is 111", default="111")
	
	(options, args) = parser.parse_args()	
	logid=E2init(sys.argv)
	
	fname=args[0]
	oldpath=os.path.dirname(fname)
	olditer=int(fname[:-4].split('_')[-1])
	p2d=f"{oldpath}/aliptcls2d_{olditer:02d}.lst"
	print(f"Input path {oldpath}, iter {olditer}, alignment file {p2d}")
	
	
	#### parse particles
	alipm=load_lst_params(p2d)
	pids=np.array([a["ptcl3d_id"] for a in alipm])
	uid=np.unique(pids)
	p3did=[np.where(pids==u)[0] for u in uid]

	##   keep track of 3D particles from the 2D ones
	idx=np.unique(pids)
	idx3d=[]
	for i in idx:
		idx3d.append(np.where(pids==i)[0])

	ln=np.mean([len(i) for i in idx3d])
	print("{} 3D particles, each contain {:.1f} 2D particles".format(len(idx3d), ln))
	
	if options.ntilt>0:
		print(f"Keeping {options.ntilt} 2D tilt per 3D particle")
		idx3d_keep=[]
		for idx in idx3d:
			tid=np.array([alipm[i]["tilt_id"] for i in idx])
			tid=np.argsort(abs(tid-np.mean(tid)))[:options.ntilt]
			idx3d_keep.append(idx[np.sort(tid)])
			
		idx3d=idx3d_keep
		
	if options.path:
		path=options.path
		print(f"writing in existing directory {path}")
	else:
		path=options.path=num_path_new("gmm_")
		print(f"writing in directory {path}")
	
	lout=[]
	for idx in idx3d:
		if options.randsym==None:
			lout.extend([alipm[i] for i in idx])
		else:
			apm=[alipm[i] for i in idx]
			for a in apm:
				x=a["xform.projection"]
				si=np.random.randint(x.get_nsym(options.randsym))
				xt=x.get_sym(options.randsym, si)
				a["xform.projection"]=xt
			lout.extend(apm)
			
	save_lst_params(lout, f"{path}/aliptcls2d_00.lst")
	
	rr=np.random.randint(0,len(idx3d),500)
	itest=[idx3d[i] for i in rr]
	lout=[]
	for idx in itest:
		lout.extend([alipm[i] for i in idx])
	save_lst_params(lout, f"{path}/aliptcls2d_test.lst")
	
	if options.model:
		m=np.loadtxt(options.model)
		np.savetxt(f"{path}/model_00.txt", m)
	else:
		print("need model input for now")
		return
	
	if options.mask==None:
		options.mask=f"{oldpath}/mask.hdf"
		
	e=EMData(fname,0,True)
	nx=e["nx"]

	
	run(f"e2gmm_spt_heterg.py --ptclsin {path}/aliptcls2d_test.lst --model {path}/model_00.txt --mask {options.mask} --clip {nx} --midout {path}/midout_test.txt --encoderout {path}/encoder.h5 --decoderout {path}/decoder.h5 --anchor {path}/anchor_00.txt --maxres {options.maxres} --minres {options.minres} --learnrate {options.learnrate} --niter {options.niter} --pas {options.pas}")
	
	nn=100
	run(f"e2gmm_eval.py --pts {path}/midout_test.txt --pcaout {path}/mid_pca.txt --ptclsin {path}/aliptcls2d_test.lst --ptclsout {path}/ptcls_cls_test.lst --mode regress --ncls 4 --nptcl {nn} --axis 0 --spt --outsize {nx}")
	
	run(f"e2proc3d.py {path}/ptcls_cls_test.hdf {path}/ptcls_cls_test.hdf --process filter.lowpass.gauss:cutoff_freq={1./options.maxres} --process normalize.edgemean")
	
	
	run(f'e2gmm_batch.py "e2gmm_spt_heterg.py --ptclsin {path}/aliptcls2d_00.lst --model {path}/model_00.txt --mask {options.mask} --clip {nx} --midout {path}/midout_all_00.txt --encoderout {path}/encoder.h5 --decoderout {path}/decoder.h5 --anchor {path}/anchor_00.txt --maxres {options.maxres} --minres {options.minres} --learnrate {options.learnrate} --niter {options.niter} --pas {options.pas}" --ptcl3d --batch 500 --niter 1 --load')
	
	nn=500
	run(f"e2gmm_eval.py --pts {path}/midout_all_00.txt --pcaout {path}/mid_pca.txt --ptclsin {path}/aliptcls2d_00.lst --ptclsout {path}/ptcls_cls_00.lst --mode regress --ncls 4 --nptcl {nn} --axis 0 --spt --outsize {nx}")
	
	run(f"e2proc3d.py {path}/ptcls_cls_test.hdf {path}/ptcls_cls_test.hdf --process filter.lowpass.gauss:cutoff_freq={1./options.maxres} --process normalize.edgemean")
	
	E2end(logid)
	return


if __name__ == '__main__':
	main()
	

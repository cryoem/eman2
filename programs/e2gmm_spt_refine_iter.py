#!/usr/bin/env python
# Muyuan Chen 2026-06
from EMAN2 import *
import numpy as np
	
def main():

	usage="""
	Iterative orientation refinement using GMM as references for SPT. Require GPU and CUDA.
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path for refinement. default is the next gmm_xx", default=None)
	parser.add_argument("--sym", type=str,help="symmetry to apply to the map", default="c1")
	parser.add_argument("--initpts", type=str,help="initial Gaussian coordinates to initialize the GMMs.", default=None)
	parser.add_argument("--niter", type=int,help="number of iteration. default is 3.", default=3)
	parser.add_argument("--startres", type=float,help="starting resolution", default=10)
	parser.add_argument("--minres", type=float,help="min resolution", default=50)
	parser.add_argument("--mask", type=str,help="mask for focused refinement.", default=None)
	parser.add_argument("--maskpp", type=str,help="mask for the final structure.", default=None)
	parser.add_argument("--tophat", type=str,help="tophat filter for post process.", default="localwiener")
	parser.add_argument("--align_mlp", action="store_true", default=False ,help="DNN based alignment")

	parser.add_argument("--startiter", type=int,help="starting iteration number.", default=1)
	parser.add_argument("--parallel", type=str,help="for e2spa_make3d.", default="thread:64")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--chunksize", type=int,help="Number of particles in each e2gmm_batch process.", default=3000)
	parser.add_argument("--ntilt", type=int,help="Number of tilt to use per subtomogram. default is 25.", default=25)
	parser.add_argument("--batchsize", type=int,help="Number of 3d particle per batch. default is 2", default=2)
	parser.add_argument("--mergelowres", type=float,help="merge low resolution information in the initial reference. avoid two half set converging to different conformations.", default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	fname=args[0]
	er=EMData(fname,0, True)
	oldpath=os.path.dirname(fname)
	
	c=[i for i in fname if i.isdigit()]
	itr=int(''.join(c[-2:]))
	print(f"refinement path {oldpath}, iteration {itr}.")
	res=options.startres
	if options.mask==None:
		options.mask=f"{oldpath}/mask.hdf"
		globalrefine=True
	else:
		globalrefine=False

	if options.path==None: 
		path=options.path=num_path_new("gmm_")
		for tag in ["masked", "maskedtight", "unmasked"]:
			fsc=np.loadtxt(f"{oldpath}/fsc_{tag}_{itr:02d}.txt")
			np.savetxt(f"{path}/fsc_{tag}_00.txt", fsc)

		for ieo, eo in enumerate(["even", "odd"]):
			if os.path.isfile(f"{oldpath}/aliptcls2d_{itr:02d}_{eo}.lst"):
				run(f"e2proclst.py {oldpath}/aliptcls2d_{itr:02d}_{eo}.lst --create {path}/aliptcls2d_00_{eo}.lst")
			else:
				run(f"e2proclst.py {oldpath}/aliptcls2d_{itr:02d}.lst --create {path}/aliptcls2d_00_{eo}.lst --getclass {ieo}")
				
			run(f"e2proc3d.py {oldpath}/threed_{itr:02d}_{eo}.hdf {path}/threed_00_{eo}.hdf --process filter.lowpass.gauss:cutoff_freq={1./options.startres:.3f} --process normalize.edgemean")
	
		if options.mergelowres>0:
			even=EMData(f"{path}/threed_00_even.hdf",0)
			odd=EMData(f"{path}/threed_00_odd.hdf",0)
			c=(even+odd)*0.5
			rs=1./options.mergelowres
			even.process_inplace("filter.highpass.tophat", {"cutoff_freq":rs})
			odd.process_inplace("filter.highpass.tophat", {"cutoff_freq":rs})
			c.process_inplace("filter.lowpass.tophat", {"cutoff_freq":rs})
			even.add(c)
			odd.add(c)
			even.write_image(f"{path}/threed_00_even.hdf")
			odd.write_image(f"{path}/threed_00_odd.hdf")
	else:
		print(f"writing in existing path {options.path}")
		path=options.path
		
	ext=options.initpts[-3:]
	if os.path.isfile(f"{options.initpts[:-4]}_even.{ext}"):
		print("using even/odd pdb")
		for ieo, eo in enumerate(["even", "odd"]):
			os.system(f"cp {options.initpts[:-4]}_{eo}.{ext} {path}/gmm_init_{eo}.{ext}")
	else:
		for ieo, eo in enumerate(["even", "odd"]):
			os.system("cp {} {}/gmm_init_{}.{}".format(options.initpts, path,eo,ext))
		
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_gmm_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()	

	etcpp=""
	if options.maskpp:
		etcpp+=f" --mask {options.maskpp}"
		
	for itr in range(options.startiter, options.startiter+options.niter):
		it0=itr-1
		
		for ieo, eo in enumerate(["even", "odd"]):
			ref_input=f"{path}/threed_{it0:02d}_{eo}.hdf"
			model_input=f"{path}/gmm_init_{eo}.{ext}"
			run(f"e2gmm_refine_jax.py --ptclsin {ref_input} --model {model_input} --maxres {res} --modelout {path}/model_{it0:02d}_{eo}.txt --trainmodel")
	
			e=EMData(f"{path}/threed_{it0:02d}_{eo}.hdf",0,True)
			p=EMData(f"{path}/aliptcls2d_{it0:02d}_{eo}.lst",0,True)
			etcali=""; etcm3d=""; etcbatch=""
			if p["nx"]>e["nx"]:
				clip=e["nx"]
				etcali+=f" --clip {clip}"
				etcm3d+=f" --outsize {clip}"

			etcali+=f" --mask {options.mask}"
			if options.align_mlp:
				etcali+=f" --align_mlp --niter 10 --midout {path}/mid_{itr:02d}_{eo}.txt --decoderout {path}/dec_{itr:02d}_{eo}"
				# if itr==1:
				etcbatch+=" --niter 1"
				etcali+=" --no_load"
				# else:
				# 	etcbatch+=" --niter 2 --multibody"
			else:
				etcali+=" --align"
				etcbatch+=" --niter 0"

			pfile=f"{path}/aliptcls2d_{itr:02d}_{eo}.lst"

			run(f'e2gmm_batch.py "e2gmm_refine_jax.py --model {path}/model_{it0:02d}_{eo}.txt --ptclsin {path}/aliptcls2d_{it0:02d}_{eo}.lst --ptclsout {pfile} --maxres {res} --minres {options.minres} --spt --spt_ntilt {options.ntilt} --batchsz {options.ntilt*options.batchsize} {etcali}" {etcbatch} --batch {options.chunksize} --ptcl3d')

			if globalrefine==True:
				lst=load_lst_params(pfile)
				for l in lst:
					l["xform.projection"]=l.pop("xform.projection_00")
				save_lst_params(lst, pfile)
				run(f"e2spa_make3d.py --input {pfile} --output {path}/threed_{itr:02d}_{eo}.hdf --parallel {options.parallel} --keep 1 --sym {options.sym} {etcm3d}")

			else:
				if not os.path.isfile(f"{path}/threed_raw_{eo}.hdf"):
					run(f"e2spa_make3d.py --input {pfile} --output {path}/threed_raw_{eo}.hdf --parallel {options.parallel} --keep 1 --sym {options.sym} {etcm3d}")
				
				im=EMUtil.get_image_count(options.mask)
				for i in range(im):
					run(f"e2spa_make3d.py --input {pfile} --output {path}/threed_{itr:02d}_{i:02d}_{eo}.hdf --parallel {options.parallel} --keep 1 --sym {options.sym} {etcm3d} --xform_key xform.projection_{i:02d}")
			
				run(f"e2gmm_merge_patch.py {path} --iter {itr} --masks {options.mask}")
			
		run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf --res {res} --tophat {options.tophat} --sym {options.sym} --thread 32 --align {etcpp}")
		
		fsc=np.loadtxt("{}/fsc_masked_{:02d}.txt".format(options.path, itr))
		fi=fsc[:,1]<0.03
		nyq=2.0*er["apix_x"]
		try: 
			res=1/fsc[fi, 0][0]
		except:
			print("resolution approaching Nyquist !?")
			res=nyq
		
		if res<nyq: res=nyq
		
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	


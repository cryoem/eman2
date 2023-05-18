#!/usr/bin/env python
# Muyuan Chen 2023-04
from EMAN2 import *
import numpy as np
import time

def main():
	
	usage="""
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="Directory of the refinement.", default=None)
	parser.add_argument("--startres", type=float,help="starting resolution", default=-1)
	parser.add_argument("--npt", type=int,help="number of Gaussian function in the GMM. Ignored when --initpts is provided", default=2000)
	parser.add_argument("--initpts", type=str,help="initial Gaussian coordinates to initialize the GMMs. Take pdb or txt files.", default=None)
	parser.add_argument("--iters", type=str,help="Iteration information. Input types of refinement separated by comma. p - 3d particle translation-rotation. r - subtilt translation-rotation.", default="r,p,r")
	parser.add_argument("--setsf", type=str,help="structure factor for sharpening", default=None)
	parser.add_argument("--keep", type=float, help="",default=.9)
	parser.add_argument("--startiter", type=int, help="",default=0)
	parser.add_argument("--sym", type=str,help="sym", default="c1")

	parser.add_argument("--mask", type=str,help="mask file applied to the GMM after each iteration. The mask can apply to amplitude or sigma depending on the --maskamp or --masksigma options.", default=None)
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:32")

	parser.add_argument("--threads", type=int, help="",default=32)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	
	fname=args[0]
	print(f"Loading refinement info from {fname}")
	oldpath=os.path.dirname(fname)
	c=[i for i in fname if i.isdigit()]
	itr=int(''.join(c[-2:]))
	print(f"old refinement path {oldpath}, iteration {itr}.")
	
	
	if options.path==None: options.path=num_path_new("gmm_")
	path=options.path
	print(f"Writing new refinement in {path}...")
	
	options.cmd=' '.join(sys.argv)
	fm=f"{path}/0_spt_gmm_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	res=options.startres
	if res<=0:
		print("ERROR: Starting resolution (resolution of the given reference) needed.")
		exit()
		
	info2dname=f"{path}/particle_info_2d.lst"
	info3dname=f"{path}/particle_info_3d.lst"
	
	if options.startiter==0:
		
		#### the particle_info_xd.lst files will be used for every spt/subtilt program runs called by the refinement loop
		info2dold=f"{oldpath}/particle_info_2d.lst"
		info3dold=f"{oldpath}/particle_info_3d.lst"
		
		if os.path.isfile(info2dold) and os.path.isfile(info3dold):
			print("Using existing metadata files from the old refinement...")
			run(f"e2proclst.py {info3dold} --create {info3dname}")
			run(f"e2proclst.py {info2dold} --create {info2dname}")
		else:
			print("Error. Require valid metadata in old refinement path.")
			exit()
			
		run(f"e2proc3d.py {fname} {path}/threed_00.hdf")
		for ieo,eo in enumerate(["even", "odd"]):
			run(f"e2proc3d.py {oldpath}/threed_{itr:02d}_{eo}.hdf {path}/threed_00_{eo}.hdf")
			
			if os.path.isfile(f"{oldpath}/aliptcls2d_{itr:02d}_{eo}.lst"):
				run(f"e2proclst.py {oldpath}/aliptcls2d_{itr:02d}_{eo}.lst --create {path}/aliptcls2d_00_{eo}.lst")
			else:
				run(f"e2proclst.py {oldpath}/aliptcls2d_{itr:02d}.lst --create {path}/aliptcls2d_00_{eo}.lst --getclass {ieo}")
				
		for m in ['masked', 'unmasked', 'maskedtight']:
			fsc=np.loadtxt(f"{oldpath}/fsc_{m}_{itr:02d}.txt")
			np.savetxt(f"{path}/fsc_{m}_00.txt", fsc)
		
		if options.initpts:
			ext=options.initpts[-3:]
			if os.path.isfile(f"{options.initpts[:-4]}_even.{ext}"):
				print("using even/odd pdb")
				for ieo, eo in enumerate(["even", "odd"]):
					os.system(f"cp {options.initpts[:-4]}_{eo}.{ext} {path}/gmm_init_{eo}.{ext}")
			else:
				for ieo, eo in enumerate(["even", "odd"]):
					os.system("cp {} {}/gmm_init_{}.{}".format(options.initpts, path,eo,ext))
			
				
	#### some metadata from 2d/3d particles
	ep=EMData(info3dname,0,True)
	boxsize=ep["ny"]
	p2=EMData(info2dname,0,True)
	padsize=p2["ny"]
	
	#### parse iteration string
	it0=options.iters.split(',')
	iters=[]
	for i in it0:
		if len(i)>1:
			r=int(i[1:])
			iters.extend([i[0]]*r)
		else:
			iters.append(i)
			
	keydic={'p':"Subtomogram alignment", 't': "Subtilt translational refinement", 'T': "Subtilt translational CCF alignment", 'r': "Subtilt rotational refinement", 'd':"Defocus refinement", 'x':"Skipping alignment",'z':"Stop"}
	
		
	ppmask,setsf,tophat="","",""
	if options.setsf:
		setsf=f" --setsf {options.setsf}"
	tophat=" --tophat localwiener"
	
	#startiter=0
	minres=200
	#### now start the actual refinement loop
	#for ii,itype in enumerate(iters):
	for ii in range(options.startiter, len(iters)):
		starttime=time.time()
		itype=iters[ii]
		itr=ii+1
		ref=f"{path}/threed_{ii:02d}.hdf"
		print(f"######## iter {itr} ##########")
		print("### {}....".format(keydic[itype]))
		
		if itype=='z':break
			
		for ieo,eo in enumerate(["even", "odd"]):
			
			run(f"e2project3d.py {path}/threed_{itr-1:02d}_{eo}.hdf --outfile {path}/projections_{eo}.hdf --orientgen=eman:delta=4 --parallel=thread:12")
			
			if itr==1:
				
				if options.initpts:
					run(f"e2gmm_refine_new.py --ptclsin {path}/projections_{eo}.hdf --model {path}/gmm_init_{eo}.{ext} --maxres {res} --modelout {path}/model_{itr:02d}_{eo}.txt --niter 40 --trainmodel --learnrate 1e-6")
					
				else:
					run(f"e2segment3d.py {path}/threed_{itr-1:02d}_{eo}.hdf --pdb {path}/model_{itr:02d}_{eo}.pdb --process=segment.kmeans:nseg={options.npt}:thr=4")
					
					run(f"e2gmm_refine_new.py --ptclsin {path}/projections_{eo}.hdf --model {path}/model_{itr:02d}_{eo}.pdb --maxres {res} --modelout {path}/model_{itr:02d}_{eo}.txt --niter 40 --trainmodel --learnrate 1e-6")
		
			else:
				run(f"e2gmm_refine_new.py --ptclsin {path}/projections_{eo}.hdf --model {path}/model_{itr-1:02d}_{eo}.txt --maxres {res} --modelout {path}/model_{itr:02d}_{eo}.txt --niter 40 --trainmodel --learnrate 1e-6")

			#### subtomogram alignment. 
			if itype=='p':
				opt=""
				ptcls=info3dname
				cmd=f'e2gmm_batch.py "e2gmm_spt_align.py --ptclsin {path}/aliptcls2d_{itr-1:02d}_{eo}.lst --ptclsout {path}/aliptcls2d_{itr:02d}_{eo}.lst --model {path}/model_{itr:02d}_{eo}.txt --maxres {res} --clip {boxsize} --minres {minres} --niter 40" --batch 500 --ptcl3d --niter 0'
				
				run(cmd)
							
			#### subtilt alignment, either including the rotation or not
			##   note a subtomogram alignment need to exist first
			if itype=='r':
					
				cmd=f'e2gmm_batch.py "e2gmm_spt_subtlt.py --ptclsin {path}/aliptcls2d_{itr-1:02d}_{eo}.lst --ptclsout {path}/aliptcls2d_{itr:02d}_{eo}.lst --model {path}/model_{itr:02d}_{eo}.txt --maxres {res} --clip {boxsize} --minres {minres} --niter 20 --info3d {path}/particle_info_3d.lst" --batch 40000 --niter 0 --subtilt'
				run(cmd)
					
			#### always reconstruct 3d maps from 2d particles
			m3dpar=f" --parallel {options.parallel}"
				
			run(f"e2spa_make3d.py --input {path}/aliptcls2d_{itr:02d}_{eo}.lst --output {path}/threed_{itr:02d}_{eo}.hdf --keep {options.keep} --outsize {boxsize} --sym {options.sym} {m3dpar}")
			run(f"e2proc3d.py {path}/threed_{itr:02d}_{eo}.hdf {path}/threed_raw_{eo}.hdf --compressbits 12")
			
		run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf {setsf} {tophat} --threads {options.threads} --restarget {res:.2f} --align --sym {options.sym} {ppmask}")

		r=calc_resolution(f"{path}/fsc_masked_{itr:02d}.txt")
		res=min(r,res*1.1)		# resolution can't take too large a step in the wrong direction

		elapsed=int(time.time()-starttime)
		print(f"Iteration {ii} complete, {elapsed//3600:d}:{(elapsed%3600)//60:02d}:{elapsed%60:02d}")
	
	E2end(logid)
	
	
#### here we use 0.2 cutoff as the 0.143 one is sometimes hard to find.
##   return a slightly higher resolution to use as the maxres for the next iteration
def calc_resolution(fscfile):
	fsc=np.loadtxt(fscfile)
	
	fi=fsc[:,1]<0.2
	if np.sum(fi)==0:
		print(f"WARNING: something wrong with the FSC curve ({fscfile}). Cannot estimate resolution. Please check.")
		rs=1./fsc[-2,0]
	else:
		rs=1./fsc[fi, 0][0]
		print("Resolution (FSC<0.2) is ~{:.1f} A".format(rs))
	
	return rs*.8
	
	
if __name__ == '__main__':
	main()
	

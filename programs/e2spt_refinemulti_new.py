#!/usr/bin/env python
# Muyuan Chen 2021-04
from EMAN2 import *
import numpy as np
from e2spt_refine_new import gather_metadata

def main():
	
	usage="""
	Multi-reference classification for the new (2021) SPT refinement protocol. Since gold-standard validation is not used here, setting a --maxres is necessary. 
	e2spt_refinemulti_new.py ref1.hdf ref2.hdf --ptcls sets/ptcls.lst --niter 5 --maxres 20
	
	Alternatively, specify a number of classes instead of providing multiple references
	e2spt_refinemulti_new.py ref.hdf --ptcls sets/ptcls.lst --niter 5 --maxres 20 --nref 3
	
	If an existing single model refinement exist, provide the aliptcls3d_xx.lst file as particles and use --loadali3d. reference maps are not necessary in this case.
	e2spt_refinemulti_new.py --ptcls spt_xx/aliptcls3d_yy.lst --niter 5 --maxres 20 --nref 3 --loadali3d
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="Particle input. required", default=None)
	parser.add_argument("--path", type=str,help="Path for the refinement", default=None)
	parser.add_argument("--nref", type=int,help="Number of classes. Without --loadali3d, it duplicate the first ref N times with phase randomization at 2 x maxres. With --loadali3d, the particles are classified to N random classes at the begining.", default=-1)
	parser.add_argument("--maskalign", type=str,default=None,help="Mask file applied to 3D alignment reference in each iteration. Not applied to the average, which will follow normal masking routine.")
	parser.add_argument("--maxres",type=float,help="Maximum resolution (the smaller number) to consider in alignment (in A, not 1/A). Default is 20A",default=20.)
	parser.add_argument("--minres",type=float,help="Minimum resolution (the larger number) to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--niter", type=int,help="number of iterations. default is 5.", default=5)
	parser.add_argument("--loadali3d",help="load previous 3d alignment from --ptcls input.",action="store_true",default=False)
	#parser.add_argument("--res", type=float,help="target resolution", default=20.)
	parser.add_argument("--skipali",action="store_true",help="Skip alignment entirely when --loadali3d is provided. Otherwise a local orientation search will still be performed.",default=False)
	#parser.add_argument("--threads", type=int,help="", default=12)
	parser.add_argument("--parallel", type=str,help="parallel options", default="thread:12")
	parser.add_argument("--sym", type=str,help="symmetry to apply to the average structure", default="c1")
	parser.add_argument("--breaksym", type=str,help="Break specified symmetry. Only used when --loadali3d is on.", default=None)
	parser.add_argument("--setsf", type=str,help="set structure factor from text file", default=None)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### most of the pre-processing are copied from e2spt_refine_new
	if options.path==None: options.path=num_path_new("sptcls_")
	path=options.path
	print(f"Writing in {path}...")
	
	options.cmd=' '.join(sys.argv)
	fm=f"{path}/0_sptcls_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	info2dname=f"{path}/particle_info_2d.lst"
	info3dname=f"{path}/particle_info_3d.lst"
	
	if os.path.isfile(info2dname) and os.path.isfile(info3dname):
		print("Using existing metadata files within the directory...")
		info2d=load_lst_params(info2dname)
		info3d=load_lst_params(info3dname)
		
	else:
		info2d, info3d = gather_metadata(options.ptcls)
		save_lst_params(info3d, info3dname)
		save_lst_params(info2d, info2dname)
		
	ep=EMData(info3dname,0,True)
	boxsize=ep["ny"]
	p2=EMData(info2dname,0,True)
	padsize=p2["ny"]
		
	if options.maskalign!=None: options.maskalign=EMData(options.maskalign)
	if options.setsf!=None:
		setsf=" --setsf {}".format(options.setsf)
	else:
		setsf=""
	
	#### generating references...
	if options.loadali3d and options.nref>0:
		print("Generating initial references by random classification...")
		nref=options.nref
		ali3d=load_lst_params(options.ptcls)
		ali2d,ali3d=classify_ptcls(ali3d, info2d, options)
		save_lst_params(ali2d, f"{path}/aliptcls2d_00.lst")
		save_lst_params(ali3d, f"{path}/aliptcls3d_00.lst")
		options.ptcls=f"{path}/aliptcls3d_00.lst"
		for i in range(options.nref):
			threed=f"{path}/threed_00_{i:02d}.hdf"
			run(f"e2spa_make3d.py --input {path}/aliptcls2d_00.lst --output {threed} --keep 1 --parallel {options.parallel} --outsize {boxsize} --pad {padsize} --sym {options.sym} --clsid {i}")
			run(f"e2proc3d.py {threed} {threed} {setsf} --process filter.lowpass.gauss:cutoff_freq={1./options.maxres} --process normalize.edgemean")
	else:
		print("Loading references...")
		if options.nref>0:
			r=EMData(args[0])
			nref=options.nref
			refs=[]
			for i in range(nref):
				e=r.process("filter.lowpass.randomphase",{"cutoff_freq":1./(options.maxres*2)})
				refs.append(e)
		else:
			refs=[EMData(a) for a in args]
			nref=len(refs)
			
		for i,r in enumerate(refs):
			r.write_image(f"{path}/threed_00_{i:02d}.hdf")
		
		
	opt=""
	if options.skipali:
		opt+=" --skipali"
		
	if options.loadali3d:
		ptcls=options.ptcls
	else:
		ptcls=info3dname
		opt+=" --fromscratch"
	if options.breaksym!=None:
		opt+=" --breaksym {}".format(options.breaksym)
		
	if options.minres>0: opt+=f" --minres {options.minres}"
	if options.maxres>0: opt+=f" --maxres {options.maxres}"
		
	#### refinement loop. most of the work is dealt with in e2spt_align_subtlt. here we just parse the input/output
	for itr in range(1,1+options.niter):
		ali2d=[]
		ali3d=[]
		
		#### run alignment for each reference map
		for ir in range(nref):
			oref=f"{path}/threed_{itr-1:02d}_{ir:02d}.hdf"
			ref=f"{path}/aliref_{ir:02d}.hdf"		# overwrite each iteration

			modref=EMData(oref)
			if options.maskalign!=None: 
				# These initial filters are to reduce the artifacts produced by masking
				if options.maxres>0: modref.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.maxres})
				if options.minres>0: modref.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.minres})
				modref.mult(options.maskalign)
			modref.write_compressed(ref,0,12,erase=True)
			
			run(f"e2spt_align_subtlt.py {ptcls} {ref} --path {path} --iter {itr} --parallel {options.parallel} {opt}")
			
			ali2d.append(f"{path}/aliptcls2d_{itr:02d}_{ir:02d}.lst")
			os.rename(f"{path}/aliptcls2d_{itr:02d}.lst", ali2d[-1])
			ali3d.append(f"{path}/aliptcls3d_{itr:02d}_{ir:02d}.lst")
			os.rename(f"{path}/aliptcls3d_{itr:02d}.lst", ali3d[-1])
		
		#### parse alignment output and classify particles
		ali3dpms=[load_lst_params(a) for a in ali3d]
		ali2dpms=[load_lst_params(a) for a in ali2d]
		score=[]
		for ali in ali3dpms:
			scr=[a["score"] for a in ali]
			score.append(scr)
		score=np.array(score)
		clsid=np.argmin(score, axis=0)
		for i in np.unique(clsid):
			print("  class {} - {} particles".format(i, np.sum(clsid==i)))
		
		for ia,ali in enumerate(ali3dpms):
			for i,a in enumerate(ali):
				a["class"]=clsid[i]
				idx=info3d[i]["idx2d"]
				for x in idx:
					ali2dpms[ia][x]["class"]=clsid[i]
			
			ali3dpms[ia]=[a for a in ali if a["class"]==ia]
					
		#### save two lists (2d & 3d) for each classi in ieach iteration
		for i in range(nref):
			save_lst_params(ali2dpms[i], ali2d[i])
			save_lst_params(ali3dpms[i], ali3d[i])
		
		#### reconstruct averaged maps and filter to target resolution.
		for ir in range(nref):
			threed=f"{path}/threed_{itr:02d}_{ir:02d}.hdf"
			a2=ali2d[ir]
			run(f"e2spa_make3d.py --input {a2} --output {threed} --keep 1 --parallel {options.parallel} --outsize {boxsize} --pad {padsize} --sym {options.sym} --clsid {ir}")
			
			run(f"e2proc3d.py {threed} {threed} {setsf} --process filter.lowpass.gauss:cutoff_freq={1./options.maxres} --process normalize.edgemean")
		

	E2end(logid)
	
#### for the random particle classification at the begining of program.
def classify_ptcls(ali3d, info2d, options):
	ali2d=[]
	cls=np.arange(len(ali3d))%options.nref
	np.random.shuffle(cls)
	if options.breaksym!=None:
		xf=Transform()
		nsym=xf.get_nsym(options.breaksym)
		symidx=np.arange(len(ali3d))%nsym
		np.random.shuffle(symidx)
		for i,a in enumerate(ali3d):
			xf3d=a["xform.align3d"].inverse()
			xf3d=xf3d.get_sym(options.breaksym, int(symidx[i]))
			a["xform.align3d"]=xf3d.inverse()
		
	for p in info2d:
		pid=p["idx3d"]
		a={"src":p["src"], "idx":p["idx"],
			"tilt_id":p["tilt_id"], "ptcl3d_id":pid, "score":-1}
		xf3d=ali3d[pid]["xform.align3d"].inverse()
		
		pjxf=p["xform.projection"]*xf3d
		a["xform.projection"]=pjxf
		a["class"]=cls[pid]
		ali2d.append(a)
		
	return ali2d,ali3d
	
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

if __name__ == '__main__':
	main()
	

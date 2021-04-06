#!/usr/bin/env python
# Muyuan Chen 2021-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="3d particle input", default=None)
	parser.add_argument("--ref", type=str,help="reference map", default=None)
	parser.add_argument("--goldstandard", type=float,help="starting resolution for gold standard refinement. default 50", default=50)
	parser.add_argument("--restarget", default=0, type=float,help="The resolution you reasonably expect to achieve in the current refinement run in A.")
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--iters", type=str,help="iterations. Types of refinement separated by comma. p - 3d particle translation-rotation. t - subtilt translation. r - subtilt translation-rotation. d - subtilt defocus. Default is p,p,p,t,r,p,r,d", default="p,p,p,t,r,p,r,d")
	parser.add_argument("--keep", type=float,help="fraction to keep", default=0.95)
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:10")
	parser.add_argument("--setsf", type=str,help="structure factor for sharpening", default=None)
	parser.add_argument("--tophat", type=str,help="tophat filter options", default=None)
	parser.add_argument("--threads", type=int,help="thread for make3d and post process", default=10)
	parser.add_argument("--startres", type=float,help="starting maximum resolution. required when goldstandard is not specified", default=-1)
	parser.add_argument("--ssnrwt", action="store_true", default=False ,help="weight particles by SSNR accroding to references")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="ues the _even/_odd version of the given reference")
	parser.add_argument("--curve", action="store_true", default=False ,help="curve refinement. still testing")
	parser.add_argument("--localrefine", action="store_true", default=False ,help="only perform local search around the solution from the last iteration")
	parser.add_argument("--loadali2d", type=str,help="load previous 2d alignment", default=None)
	parser.add_argument("--loadali3d", type=str,help="load previous 3d alignment", default=None)
	parser.add_argument("--mask", type=str,help="Mask applied to the results (instead of automasking)", default=None)
	parser.add_argument("--breaksym", type=str,help="symmetry to break", default=None)
	parser.add_argument("--maskalign", type=str,help="Mask file applied to 3D alignment reference in each iteration. Not applied to the average, which will follow normal masking routine.", default=None)
	parser.add_argument("--maxshift", type=int, help="maximum shift. default box size/6",default=-1)
	parser.add_argument("--maxang", type=int, help="maximum angle difference from starting point for localrefine. ",default=30)
	parser.add_argument("--smooth",type=float,help="smooth local motion by this factor. smoother local motion with larger numbers. default 100",default=100)
	parser.add_argument("--smoothN",type=int,help="number of neighboring particles used for smoothing. default 15",default=15)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.path==None: options.path=num_path_new("spt_")
	path=options.path
	print(f"Writing in {path}...")
	
	info2dname=f"{path}/particle_info_2d.lst"
	info3dname=f"{path}/particle_info_3d.lst"
	
	if os.path.isfile(info2dname) and os.path.isfile(info3dname):
		print("Using existing metadata files within the directory...")
		
	else:
		info2d, info3d = gather_metadata(options.ptcls)
		save_lst_params(info3d, info3dname)
		save_lst_params(info2d, info2dname)
		
	options.cmd=' '.join(sys.argv)
	fm=f"{path}/0_spt_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	print("Preparing references...")
	opt=""
	if options.goldstandard>0 and options.goldcontinue==False:
		res=options.goldstandard
		refs={"even":options.ref, "odd":options.ref}
		opt+="--process filter.lowpass.gauss:cutoff_freq={r:.4f} --process filter.lowpass.randomphase:cutoff_freq={r:.4f}".format(r=1./res)
		
	else:
		res=options.startres
		if res<=0:
			print("ERROR: Starting resolution (resolution of the given reference) needed.")
			exit()
			
		if options.goldcontinue:
			refs={eo:options.ref[:-4]+f"_{eo}.hdf" for eo in ["even","odd"]}
		else:
			print("WARNING: running without even/odd spliting. this could introduce model bias...")
			refs={"even":options.ref, "odd":options.ref}
	
	
	er=EMData(refs["even"],0,True)
	ep=EMData(info3dname,0,True)
	boxsize=ep["ny"]

	if abs(1-ep["apix_x"]/er["apix_x"])>0.01 or ep["ny"]!=er["ny"]:
		print("Reference-particle apix or box size mismatch. will scale/clip reference to match particles")
		if options.mask:
			print(" Please note the given mask will NOT be rescaled. Please scale/clip manually if needed")
			
		rs=er["apix_x"]/ep["apix_x"]
		if rs>1.:
			opt+=" --clip {} --scale {} --process mask.soft:outer_radius=-1".format(ep["nx"], rs)
		else:
			opt+=" --scale {} --clip {} --process mask.soft:outer_radius=-1".format(rs, ep["nx"])
		
	if options.maskalign!=None: opt+=f" --multfile {options.maskalign}"
	
	for eo in ["even", "odd"]:
		rf=refs[eo]
		run(f"e2proc3d.py {rf} {path}/threed_00_{eo}.hdf {opt}")
		
	p2=EMData(info2dname,0,True)
	padsize=p2["ny"]
	
	last3d=last2d=None
	
	if options.loadali3d:
		fout=f"{path}/aliptcls3d_00.lst"
		run(f"e2proclst.py {options.loadali3d} --create {fout} --force")
		last3d=fout
		
	if options.loadali2d:
		fout=f"{path}/aliptcls2d_00.lst"
		run(f"e2proclst.py {options.loadali2d} --create {fout} --force")
		last2d=fout
		
	ppmask=setsf=tophat=""
	if options.setsf:
		setsf=f" --setsf {options.setsf}"
	if options.tophat:
		tophat=f" --tophat {options.tophat}"
	if options.mask:
		ppmask=f" --mask {options.mask}"
	if options.maskalign!=None: maskalign=EMData(options.maskalign,0)
	else: maskalign=None
	
	iters=options.iters.split(',')
	keydic={'p':"Subtomogram alignment", 't': "Subtilt translational refinement", 'r': "Subtilt rotational refinement", 'd':"Defocus refinement"}
	for ii,itype in enumerate(iters):
		
		itr=ii+1
		ref=f"{path}/threed_{ii:02d}.hdf"
		print(f"######## iter {itr} ##########")
		print("### {}....".format(keydic[itype]))
		
		# Ok, this is a hack to avoid adding a new option to each subprogram. May be a little confusing if the program gets interrupted
		if options.maskalign!=None:
			for eo in ("even","odd"):
				tmp=f"{path}/tmpmsk_{eo}.hdf"
				refeo=f"{path}/threed_{ii:02d}_{eo}.hdf"
				refv=EMData(refeo,0)
				refv.write_compressed(tmp,0,12)
				refv.mult(maskalign)
				refv.write_compressed(refeo,0,12)
	
		if itype=='p':
			opt=""
			if options.localrefine and last3d:
				ptcls=last3d
				opt+=f" --maxshift {options.maxshift} --maxang {options.maxang}"
			elif options.curve:
				ptcls=info3dname
				opt+=" --curve"
			else:
				ptcls=info3dname
				opt+=" --fromscratch"
			
			if last2d:
				opt+=f" --plst {last2d}"
			if options.breaksym:
				opt+=f" --breaksym {options.breaksym}"
				
			cmd=f"e2spt_align_subtlt.py {ptcls} {ref} --path {path} --iter {itr} --goldcontinue --maxres {res:.2f} --parallel {options.parallel} {opt}"
			run(cmd)
			
			last3d=f"{path}/aliptcls3d_{itr:02d}.lst"
			
		if itype=='t' or itype=='r':
			if last3d==None:
				print("Need 3D particle alignment before subtilt refinement. exit.")
				exit()
				
			cmd=f"e2spt_subtlt_local.py --ref {ref} --path {path} --iter {itr} --maxres {res} --parallel {options.parallel} --goldcontinue --refine_trans --aliptcls3d {last3d} --smooth {options.smooth} --smoothN {options.smoothN}"
			if itype=='r':
				cmd+=" --refine_rot"
			if options.maxshift>0:
				cmd+=f" --maxshift {options.maxshift}"
				
			run(cmd)
			last2d=f"{path}/aliptcls2d_{itr:02d}.lst"
				
		if itype=='d':
			if last2d==None or last3d==None:
				print("Need 3D and 2D particle alignment before defocus refinement. exit.")
				exit()
				
			cmd=f"e2spt_subtlt_local.py --ref {ref} --path {path} --iter {itr} --maxres {res} --parallel {options.parallel} --goldcontinue --refine_defocus --aliptcls3d {last3d} --aliptcls2d {last2d}  --smooth {options.smooth} --smoothN {options.smoothN}"
			
			run(cmd)
			last2d=f"{path}/aliptcls2d_{itr:02d}.lst"
			
			
		for eo in ["even", "odd"]:
			run(f"e2spa_make3d.py --input {path}/aliptcls2d_{itr:02d}.lst --output {path}/threed_{itr:02d}_{eo}.hdf --keep {options.keep} --clsid {eo} --parallel thread:{options.threads} --outsize {boxsize} --pad {padsize} --sym {options.sym}")
			
		if options.ssnrwt and itr==len(iters):
			run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf {setsf} --threads {options.threads} {ppmask}")
			res=calc_resolution(f"{path}/fsc_masked_{itr:02d}.txt")
			for eo in ["even", "odd"]:
				run(f"e2spa_make3d.py --input {path}/aliptcls2d_{itr:02d}.lst --output {path}/threed_{itr:02d}_{eo}.hdf --keep {options.keep} --clsid {eo} --parallel thread:{options.threads} --outsize {boxsize} --pad {padsize} --ref {path}/threed_{itr:02d}_{eo}.hdf --maxres {res} --sym {options.sym}")
			
		run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf {setsf} {tophat} --threads {options.threads} {ppmask} --sym {options.sym}")
		res=calc_resolution(f"{path}/fsc_masked_{itr:02d}.txt")

		# put the unmasked file back again once we finish the iteration
		if options.maskalign!=None:
			for eo in ("even","odd"):
				tmp=f"{path}/tmpmsk_{eo}.hdf"
				refeo=f"{path}/threed_{ii:02d}_{eo}.hdf"
				os.unlink(refeo)
				os.rename(tmp,refeo)
	
	E2end(logid)
	
	
def gather_metadata(pfile):
	print("Gathering metadata...")
	params=[]
	if pfile.endswith(".lst"):
		lst=LSXFile(pfile, True)
		nptcl=lst.n
		for i in range(nptcl):
			l=lst.read(i)
			params.append([l[1], l[0]])
			
		lst.close()
	
	else:
		nptcl=EMUtil.get_image_count(pfile)
		params=[[pfile, i] for i in range(nptcl)]
	
	info3d=[]
	info2d=[]
	for ii,pm in enumerate(params):
		img=EMData(pm[0], pm[1], True)
		imgsrc=img["class_ptcl_src"]
		imgidx=img["class_ptcl_idxs"]
		coord=img["ptcl_source_coord"]
		
		try: rhdrs=EMData.read_images(imgsrc,imgidx,IMAGE_UNKNOWN,True)
		except:
			print(f"couldn't read {imgidx} from {imgsrc}")
			sys.exit(1)
			
		idx2d=[]
		for k,i in enumerate(imgidx): 
			e=rhdrs[k]
			dc={"src":imgsrc,"idx":i,
				"idx3d":ii, "xform.projection":e["xform.projection"], "tilt_id":e["tilt_id"]}
			idx2d.append(len(info2d))
			info2d.append(dc)
		
		dc={"src":pm[0], "idx":pm[1], "coord":coord, "idx2d":idx2d}
		
		info3d.append(dc)

		sys.stdout.write("\r {}/{}".format(ii, len(params)))
		sys.stdout.flush()
	print()
		
	return info2d, info3d
	
	
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
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	

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
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--iters", type=str,help="iterations. Types of refinement separated by comma. p - 3d particle translation-rotation. t - subtilt translation. r - subtilt translation-rotation. d - subtilt defocus. Default is p,p,p,t,r,p,r,d", default="p,p,p,t,r,p,r,d")
	parser.add_argument("--keep", type=float,help="fraction to keep", default=0.95)
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:10")
	parser.add_argument("--setsf", type=str,help="structure factor for sharpening", default=None)
	parser.add_argument("--tophat", type=str,help="tophat filter options", default=None)
	parser.add_argument("--threads", type=int,help="thread for make3d and post process", default=10)
	parser.add_argument("--startres", type=float,help="starting maximum resolution. required when goldstandard is not specified", default=-1)
	parser.add_argument("--ssnrwt", action="store_true", default=False ,help="weight particles by SSNR accroding to references")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="ues the _even/_odd version of the given reference")
	parser.add_argument("--loadali2d", type=str,help="load previous 2d alignment", default=None)
	parser.add_argument("--loadali3d", type=str,help="load previous 3d alignment", default=None)

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
	if options.goldstandard>0:
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
		rs=er["apix_x"]/ep["apix_x"]
		if rs>1.:
			opt+=" --clip {} --scale {} --process mask.soft:outer_radius=-1".format(ep["nx"], rs)
		else:
			opt+=" --scale {} --clip {} --process mask.soft:outer_radius=-1".format(rs, ep["nx"])
		
	for eo in ["even", "odd"]:
		rf=refs[eo]
		run(f"e2proc3d.py {rf} {path}/threed_00_{eo}.hdf {opt}")
		
	p2=EMData(info2dname,0,True)
	padsize=p2["ny"]
	
	last3d=options.loadali3d
	last2d=options.loadali2d
	
	if options.setsf:
		setsf=f" --setsf {options.setsf}"
	else:
		setsf=""
		
	if options.tophat:
		tophat=f" --tophat {options.tophat}"
	else:
		tophat=""
	
	iters=options.iters.split(',')
	keydic={'p':"Subtomogram alignment", 't': "Subtilt translational refinement", 'r': "Subtilt rotational refinement", 'd':"Defocus refinement"}
	for ii,itype in enumerate(iters):
		
		itr=ii+1
		ref=f"{path}/threed_{ii:02d}.hdf"
		print(f"######## iter {itr} ##########")
		print("### {}....".format(keydic[itype]))
		
		if itype=='p':
			ptcls=f"{path}/particle_info_3d.lst"
			cmd=f"e2spt_align_subtlt.py {ptcls} {ref} --path {path} --iter {itr} --goldcontinue --fromscratch --maxres {res} --parallel {options.parallel}"
			if last2d:
				cmd+=f" --plst {last2d}"
				
			run(cmd)
			last3d=f"{path}/aliptcls3d_{itr:02d}.lst"
			
		if itype=='t' or itype=='r':
			if last3d==None:
				print("Need 3D particle alignment before subtilt refinement. exit.")
				exit()
				
			cmd=f"e2spt_subtlt_local.py --ref {ref} --path {path} --iter {itr} --maxres {res} --parallel {options.parallel} --goldcontinue --refine_trans --aliptcls3d {last3d}"
			if itype=='r':
				cmd+=" --refine_rot"
				
			run(cmd)
			last2d=f"{path}/aliptcls2d_{itr:02d}.lst"
				
		if itype=='d':
			if last2d==None or last3d==None:
				print("Need 3D and 2D particle alignment before defocus refinement. exit.")
				exit()
				
			cmd=f"e2spt_subtlt_local.py --ref {ref} --path {path} --iter {itr} --maxres {res} --parallel {options.parallel} --goldcontinue --refine_defocus --aliptcls3d {last3d} --aliptcls2d {last2d}"
			
			run(cmd)
			last2d=f"{path}/aliptcls2d_{itr:02d}.lst"
			
			
		for eo in ["even", "odd"]:
			run(f"e2spa_make3d.py --input {path}/aliptcls2d_{itr:02d}.lst --output {path}/threed_{itr:02d}_{eo}.hdf --keep {options.keep} --clsid {eo} --parallel thread:{options.threads} --outsize {boxsize} --pad {padsize}")
			
		if options.ssnrwt and itr==len(iters):
			run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf {setsf} --threads {options.threads}")
			res=calc_resolution(f"{path}/fsc_masked_{itr:02d}.txt")
			for eo in ["even", "odd"]:
				run(f"e2spa_make3d.py --input {path}/aliptcls2d_{itr:02d}.lst --output {path}/threed_{itr:02d}_{eo}.hdf --keep {options.keep} --clsid {eo} --parallel thread:{options.threads} --outsize {boxsize} --pad {padsize} --ref {path}/threed_{itr:02d}_{eo}.hdf --maxres {res}")
			
		run(f"e2refine_postprocess.py --even {path}/threed_{itr:02d}_even.hdf {setsf} {tophat} --threads {options.threads}")
		res=calc_resolution(f"{path}/fsc_masked_{itr:02d}.txt")
	
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
		
		idx2d=[]
		for i in imgidx: 
			e=EMData(imgsrc, i, True)
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
	
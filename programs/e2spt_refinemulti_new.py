#!/usr/bin/env python
# Muyuan Chen 2021-04
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="path", default=None)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--nref", type=int,help="duplicate the first ref N times with phase randomization at 2xres", default=-1)
	parser.add_argument("--niter", type=int,help="number of iterations", default=5)
	parser.add_argument("--loadali3d", type=str,help="load previous 3d alignment", default=None)
	parser.add_argument("--res", type=float,help="target resolution", default=20.)
	parser.add_argument("--skipali",action="store_true",help=".",default=False)
	parser.add_argument("--threads", type=int,help="", default=12)
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:12")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if options.path==None: options.path=num_path_new("sptcls_")
	path=options.path
	print(f"Writing in {path}...")
	
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
		

	options.cmd=' '.join(sys.argv)
	fm=f"{path}/0_spt_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	if options.nref>0:
		r=EMData(args[0])
		nref=options.nref
		refs=[]
		for i in range(nref):
			e=r.process("filter.lowpass.randomphase",{"cutoff_freq":1./(options.res*2)})
			refs.append(e)
	else:
		refs=[EMData(a) for a in args]
		nref=len(refs)
		
	for i,r in enumerate(refs):
		r.write_image(f"{path}/threed_00_{i:02d}.hdf")
		
	ep=EMData(info3dname,0,True)
	boxsize=ep["ny"]
	p2=EMData(info2dname,0,True)
	padsize=p2["ny"]
		
	opt=""
	if options.skipali:
		opt+=" --skipali"
		
	if options.loadali3d:
		ptcls=options.loadali3d
	else:
		ptcls=info3dname
		opt+=" --fromscratch"
		
		
	for itr in range(1,1+options.niter):
		ali2d=[]
		ali3d=[]
		
		for ir in range(nref):
			ref=f"{path}/threed_{itr-1:02d}_{ir:02d}.hdf"
			run(f"e2spt_align_subtlt.py {ptcls} {ref} --path {path} --iter {itr} --maxres {options.res:.2f} --parallel {options.parallel} {opt}")
			
			ali2d.append(f"{path}/aliptcls2d_{itr:02d}_{ir:02d}.lst")
			os.rename(f"{path}/aliptcls2d_{itr:02d}.lst", ali2d[-1])
			ali3d.append(f"{path}/aliptcls3d_{itr:02d}_{ir:02d}.lst")
			os.rename(f"{path}/aliptcls3d_{itr:02d}.lst", ali3d[-1])
		
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
					
		for i in range(nref):
			save_lst_params(ali2dpms[i], ali2d[i])
			save_lst_params(ali3dpms[i], ali3d[i])
		
		
		for ir in range(nref):
			threed=f"{path}/threed_{itr:02d}_{ir:02d}.hdf"
			a2=ali2d[ir]
			run(f"e2spa_make3d.py --input {a2} --output {threed} --keep 1 --parallel thread:{options.threads} --outsize {boxsize} --pad {padsize} --sym {options.sym} --clsid {ir}")
			
			run(f"e2proc3d.py {threed} {threed} --process filter.lowpass.gauss:cutoff_freq={1./options.res} --process normalize.edgemean")
		

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
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


if __name__ == '__main__':
	main()
	
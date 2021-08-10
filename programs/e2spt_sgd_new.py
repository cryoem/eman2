#!/usr/bin/env python
# Muyuan Chen 2021-08
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--niter", type=int,help="iterations", default=50)
	parser.add_argument("--parallel", type=str,help="parallel", default="thread:12")
	parser.add_argument("--shrink", type=int,help="shrink", default=1)
	parser.add_argument("--batch", type=int,help="batch size", default=12)
	parser.add_argument("--learnrate", type=float,help="learning rate", default=.2)
	parser.add_argument("--res", type=float,help="resolution", default=50)
	parser.add_argument("--ref", type=str,help="reference", default=None)
	#parser.add_argument("--sym", type=str,help="sym", default='c1')
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	ptcls=args[0]
	#ref=args[1]
			
	sys.path.insert(0, "/home/muyuanc/eman2source/eman2/programs/")
	from e2spt_refine_new import gather_metadata
	
	if options.path==None: options.path=num_path_new("sptsgd_")
	path=options.path
	
	info2dname=f"{path}/particle_info_2d.lst"
	info3dname=f"{path}/particle_info_3d.lst"
	info2d, info3d = gather_metadata(ptcls)
	save_lst_params(info3d, info3dname)
	save_lst_params(info2d, info2dname)

	res=options.res
	lr=options.learnrate
	shrink=options.shrink
	batch=options.batch
	
	e=EMData(info3dname,0)
	if shrink>1:
		e.process_inplace("math.meanshrink",{"n":shrink})
	hdr=e.get_attr_dict()
	sz=e["nx"]#//shrink
	pad=EMData(info2dname,0,True)["nx"]//shrink
	
	fname=f"{path}/output_all.hdf"
	if os.path.isfile(fname):
		os.remove(fname)
		
	def make_3d(ali2d):
		normvol=EMData(pad//2+1, pad, pad)
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"trilinear", "normout":normvol})
		recon.setup()
		for a in ali2d:
			e=EMData(a["src"],a["idx"])
			xf=Transform(a["xform.projection"])
			
			if shrink>1:
				e.process_inplace("math.meanshrink",{"n":shrink})
				xf.set_trans(-xf.get_trans()/shrink)
			
			ep=recon.preprocess_slice(e, xf)
			recon.insert_slice(ep,xf,1)

		threed=recon.finish(False)

		return normvol, threed

	def post_process(threed):
		avg=threed.do_ift()
		avg.depad()
		avg.process_inplace("xform.phaseorigin.tocenter")
		avg=avg.get_clip(Region(pad//2-sz//2,pad//2-sz//2,pad//2-sz//2, sz,sz,sz), fill=0)
		avg.set_attr_dict(hdr)
		avg.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./res})
		avg.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
		
		avg.process_inplace("math.fft.resample",{"n":1/shrink})
		avg.process_inplace("normalize.edgemean")
		avg.process_inplace("mask.soft",{"outer_radius":-10,"width":10})
		return avg

	npt=len(info3d)
	
	if options.ref==None:
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:batch])
		tt=parsesym("c1")
		xfs=tt.gen_orientations("rand",{"n":batch,"phitoo":True})
		ali2d=[]
		for ii,xf in zip(idx,xfs):
			i2d=info3d[ii]["idx2d"]
			i2d=[info2d[i] for i in i2d]
			for i in i2d:
				d={"src":i["src"],"idx":i["idx"]}
				d["xform.projection"]=i["xform.projection"]*xf
				ali2d.append(d)
			
		norm0, thrd0=make_3d(ali2d)
		avg0=post_process(thrd0)
		avg0.write_image(f"{path}/output.hdf")
		avg0.write_image(fname, -1)
		
	else:
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:batch])
		ali3d=[info3d[i] for i in idx]
		save_lst_params(ali3d, info3dname)
		
		launch_childprocess(f"e2spt_align_subtlt.py {path}/particle_info_3d.lst {options.ref} --path {path} --maxres {res} --parallel {options.parallel} --fromscratch --iter 0")
		ali2d=load_lst_params(f"{path}/aliptcls2d_00.lst")
		norm0, thrd0=make_3d(ali2d)
		avg0=post_process(thrd0)
		avg0.write_image(f"{path}/output.hdf")
		avg0.write_image(fname, -1)
	
	for itr in range(options.niter):
		print(itr)
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:batch])
		ali3d=[info3d[i] for i in idx]
		save_lst_params(ali3d, info3dname)
		
		launch_childprocess(f"e2spt_align_subtlt.py {path}/particle_info_3d.lst {path}/output.hdf --path {path} --maxres {res} --parallel {options.parallel} --fromscratch --iter 0")
		ali2d=load_lst_params(f"{path}/aliptcls2d_00.lst")
		
		
		norm1, thrd1=make_3d(ali2d)
		#avg1=post_process(thrd1)
		
		#avg1.write_image(f"{path}/output1.hdf")
		
		out=thrd0*(1-lr)+thrd1*lr
		avg=post_process(out)
		avg.write_image(f"{path}/output.hdf")
		avg.write_image(fname, -1)

		thrd0=out.copy()
		
	E2end(logid)
	
	
	
if __name__ == '__main__':
	main()
	

#!/usr/bin/env python
# Muyuan Chen 2022-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default="sgd_00")
	parser.add_argument("--res", type=int,help="res", default=30)
	parser.add_argument("--batch", type=int,help="batch", default=32)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	fname=args[0]
	ptcls=load_lst_params(fname)
	npt=len(ptcls)
	batch=options.batch
	
	
	e=EMData(fname,0, True)
	options.hdr=hdr=e.get_attr_dict()
	options.sz=sz=e["nx"]
	options.pad=pad=good_size(sz*1.5)
	options.sym='c1'
	options.shrink=1
	path=options.path
	options.parallel="thread:32"
	
	idx=np.arange(npt)
	np.random.shuffle(idx)
	idx=np.sort(idx[:batch])
	
	tt=parsesym("c1")
	xfs=tt.gen_orientations("rand",{"n":batch,"phitoo":True})
		
	ali2d=[]
	for ii,xf in zip(idx,xfs):
		i2d=ptcls[ii]
		d={"src":i2d["src"],"idx":i2d["idx"]}
		d["xform.projection"]=xf
		ali2d.append(d)
		
	thrd0=make_3d(ali2d, options)
	avg0=post_process(thrd0, options)
	avg0.write_image(f"{path}/output.hdf")
	avg0.write_compressed(f"{path}/output_all.hdf", -1, 12, nooutliers=True)
	lr=.1
	
	for itr in range(200):
		print(itr)
		
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:batch])
		pt=[ptcls[i] for i in idx]
		save_lst_params(pt, f"{path}/ptcls.lst")
		
		launch_childprocess(f"e2spa_align.py --ptclin {path}/ptcls.lst --ptclout {path}/ptcls_ali.lst --ref {path}/output.hdf --maxres {options.res} --parallel {options.parallel}")
	
		ali2d=load_lst_params(f"{path}/ptcls_ali.lst")
		
		thrd1=make_3d(ali2d, options)
		
		out=thrd0*(1-lr)+thrd1*lr
		avg=post_process(out, options)
		avg.write_image(f"{path}/output.hdf")
		avg.write_compressed(f"{path}/output_all.hdf", -1, 12, nooutliers=True)
		thrd0=out.copy()
	
	E2end(logid)
	

def make_3d(ali2d, options):
	#normvol=EMData(pad//2+1, pad, pad)
	pad=options.pad
	recon=Reconstructors.get("fourier", {"sym":options.sym,"size":[pad,pad,pad], "mode":"trilinear"})
	recon.setup()
	
	thrds=[threading.Thread(target=do_insert,args=(recon, a, options.shrink)) for a in ali2d]
	for t in thrds:  t.start()
	for t in thrds:  t.join()
		
	threed=recon.finish(False)

	return threed

def do_insert(recon, a, shrink):
	
	e=EMData(a["src"],a["idx"])
	xf=Transform(a["xform.projection"])
	xf.set_trans(-xf.get_trans())
	
	if shrink>1:
		e.process_inplace("math.meanshrink",{"n":shrink})
		xf.set_trans(xf.get_trans()/shrink)
	
	ep=recon.preprocess_slice(e, xf)
	recon.insert_slice(ep,xf,1)
	
	return

def post_process(threed, options):
	pad=options.pad
	sz=options.sz
	
	avg=threed.do_ift()
	avg.depad()
	avg.process_inplace("xform.phaseorigin.tocenter")
	avg=avg.get_clip(Region(pad//2-sz//2,pad//2-sz//2,pad//2-sz//2, sz,sz,sz), fill=0)
	avg.set_attr_dict(options.hdr)
	avg.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./options.res})
	avg.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
	
	if options.shrink>1:
		avg.process_inplace("math.fft.resample",{"n":1/options.shrink})
	avg.process_inplace("normalize.edgemean")
	avg.process_inplace("mask.soft",{"outer_radius":-10,"width":10})
	
	return avg
	
if __name__ == '__main__':
	main()
	

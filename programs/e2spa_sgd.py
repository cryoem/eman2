#!/usr/bin/env python
# Muyuan Chen 2022-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--breaksym", type=str,help="symmetry breaking", default="c1")
	parser.add_argument("--ref", type=str,help="reference file", default=None)
	parser.add_argument("--res", type=int,help="res", default=30)
	parser.add_argument("--learnrate", type=float,help="learning rate", default=0.1)
	parser.add_argument("--batch", type=int,help="batch", default=128)
	parser.add_argument("--niter", type=int,help="number of iterations", default=100)
	parser.add_argument("--ncls", type=int,help="number of classes", default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--curve", action="store_true", default=False ,help="curve mode")
	parser.add_argument("--fromscratch", action="store_true", default=False ,help="start alignment from scratch even if the lst input contains orientation")
	parser.add_argument("--skipali", action="store_true", default=False ,help="load transform from image header and skip alignment here")
	parser.add_argument("--sym", type=str,help="", default="c1")
	parser.add_argument("--parallel", type=str,help="", default="thread:32")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	fname=args[0]
	ptcls=load_lst_params(fname)
	npt=len(ptcls)
	batch=options.batch
	if options.path==None: options.path=num_path_new("sgd_")
	
	e=EMData(fname,0, True)
	options.hdr=hdr=e.get_attr_dict()
	options.sz=sz=e["nx"]
	options.pad=pad=good_size(sz*1.5)
	options.shrink=1
	path=options.path
	if options.breaksym!="c1":
		options.skipali=True
	
	save_lst_params(ptcls, f"{path}/ptcls_input.lst")
	
	idx=np.arange(npt)
	np.random.shuffle(idx)
	idx=np.sort(idx[:batch])
	
	tt=parsesym("c1")
	xfcrs=tt.gen_orientations("saff",{"delta":7.4,"phitoo":7.4,"inc_mirror":1})
	if options.curve:
		xfcrs=[x for x in xfcrs if abs(x.get_params("eman")['alt']-90)<10]
	#xfs=tt.gen_orientations("rand",{"n":batch,"phitoo":True})
	xfs=np.random.choice(xfcrs, batch)
	cls=np.arange(npt)%options.ncls
	np.random.shuffle(cls)
		
	ali2d=[]
	for ii,xf in zip(idx,xfs):
		i2d=ptcls[ii]
		d={"src":i2d["src"],"idx":i2d["idx"]}
		if options.fromscratch==False and ("xform.projection" in i2d):
			d["xform.projection"]=i2d["xform.projection"]
		else:
			d["xform.projection"]=xf
		
		d["class"]=cls[ii]
		ali2d.append(d)
		
		
	threeds=[]
	for ic in range(options.ncls):
		al=[a for a in ali2d if a["class"]==ic]
		if options.ref:
			save_lst_params(al, f"{path}/ptcls.lst")
			launch_childprocess(f"e2spa_align.py --ptclin {path}/ptcls.lst --ptclout {path}/ptcls_ali_cls{ic}.lst --ref {options.ref} --maxres {options.res} --parallel {options.parallel}")
			al=load_lst_params(f"{path}/ptcls_ali_cls{ic}.lst")
		thrd0=make_3d(al, options)
		threeds.append(thrd0)
		avg0=post_process(thrd0, options)
		avg0.write_image(f"{path}/output_cls{ic}.hdf")
		avg0.write_compressed(f"{path}/output_all_cls{ic}.hdf", -1, 8, nooutliers=True)
		
	lr=options.learnrate
	save_lst_params(ali2d, f"{path}/ptcls.lst")
	
	for itr in range(options.niter):
		print(itr)
		
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:batch])
		pt=[ptcls[i] for i in idx]
		save_lst_params(pt, f"{path}/ptcls.lst")
		
		etc=""
		if options.curve:
			etc+=" --curve"
		if options.skipali:
			etc+=" --skipali"
		if options.breaksym!="c1":
			etc+=f" --breaksym {options.breaksym}"
			
		alis=[]
		score=[]
		for ic in range(options.ncls):
			launch_childprocess(f"e2spa_align.py --ptclin {path}/ptcls.lst --ptclout {path}/ptcls_ali_cls{ic}.lst --ref {path}/output_cls{ic}.hdf --maxres {options.res} --parallel {options.parallel} {etc}")
	
			ali2d=load_lst_params(f"{path}/ptcls_ali_cls{ic}.lst")
			scr=[a["score"] for a in ali2d]
			alis.append(ali2d)
			score.append(scr)
			
		score=np.array(score)
		cls=np.argmin(score, axis=0)
		print([np.sum(cls==i) for i in range(options.ncls)])
		
		thdnew=[]
		for ic in range(options.ncls):
			ali2d=[a for a,c in zip(alis[ic], cls) if c==ic]
			print(len(ali2d))
			
		
			thrd1=make_3d(ali2d, options)
			thrd0=threeds[ic]
		
			out=thrd0*(1-lr)+thrd1*lr
			avg=post_process(out, options)
			avg.write_image(f"{path}/output_cls{ic}.hdf")
			avg.write_compressed(f"{path}/output_all_cls{ic}.hdf", -1, 8, nooutliers=True)
			thrd0=out.copy()
			thdnew.append(thrd0)
			
		threeds=thdnew
	
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
	

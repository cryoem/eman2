#!/usr/bin/env python
# Muyuan Chen 2021-08
from EMAN2 import *
import numpy as np
from e2spt_refine_new import gather_metadata

def main():
	
	usage="""prog <particle stack>
	Generate initial model from subtomogram particles using stochastic gradient descent.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=1,rowspan=1, colspan=2, mode="model")
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--res", type=float,help="target resolution", default=50,guitype='floatbox',row=2, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--niter", type=int,help="iterations", default=100, guitype='intbox',row=2, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--shrink", type=int,help="shrink", default=1, guitype='intbox',row=4, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>=<value>. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel",default="thread:17", guitype='strbox', row=6, col=1, rowspan=1, colspan=2, mode="model[thread:8]")
	parser.add_argument("--ncls", type=int,help="number of classes", default=1,guitype='intbox',row=8, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--batch", type=int,help="batch size", default=12,guitype='intbox',row=8, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--keep", type=float,help="Fraction of particles to keep. will actually align more particles and use the number of particles specified by batch", default=.7)
	parser.add_argument("--learnrate", type=float,help="learning rate, default 0.2", default=.2,guitype='floatbox',row=10, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--ref", type=str,help="Reference volume, not required. default=none", default=None,guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=11, col=1,rowspan=1, colspan=2, mode="model")
	parser.add_argument("--sym", type=str,help="symmetry. Only c1 unless --ref used", default="c1",guitype='strbox',row=12, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--classify",action="store_true",help="classify particles to the best class. there is the risk that some classes may end up with no particle. by default each class will include the best batch particles, and different classes can overlap.",default=False)
	parser.add_argument("--curve",action="store_true",help="Mode for filament structure refinement.",default=False)
	parser.add_argument("--vector",action="store_true",help="similar to --curve but keep vector direction as well.",default=False)
	parser.add_argument("--refine",action="store_true",help="only refine from existing orientations.",default=False)
	parser.add_argument("--breaksym", type=str,help="require --refine", default=None)
	parser.add_argument("--skipali",action="store_true",help="require --breaksym.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int,help="ppid", default=-2)
	parser.add_argument("--threads", type=int,help="threads", default=24)

	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	ptcls=args[0]
	
	res=options.res
	lr=options.learnrate
	shrink=options.shrink
	batch=options.batch
	ncls=options.ncls
	
	if options.path==None: options.path=num_path_new("sptsgd_")
	path=options.path
	
	nptcl=int(batch*ncls/options.keep)
	nthread=int(options.parallel.split(':')[1])
	if nthread>nptcl:
		print("Cannot use more threads than batch. Using {} threads now".format(nptcl))
		pl=options.parallel.split(':')
		pl[1]=str(nptcl)
		options.parallel=':'.join(pl)
		print(options.parallel)
	
	options.cmd=' '.join(sys.argv)
	fm=f"{path}/0_spt_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()
	
	info2dname=f"{path}/particle_info_2d.lst"
	info3dname=f"{path}/particle_info_3d.lst"
	info2d, info3d = gather_metadata(ptcls)
	save_lst_params(info3d, info3dname)
	save_lst_params(info2d, info2dname)

	
	e=EMData(info3dname,0)
	if shrink>1:
		e.process_inplace("math.meanshrink",{"n":shrink})
	options.hdr=hdr=e.get_attr_dict()
	options.sz=sz=e["nx"]#//shrink
	#options.pad=pad=EMData(info2dname,0,True)["nx"]//shrink
	options.pad=pad=good_size(sz*1.4)

	npt=len(info3d)
	
	fnames=[f"{path}/output_all_cls{ic}.hdf" for ic in range(ncls)]
	for fname in fnames:
		if os.path.isfile(fname):
			os.remove(fname)
	
	thrd0s=[]
	if options.ref==None:
		for ic in range(ncls):
			idx=np.arange(npt)
			np.random.shuffle(idx)
			idx=np.sort(idx[:batch])
			if options.curve or options.vector or options.refine:
				xfs=[]
				for ii in idx:
					s=info3d[ii]
					e=EMData(ptcls, int(ii), True)
					xf=Transform(e["xform.align3d"])
					if options.curve or options.vector:
						r=Transform({"type":"eman", "phi":np.random.rand()*360})
						xf=r*xf
					xfs.append(xf.inverse())
				
			else:
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
				
			thrd0=make_3d(ali2d, options)
			avg0=post_process(thrd0, options)
			avg0.write_image(f"{path}/output_cls{ic}.hdf")
			avg0.write_compressed(fnames[ic], -1, 12, nooutliers=True)
			thrd0s.append(thrd0)
		
	else:
		for ic in range(ncls):
			idx=np.arange(npt)
			np.random.shuffle(idx)
			idx=np.sort(idx[:nptcl])
			ali3d=[info3d[i] for i in idx]
			save_lst_params(ali3d, info3dname)
			
			launch_childprocess(f"e2spt_align_subtlt.py {path}/particle_info_3d.lst {options.ref} --path {path} --maxres {res} --parallel {options.parallel} --fromscratch --iter 0  --sym {options.sym}")
			ali2d=load_lst_params(f"{path}/aliptcls2d_00.lst")
			thrd0=make_3d(ali2d, options)
			thrd0s.append(thrd0)
			avg0=post_process(thrd0, options)
			avg0.write_image(f"{path}/output_cls{ic}.hdf")
			avg0.write_compressed(fnames[ic], -1, 12, nooutliers=True)
	
	for itr in range(options.niter):
		#print(itr)
		idx=np.arange(npt)
		np.random.shuffle(idx)
		idx=np.sort(idx[:nptcl])
		ali3d=[info3d[i] for i in idx]
		save_lst_params(ali3d, info3dname)
		
		a3dout=[]
		a2dout=[]
		for ic in range(ncls):
			print(f"iter {itr}, class {ic}: ")
			cmd=f"e2spt_align_subtlt.py {path}/particle_info_3d.lst {path}/output_cls{ic}.hdf --path {path} --maxres {res} --parallel {options.parallel} --iter 0 --sym {options.sym}"
			if options.curve:
				cmd+=" --curve"
			elif options.vector:
				cmd+=" --vector"
			elif options.refine:
				if options.breaksym:
					cmd+=f" --breaksym {options.breaksym}"
				if options.skipali:
					cmd+=" --skipali"

			else:
				cmd+=" --fromscratch"
				
			launch_childprocess(cmd)
			
			a3dout.append(load_lst_params(f"{path}/aliptcls3d_00.lst"))
			a2dout.append(load_lst_params(f"{path}/aliptcls2d_00.lst"))
		
		
		score=[]
		for ali in a3dout:
			scr=[a["score"] for a in ali]
			score.append(scr)
		score=np.array(score)
		np.savetxt(f"{path}/score.txt", score.T)
		clsid=np.argmin(score, axis=0)
		for ic in range(ncls):
			if options.classify:
				scrid=np.where(clsid==ic)[0]
				print("  class {} - {} particles".format(ic, np.sum(clsid==ic)))
			else:
				scr=score[ic].copy()
				scrid=np.argsort(scr)[:batch]
			ali2d=[]
			for a in a2dout[ic]:
				if a["ptcl3d_id"] in scrid:
					ali2d.append(a)
					
			thrd1=make_3d(ali2d, options)
			thrd0=thrd0s[ic]
			
			out=thrd0*(1-lr)+thrd1*lr
			avg=post_process(out, options)
			avg.write_image(f"{path}/output_cls{ic}.hdf")
			avg.write_compressed(fnames[ic], -1, 12, nooutliers=True)

			thrd0s[ic]=out.copy()
			
	E2end(logid)
	

def make_3d(ali2d, options):
	#normvol=EMData(pad//2+1, pad, pad)
	pad=options.pad
	recon=Reconstructors.get("fourier", {"sym":options.sym,"size":[pad,pad,pad], "mode":"trilinear"})
	recon.setup()
	
	thrds=[threading.Thread(target=do_insert,args=(recon, ali2d[i::options.threads], options.shrink)) for i in range(options.threads)]

	for t in thrds:  t.start()
	for t in thrds:  t.join()
		
	#for a in ali2d:
		#e=EMData(a["src"],a["idx"])
		#xf=Transform(a["xform.projection"])
		#xf.set_trans(-xf.get_trans())
		
		#if options.shrink>1:
			#e.process_inplace("math.meanshrink",{"n":options.shrink})
			#xf.set_trans(xf.get_trans()/options.shrink)
		
		#ep=recon.preprocess_slice(e, xf)
		#recon.insert_slice(ep,xf,1)
	
	threed=recon.finish(False)

	return threed

def do_insert(recon, a2d, shrink):
	
	for a in a2d:
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
	

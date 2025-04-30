#!/usr/bin/env python
from EMAN2 import *
from EMAN2jsondb import JSTask
import numpy as np
from scipy.optimize import minimize

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--path",type=str,help="Path to a folder where results should be stored, following standard naming conventions",default="subtlt_00")
	parser.add_argument("--ref",type=str,help="reference map",default=None)
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--maxres",type=float,help="Maximum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--minres",type=float,help="Minimum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--smooth",type=float,help="smooth local motion by this factor. smoother local motion with larger numbers. default 100",default=100)
	parser.add_argument("--smoothN",type=int,help="number of neighboring particles used for smoothing. default 15",default=15)
	parser.add_argument("--maxshift",type=float,help="max shift in pixel. default default box size/6",default=-1)
	parser.add_argument("--refine_trans",action="store_true",help="do translational alignment.",default=False)
	parser.add_argument("--refine_trans_ccf",action="store_true",help="do translational alignment using simple ccf.",default=False)
	parser.add_argument("--refine_rot",action="store_true",help="do translational-rotational alignment. better to start from an existing translational alignment.",default=False)
	parser.add_argument("--refine_defocus",action="store_true",help="do defocus refinement. need aliptcls input. doesn't work with refine_trans or rot yet..",default=False)
	parser.add_argument("--use3d",action="store_true",help="use projection of 3d particles instead of 2d ones..",default=False)
	parser.add_argument("--preprocess", metavar="processor_name:param1=value1:param2=value2", type=str, default=None, help="Preprocess each 2-D subtilt while loading (alignment only)")
	
	parser.add_argument("--aliptcls2d",type=str,help="optional aliptcls input. the program can start search from the position from last run.",default="")
	parser.add_argument("--aliptcls3d",type=str,help="optional aliptcls input.",default="")

	parser.add_argument("--range", type=str,help="process a range or 2d particles instead of the full set", default=None)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default="thread:4")
	parser.add_argument("--debug",action="store_true",help="for testing.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv, options.ppid)
	
	options.info2dname="{}/particle_info_2d.lst".format(options.path)
	options.info3dname="{}/particle_info_3d.lst".format(options.path)
	n=EMUtil.get_image_count(options.info2dname)
	if options.range:
		tasks=eval(f"range({options.range})")
		tasks=list(tasks)
	else:
		tasks=list(range(n))

	if options.preprocess!=None: options.preprocess = parsemodopt(options.preprocess)

	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel, module="e2spt_subtlt_local.SptAlignTask")
	
	num_cpus = etc.cpu_est()
	if options.debug:
		tasks=tasks[:num_cpus*4]
	print("{} jobs on {} CPUs".format(len(tasks), num_cpus))
	njob=num_cpus
	
	tids=[]
	for i in range(njob):
		t=tasks[i::njob]
		task=SptAlignTask(t, options)
		if options.debug:
			ret=task.execute(print)
			return 
		tid=etc.send_task(task)
		tids.append(tid)

	while 1:
		st_vals = etc.check_task(tids)
		if -100 in st_vals:
			print("Error occurs in parallelism. Exit")
			return
		E2progress(logid, np.mean(st_vals)/100.)
		
		if np.min(st_vals) == 100: break
		time.sleep(5)
	
	output=[None]*n
	for i in tids:
		rets=etc.get_results(i)[1]
		for r in rets:
			output[r[0]]=r[1]
		
	del etc
	
	if options.range:
		r=options.range.split(',')
		fm="{}/aliptcls2d_{:02d}_{:09d}_{:09d}.lst".format(options.path, options.iter, int(r[0]),int(r[1]))
		out=[o for o in output if o!=None]
		save_lst_params(out, fm)
		
	else:
		fm="{}/aliptcls2d_{:02d}.lst".format(options.path, options.iter)
		save_lst_params(output, fm)
	
	E2end(logid)

class SptAlignTask(JSTask):
	
	
	def __init__(self, data, options):
		
		JSTask.__init__(self,"SptAlign",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		## optimize translation for scipy minimizer
		def test_trans(tx):
			sft=np.linalg.norm(np.array(tx)-xout[3:])
			if sft>options.maxshift: return 1
			pjsft=[p.process("xform", {"tx":tx[0], "ty":tx[1]}) for p in pjsmall]

			fsc=[m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmall, pjsft)]
			fsc=np.array(fsc).reshape((len(imgsmall), 3,-1))[:,1]
			scr=np.mean(fsc[:, minp:maxp], axis=1)
			scr=np.sum(scr*wt)/np.sum(wt)
			return -scr
		
		## optimize full transform
		## consider translation of neighboring particles during rotation of center particle 
		def test_xform(xx):
			
			sft=np.linalg.norm(np.array(xx[3:])-xout[3:])
			if sft>options.maxshift: return 1
		
			dxf=Transform({"type":"xyz","xtilt":xx[0],"ytilt":xx[1],"ztilt":xx[2],
					"tx":xx[3], "ty":xx[4]})
			xfs=[dxf*x for x in xfrawsel] ## xform of center particle
			
			## now calculate extra translation of neighbor ptcls 
			dxf.set_trans(0,0,0)
			cdrot=np.array([dxf.transform(c.tolist()) for c in cdsel])
			crpj=np.array([xp.transform(c.tolist()) for c in cdrot])
			dtr=crpj-cdpj
			dtr=dtr*ss/ny
			for d,x in zip(dtr, xfs): 
				x.translate(d[0], d[1])
				
			pjsmall=make_projs(refsmall, refid, xfs)
			
			fsc=[m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmall, pjsmall)]
			fsc=np.array(fsc).reshape((len(imgsmall), 3,-1))[:,1]

			scr=np.mean(fsc[:, minp:maxp], axis=1)
			scr=np.sum(scr*wt)/np.sum(wt)
			
			return -scr
		
		callback(0)
		options=self.options
		
		### load metadata of all 2d and 3d particles
		if options.debug:
			print("loading metadata from {}, {}...".format(options.info2dname, options.info3dname))
		info3d=load_lst_params(options.info3dname)
		#info2d=load_lst_params(options.info2dname)
		if options.aliptcls3d!="":
			alipm=load_lst_params(options.aliptcls3d)
			for i,a in zip(info3d, alipm):
				i["xform.align3d"]=a["xform.align3d"]
				if "score" in a:
					i["score"]=a["score"]
				if "orig_idx" in a:
					i["orig_idx"]=a["orig_idx"]
				if "class" in a:
					i["orig_class"]=a["class"]
		
		frompast=(options.aliptcls2d!="")
		
		### load references
		if options.goldcontinue:
			refnames=[options.ref[:-4]+f"_{eo}.hdf" for eo in ["even", "odd"]]
		else:
			refnames=[options.ref]*2
			
		ref=EMData(refnames[0],0,True)
		by=ref["ny"]
		ny=by	 #ny=good_size(by*1.2)
		apix=ref["apix_x"]
		
		refs=[]
		for r in refnames:
			ref=EMData(r,0)
			ref=ref.get_clip(Region((by-ny)/2, (by-ny)/2,(by-ny)/2, ny, ny,ny))
			refft=ref.do_fft()
			refft.process_inplace("xform.phaseorigin.tocenter")
			refft.process_inplace("xform.fourierorigin.tocenter")
			refs.append(refft)
				
		#### resolution range
		if options.minres>0:
			minp=ceil(ny*ref["apix_x"]/options.minres)
		else:
			minp=2
			
		if options.maxres>0:
			maxp=ceil(ny*ref["apix_x"]/options.maxres)
			maxy=good_size(maxp*3)
			maxy=int(min(maxy, ny))
		else:
			maxy=ny
			maxp=ny//2
			
		ss=maxy
		refsmall=[]
		for r in refs:
			a=r.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
			refsmall.append(a)	
		
		if options.maxshift<0:
			options.maxshift=ss//6
		
		if options.debug: print("max res: {:.2f}, max box size {}".format(options.maxres, maxy))
		
		rets=[]
		info2dlst=LSXFile(options.info2dname,True)
		def readi(i):
			n,p,d=info2dlst.read(i)
			return {"idx":n,"src":p,"readi":i,**d}
		
		for di in self.data:
			## prepare metadata
			dc=readi(di)
			#dc=load_lst_params(options.info2dname, [di])[0]
			#dc=info2d[di] ## dictionary for the current 2d particle
			fname=dc["src"].replace("particles", "particles3d", 1)
			#fname=fname.replace("ptcls2d", "ptcls3d", 1)[:-4]+".hdf"
			tid=dc["tilt_id"]
			
			## all 3d particles from the tomogram
			d3d0=[d for d in info3d if d["src"]==fname]
			d2d=[] ## all 2d particles from the selected tilt of that tomogram
			d3d=[] ## same as d3d0, but excluding the particles that are missing on the selected tilt
			for d3 in d3d0:
				try:
					d2=[readi(i) for i in d3["idx2d"]]
					
				except:
					print("!!!!!!!!!",d3)
					continue
				#d2=[info2d[d] for d in d3["idx2d"]]
				d2=[d for d in d2 if d["tilt_id"]==tid]
				if len(d2)==0:
					## sometimes a subtilt does not exist
					continue
				d2d.append(d2[0])
				d3d.append(d3)
			
			coord=np.array([d["coord"] for d in d3d]).astype(float)
			txfs=[d["xform.align3d"].inverse() for d in d3d]
			coord-=np.array([t.get_trans() for t in txfs])
			xfpj=[d["xform.projection"] for d in d2d]
			xfraw=[a*b for a,b in zip(xfpj, txfs)]
			score=0
			#if options.debug: print(dc, d2d, d3d0, d3d)
			## the index of the selected particle from all particles on the same tilt
			ip=[i for i,d in enumerate(d2d) if d["idx3d"]==dc["idx3d"]][0]
			
			if frompast:
				ia=[d["readi"] for d in d2d]
				apm=load_lst_params(options.aliptcls2d, ia)
				xfali=[d["xform.projection"] for d in apm]
				score=apm[ip]["score"]
				#i["pastxf"]=a["xform.projection"]
				#i["score"]=a["score"]
				#xfali=[d["pastxf"] for d in d2d]
				pastxf=[b*a.inverse()for a,b in zip(xfraw, xfali)]
				pastxf=[p.get_params("xyz") for p in pastxf]
				#score=d2d[ip]["score"]
			
			## select neighboring particles (including self)
			c=coord[ip]
			dst=np.linalg.norm(coord-c, axis=1)
			srtid=np.argsort(dst)[:options.smoothN]
			dst=dst[srtid]
			wt=np.exp(-dst/options.smooth)
			
			## now load 2d particles
			xfrawsel=[Transform(xfraw[i]) for i in srtid]
			for x in xfrawsel: x.set_trans(x.get_trans()*ss/ny)
			
			d3dsel=[d3d[i] for i in srtid]
			d2dsel=[d2d[i] for i in srtid]
			if "orig_class" in d3dsel[0]:
				refid=[d["orig_class"] for d in d3dsel]
			elif "orig_idx" in d3dsel[0]:
				refid=[d["orig_idx"]%2 for d in d3dsel]
			else:
				refid=[d["idx3d"]%2 for d in d2dsel]
			
			if options.use3d:
				## use projection of 3d particles. a bit too slow...
				imgs=[EMData(d["src"], d["idx"]) for d in d3dsel]
				imgsmall=[]
				for i,img in enumerate(imgs):
					by=img["ny"]
					img=img.get_clip(Region((by-ny)/2, (by-ny)/2,(by-ny)/2, ny, ny,ny))
					x=d2dsel[i]["xform.projection"]
					img=img.project('standard',{"transform":x})
					a=img.do_fft()
					a.process_inplace("xform.phaseorigin.tocenter")
					a.process_inplace("xform.fourierorigin.tocenter")
					a=a.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
					a.process_inplace("xform.fourierorigin.tocenter")

					imgsmall.append(a)
				
			else:
				imgsmall=[]
				for d in d2dsel:
					e=EMData(d["src"], d["idx"])
					if options.preprocess!=None:
						e.process_inplace(options.preprocess[0],options.preprocess[1])
						if options.refine_trans_ccf: 
							e.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.minres})
							e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.maxres})
				
					
					e.clip_inplace(Region((e["nx"]-ny)//2, (e["ny"]-ny)//2, ny, ny))
					a=e.do_fft()
					a.process_inplace("xform.phaseorigin.tocenter")
					a.process_inplace("xform.fourierorigin.tocenter")
					a=a.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
					a.process_inplace("xform.fourierorigin.tocenter")
					imgsmall.append(a)
					
				### just read 2d particles. maybe unsafe for thick sample
				#imgs=[EMData(d["src"], d["idx"]) for d in d2dsel]
				#if options.preprocess!=None:
##					print(f"Applying {options.preprocess} to subtilts")
					#for i in imgs:
						#i.process_inplace(options.preprocess[0],options.preprocess[1])
						#if options.refine_trans_ccf: 
							#i.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/options.minres})
							#i.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.maxres})
						##i.process_inplace("mask.fft.peak",{"removepeaks":1,"thresh_sigma":1.8}) # STEVE, NPC TESTING
				#imgsmall=[]
				
				#for e in imgs:
					#e.clip_inplace(Region((e["nx"]-ny)//2, (e["ny"]-ny)//2, ny, ny))
					#a=e.do_fft()
					#a.process_inplace("xform.phaseorigin.tocenter")
					#a.process_inplace("xform.fourierorigin.tocenter")
					#a=a.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
					#a.process_inplace("xform.fourierorigin.tocenter")
					#imgsmall.append(a)
					
			#imgs=None
			## refine defocus value. does not seem to improve result yet
			if options.refine_defocus:
				if frompast:
					xfs=[Transform(xfali[i]) for i in srtid]
					for x in xfs: x.set_trans(x.get_trans()*ss/ny)
				else:
					xfs=xfrawsel
					
				pjsmall=make_projs(refsmall, refid, xfs)
				fsc=[m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmall, pjsmall)]
				fsc=np.array(fsc).reshape((len(imgsmall), 3,-1))[:,1,:-1]
				
				## invert phase back to before ctf correction first
				for i,m in enumerate(imgsmall):
					ctf=m["ctf"]
					ds=1./(apix*ny)
					c=ctf.compute_1d(m["nx"], ds, Ctf.CtfType.CTF_SIGN)
					fsc[i]*=c[:fsc.shape[1]]
				
				zrg=np.arange(-.2,.2,0.005)
				scr=[]
				for dz in zrg:
					s=[]
					for i,m in enumerate(imgsmall):
						ctf=EMAN2Ctf(m["ctf"])
						ctf.defocus=ctf.defocus+dz
						c=ctf.compute_1d(e["nx"], ds, Ctf.CtfType.CTF_SIGN)
						c=fsc[i]*c[:fsc.shape[1]]
						c=np.mean(c[minp:maxp])
						s.append(c)

					s=np.array(s)
					scr.append(-np.sum(s*wt))
					
				defocus=zrg[np.argmin(scr)]
				score=np.min(scr)
				if options.debug: print("{} - {} : {} -> {}".format(di,defocus, scr[len(scr)//2], np.min(scr)))
			else:
				defocus=0
				
			## alignment starting point
			if frompast:
				p=pastxf[ip].copy()
				xout=[0.,0.,0.,p["tx"]*ss/ny,p["ty"]*ss/ny]
			else:
				xout=[0.0]*5
			
			xout00=np.array(xout).copy()
			## translational alignment
			if options.refine_trans:
				pjsmall=make_projs(refsmall, refid, xfrawsel)
				res=minimize(test_trans, [xout[3], xout[4]],
						method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
				
				xout[3:]=res.x
				score=res.fun
				if options.debug: print("{} - {} -> {} : {:.5f} -> {:.5f}".format(di, np.round(xout00[3:],4),np.round(xout[3:],4), test_trans([xout00[3], xout00[4]]), res.fun))


			## translational alignment using ccf aligner
			if options.refine_trans_ccf:
				pjsmall=make_projs(refsmall, refid, xfrawsel)
				
				ali=imgsmall[0].align("translational",pjsmall[0],{"maxshift":options.maxshift})
				xf=ali["xform.align2d"]
				xout[3],xout[4]=xf.get_trans_2d()
				xout[3]*=-1
				xout[4]*=-1
				score=ali["score_align"]
			
			## now full transform alignment 
			if options.refine_rot:
				xp=xfpj[srtid[0]]
				cdsel=coord[srtid].copy()
				cdsel-=cdsel[0]
				cdpj=np.array([xp.transform(c.tolist()) for c in cdsel])
				
				y0=test_xform(xout)
				res=minimize(test_xform, xout,
						method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})
						
				xout=res.x	
				score=res.fun	
			
				if options.debug: print("{} - {} -> {} : {:.5f} -> {:.5f}".format(di, np.round(xout00,4),np.round(xout,4), y0, res.fun))
				
			## prepare output
			xf=Transform(xfraw[ip])
			dxf=Transform({"type":"xyz","xtilt":xout[0],"ytilt":xout[1],"ztilt":xout[2],
						"tx":xout[3],"ty":xout[4]})
			dxf.set_trans(dxf.get_trans()*ny/ss)
			xf=dxf*xf
			if options.debug: 
				print(dxf)
				print(xf)
			d3=info3d[dc["idx3d"]]
			if "orig_class" in d3:
				clsid=d3["orig_class"] 
			elif "orig_idx" in d3:
				clsid=d3["orig_idx"]%2
			else: 
				clsid=dc["idx3d"]%2
			
			c={
				"src":dc["src"], 
				"idx":dc["idx"],
				"xform.projection":xf, 
				"score":score,
				"class":clsid, 
				"defocus":defocus, 
				"dxf":dxf,
				"tilt_id":tid,
				"ptcl3d_id":dc["idx3d"]
				}
			
			rets.append((di,c))

			callback(len(rets)*100//len(self.data))
		
		return rets

def make_projs(refs, rid, xfs):
	pjs=[]
	for ci,xf in enumerate(xfs):
		ref=refs[rid[ci]]
		pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
		pjs.append(pj)
		
	return pjs


if __name__ == "__main__":
	main()


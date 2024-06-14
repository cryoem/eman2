#!/usr/bin/env python
# Muyuan Chen 2018-04
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2jsondb import JSTask
from EMAN2PAR import EMTaskCustomer
from scipy.optimize import minimize

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ptclin", type=str,help="particle input", default=None)
	parser.add_argument("--ptclout", type=str,help="particle output", default=None)
	parser.add_argument("--ref", type=str,help="reference input", default=None)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12")

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data")
	parser.add_argument("--curve", action="store_true", default=False ,help="curve mode")
	#parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=-1)
	parser.add_argument("--localrefine", type=int, default=-1 ,help="local refinement. larger value correspond to smaller local region")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="split even/odd subset and references.")
	parser.add_argument("--extrarot", action="store_true", default=False ,help="do an extra rotation to combat the interpolation issues.")
	parser.add_argument("--skipali",action="store_true",help="Skip alignment and only calculate the score. Incompatible with --fromscratch, but --breaksym will still be considered.",default=False)
	#parser.add_argument("--skiprefine", action="store_true", default=False ,help="coarse search only.")
	#parser.add_argument("--ctfweight", action="store_true", default=False ,help="weight by ctf. not used yet...")
	#parser.add_argument("--slow", action="store_true", default=False ,help="slow but finer search")
	parser.add_argument("--maxres", type=float,default=-1, help="max resolution for cmp")
	parser.add_argument("--minrespx", type=int,default=4, help="skip the first x pixel in fourier space")
	parser.add_argument("--sym", type=str,help="symmetry. ", default="c1")
	parser.add_argument("--breaksym", type=str,help="breaking symmetry. only works for --localrefine or --skipali", default="c1")

	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	parser.add_argument("--alt_min", type=int,help="alt min", default=-1)
	parser.add_argument("--alt_max", type=int,help="alt ax", default=-1)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#m=EMData(options.ref)
	#bxsz=m["nx"]
	#apix=m["apix_x"]
	
	options.shrink=1
	pinfo=load_lst_params(options.ptclin)
	nptcl=len(pinfo)
	#if options.maxshift<0:
		#options.maxshift=bxsz//2
	
	print("Initializing parallelism...")
	etc=EMTaskCustomer(options.parallel, module="e2spa_align.SpaAlignTask")	
	num_cpus = etc.cpu_est()
	
	print("{} particles".format(nptcl))
	print("{} total CPUs available".format(num_cpus))
	if options.debug: 
		nptcl=min(4*num_cpus, nptcl)
		print("Debugging mode. running on one thread with 8 particles")
		
	
	nbatch=min(nptcl, num_cpus)
	
	infos=[[] for i in range(nbatch)]
	for i,info in enumerate(pinfo):
		infos[i%nbatch].append([i, info])
		
	print("{} jobs, each with {:.1f} particles".format(len(infos), np.mean([len(i) for i in infos])))
	
	tids=[]
	for info in infos:
		task = SpaAlignTask(info, options)
			
		if options.debug:
			task.execute(print)
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
	
	output=[None]*nptcl
	for i in tids:
		ret=etc.get_results(i)[1]
		for r in ret:
			output[r[0]]=r[1]
	
	del etc
	
	fm=options.ptclout
	save_lst_params(output, fm)
	
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	

class SpaAlignTask(JSTask):
	
	
	def __init__(self, info, options):
		
		data={"info":info}
		JSTask.__init__(self,"SpaAlignTask",data,{},"")
		self.options=options
	
	def execute(self, callback):
		time0=time.time()
		
		def test_rot(x, returnxf=False):
			if isinstance(x, Transform):
				xf=x
			else:
				x=list(x)
				x.extend(curxf[len(x):])
				xf=Transform({"type":"eman", "az":x[0], "alt":x[1], "phi":x[2],"tx":x[3], "ty":x[4]})
			
			v=np.random.randn(3)
			v/=np.linalg.norm(v)
			tiny=Transform({"type":"spin","omega":.1,"n1":v[0],"n2":v[1],"n3":v[2]})
			xf=tiny*xf
			
			pj=refsmall.project('gauss_fft',{"transform":xf, "returnfft":1})
			x0=options.minrespx; x1=int(ss*.4)
			xf2=Transform(xf)

			ccf=imgsmall.calc_ccf(pj)
			pos=ccf.calc_max_location_wrap(mxsft, mxsft, 0)
			xf.translate(pos[0], pos[1],0)
			pj.process_inplace("xform", {"tx":pos[0], "ty":pos[1]})
		
			fsc=imgsmall.calc_fourier_shell_correlation(pj)
			fsc=np.array(fsc).reshape((3,-1))[:, x0:x1]
				
			wt=np.ones_like(fsc[1])
			if ctfwt:
				wt*=ctfcv[x0:x1]
			
			scr=-np.sum(fsc[1]*wt)/np.sum(wt)
			
			
			if options.extrarot:
				
				pj2=refrotsmall.project('gauss_fft',{"transform":xf2*r45, "returnfft":1})
				ccf=imgsmall.calc_ccf(pj2)
				pos=ccf.calc_max_location_wrap(mxsft, mxsft, 0)
				pj2.process_inplace("xform", {"tx":pos[0], "ty":pos[1]})
			
				fsc2=imgsmall.calc_fourier_shell_correlation(pj2)
				fsc2=np.array(fsc2).reshape((3,-1))[:, x0:x1]
				
				scr2=-np.sum(fsc2[1]*wt)/np.sum(wt)
				scr=(scr*.5+scr2*.5)
			
			if options.curve:
				alt=xf.get_params("eman")["alt"]
				a=abs(alt-90.)/10
				scr+=a**2
				
			if returnxf:
				return scr, xf
			else:
				return scr
		
		
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		
		#### do the even/odd split if required
		reffile=options.ref
		if options.goldcontinue:
			refnames=[reffile[:-4]+f"_{eo}.hdf" for eo in ["even", "odd"]]
		else:
			refnames=[reffile, reffile]
			
		refs=[]
		for r in refnames:
			ref=EMData(r)
			ref=ref.do_fft()
			ref.process_inplace("xform.phaseorigin.tocenter")
			ref.process_inplace("xform.fourierorigin.tocenter")
			refs.append(ref)
			
		if options.extrarot:
			refsrot=[]
			r45=Transform({"type":"eman","alt":45,"az":45})
			for r in refnames:
				ref=EMData(r)
				ref.process("xform",{"transform":r45})
				ref=ref.do_fft()
				ref.process_inplace("xform.phaseorigin.tocenter")
				ref.process_inplace("xform.fourierorigin.tocenter")
				refsrot.append(ref)
			
		
		ref=EMData(refnames[0], 0, True)
		ny0=ny=ref["ny"]
		if options.maxres>0:
			maxrescut=ceil(ny*ref["apix_x"]/options.maxres)
			maxy=good_size(maxrescut*3)
			maxy=int(min(maxy, ny))
		else:
			maxy=ny
			maxrescut=1e5
		
		ssrg=2**np.arange(4,12, dtype=int)
		ssrg[:2]=ssrg[:2]*3/2
		ssrg=ssrg.tolist()
		
		astep=7.4
		sym=Symmetries.get(options.sym)
		
		if options.curve:
			xfcrs=sym.gen_orientations("eman",{"delta":astep,"phitoo":0,"inc_mirror":1,"alt_min":80, "alt_max":100})
		elif options.alt_max>0:
			xfcrs=sym.gen_orientations("eman",{"delta":astep,"phitoo":astep,"inc_mirror":1, "alt_min":options.alt_min, "alt_max":options.alt_max})
			if options.alt_min==0:
				xs=sym.gen_orientations("eman",{"delta":astep,"phitoo":astep,"inc_mirror":1, "alt_min":180-options.alt_max, "alt_max":180})
				xfcrs+=xs
		else:
			xfcrs=sym.gen_orientations("eman",{"delta":astep,"phitoo":astep,"inc_mirror":1})
		#print(np.unique([x.get_params("eman")['alt'] for x in xfcrs]))
		
		for infoi, infos in enumerate(data["info"]):
			ii=infos[0]
			info=infos[1]
			if "class" in info:
				clsid=info["class"]
			else:
				clsid=ii%2
			
			img=EMData(info["src"],info["idx"])
			img.clip_inplace(Region((img["nx"]-ny)//2, (img["ny"]-ny)//2, ny, ny))
			img.process_inplace("mask.soft",{"outer_radius":-10,"width":10})
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			img.process_inplace("xform.fourierorigin.tocenter")
			
			npos=32
			istart=0
			initxfs=[]
				
			if "xform.projection" in info:
				initxfs.append(info["xform.projection"])
					
			if options.localrefine>0:	
				newxfs=[]
				istart=min(options.localrefine, len(ssrg)-1)
				npos=npos//(2**istart)
				npos=max(npos, 1)
				
				for ixf in initxfs:
					ix=ixf.get_params("eman")
					newxfs.append(Transform(ix))

					for i in range(npos-1):
						d={"type":"eman","tx":ix["tx"], "ty":ix["ty"]}
						for ky in ["az", "alt", "phi"]:
							d[ky]=ix[ky]+np.random.randn()*5./np.pi*2
						newxfs.append(Transform(d))
				
			if options.skipali:
				### skip alignment alltogether
				newxfs=initxfs
				istart=len(ssrg)-1
				npos=1
				
			if options.breaksym!="c1":
				xf2=[]
				for xf in newxfs:
					for i in range(Transform.get_nsym(options.breaksym)):
						x2=xf.get_sym(options.breaksym, i)
						xf2.append(x2)
				newxfs=xf2
				npos=len(newxfs)
				
			if img.has_attr("ctf"):
				ctf=img["ctf"]
				ds=1./(ny*ref["apix_x"])
				ctf.bfactor=10
				ctfcv=abs(np.array(ctf.compute_1d(ny,ds,Ctf.CtfType.CTF_AMP)))
				try:		# this may fail for some neg stain images close to focus, failsafe to disable weighting
					ci=np.where(np.diff(ctfcv)<0)[0][0]
					ctfcv[:ci]=ctfcv[ci]
					ctfwt=True
				except:
					ctfwt=False
			else:
				ctfwt=False
				
			for si in range(istart, len(ssrg)):
				ss=ssrg[si]
				if ss>=maxy: 
					ss=maxy
					
				refsmall=refs[clsid].get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
				if options.extrarot:
					refrotsmall=refsrot[clsid].get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
				imgsmall=img.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
				imgsmall.process_inplace("xform.fourierorigin.tocenter")
					
				mxsft=ss//8
				astep=89.999/floor((np.pi/(3*np.arctan(2./ss))))
				
				score=[]
				
				if options.debug:
					print(ss, npos, astep)
					
				if si==0:
					
					newxfs=[]
					for xf in xfcrs:
						scr, x=test_rot(xf, True)
						score.append(scr)
						newxfs.append(x)
				
				else:
					xfs=newxfs
					if len(initxfs)>0 and options.skipali==False:
						xfs.append(initxfs[0])
					newxfs=[]
					for xf0 in xfs:
						x=xf0.get_params("eman")
						curxf=[x["az"], x["alt"], x["phi"],x["tx"]*ss/ny,x["ty"]*ss/ny]
						x0=curxf[:3]
						if options.skipali:
							scr, x=test_rot(curxf, True)
						else:
							res=minimize(test_rot, x0, method='Powell', options={'ftol': 1e-3, 'disp': False, "maxiter":20})
							scr, x=test_rot(res.x, True)
						score.append(scr)
						newxfs.append(x)
						
				xfs=newxfs

				newxfs=[]
				newscore=[]
				idx=np.argsort(score)
				for i in idx:
					dt=[(x.inverse()*xfs[i]).get_params("spin")["omega"] for x in newxfs]
					if len(dt)==0 or np.min(dt)>astep*4:
						newxfs.append(xfs[i])
						newscore.append(score[i])
					if len(newxfs)>=npos:
						break
					
				for xi,xf in enumerate(newxfs):
					## translation between -ny//2 and ny//2
					x=np.array(xf.get_trans()*ny/ss)
					x=x%ny
					x=x-ny*(x>ny//2)
					
					xf.set_trans(x.tolist())
					
					
				if options.verbose>1:
					print("size: {}, xfs: {}".format(ss, len(newxfs)))
					for x in newxfs:
						xf=x.get_params("xyz")
						print("\t{:.1f} {:.1f} {:.1f} {:.1f} {:.1f}".format(xf["xtilt"], xf["ytilt"], xf["ztilt"], xf['tx'], xf['ty']))
				
				npos=max(1, npos//2)
				lastastep=astep
				if ss>=maxy:
					break
				
			r={	"src":info["src"], "idx":info["idx"], 
				"xform.projection":newxfs[0], "score":newscore[0],"class":clsid}
			callback(100*float(infoi/len(self.data["info"])))
			rets.append((ii, r))
			if options.debug:
				print('time',time.time()-time0)
				print(newxfs[0])
				return
			
		return rets


	
if __name__ == '__main__':
	main()

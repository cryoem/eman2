#!/usr/bin/env python
# Muyuan Chen 2018-04
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2jsondb import JSTask
import re
from EMAN2_utils import make_path
from shutil import copy2
from EMAN2PAR import EMTaskCustomer
from scipy.optimize import minimize

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ptclin", type=str,help="particle input", default=None)
	parser.add_argument("--ptclout", type=str,help="particle output", default=None)
	parser.add_argument("--ref", type=str,help="reference input", default=None)
	parser.add_argument("--threedout", type=str,help="map output", default=None)

	parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.5", default=0.5)


	parser.add_argument("--threads", type=int,help="Number of CPU threads to use. Default is 12.", default=12)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12")

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data (threads * 8 particles)")
	

	parser.add_argument("--transonly", action="store_true", default=False ,help="only refine translation")
	parser.add_argument("--savepath", action="store_true", default=False ,help="save alignment path in a json file for testing.")
	parser.add_argument("--scipytest", action="store_true", default=False ,help="test scipy optimizer")
	parser.add_argument("--fromscratch", action="store_true", default=False ,help="align from scratch and ignore previous particle transforms. for spt mostly. will include mirror")
	parser.add_argument("--refineastep", type=float,help="Mean angular variation for refine alignment", default=2.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=4)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=8)
	parser.add_argument("--localrefine", action="store_true", default=False ,help="local refinement")
	parser.add_argument("--defocus", action="store_true", default=False ,help="refine defocus. Still under development")
	parser.add_argument("--seedmap", action="store_true", default=False ,help="seed")


	parser.add_argument("--padby", type=float,default=2.0, help="pad by factor. default is 2")
	parser.add_argument("--maxres", type=float,default=-1, help="max resolution for cmp")
	parser.add_argument("--minres", type=float,default=-1, help="min resolution for cmp")
	parser.add_argument("--sym", type=str,help="symmetry. will use symmetry from spt refinement by default", default="c1")
	parser.add_argument("--smooth", type=int,help="Smooth trajectory per image based on nearby particles. Still under development", default=-1)
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	parser.add_argument("--nkeep", type=int,help="", default=1)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)


	lstname=options.ptclin
	threedname=options.ref
	lname=options.ptclout
	threedout=options.threedout
	
	lst=LSXFile(lstname, True)
	m=EMData(threedname)
	bxsz=m["nx"]
	apix=m["apix_x"]
	
	options.shrink=1
	pinfo=[]
	if options.debug: nptcl=options.threads*8
	else: nptcl=lst.n
	for i in range(nptcl):
		pinfo.append(lst.read(i))
	lst=None
	
	print("Initializing parallelism...")
	if options.scipytest:
		etc=EMTaskCustomer(options.parallel, module="e2spt_tiltrefine_oneiter.SptNewTltRefineTask")
	else:
		etc=EMTaskCustomer(options.parallel, module="e2spt_tiltrefine_oneiter.SptTltRefineTask")
		
	num_cpus = etc.cpu_est()
	
	print("{} total CPUs available".format(num_cpus))
	print("{} jobs".format(nptcl))
	infos=[[] for i in range(num_cpus)]
	for i,info in enumerate(pinfo):
		infos[i%num_cpus].append([i, info])
	
	tids=[]
	for info in infos:
		if options.scipytest:
			task = SptNewTltRefineTask(info, threedname, options)
		else:
			task = SptTltRefineTask(info, threedname, options)
			
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
	
	dics=[0]*nptcl
	for i in tids:
		ret=etc.get_results(i)[1]
		for r in ret:
			#print(r)
			ii=r.pop("idx")
			dics[ii]=r
	
	del etc
	
	if options.savepath:
		dtmp={i:dics[i] for i in range(len(dics))}
		jstmp=js_open_dict(lname.replace('.lst', '.json'))
		jstmp.update(dtmp)
		jstmp=None
	
	
	allscr=[d["score"] for d in dics]
	if options.smooth>0 or options.defocus:
		### need to add per tilt smoothing later...
		s=np.array(allscr)
		np.savetxt(lname.replace(".lst", "_score.txt"), s)
		print(s.shape)
		#return
	
	maxl=np.max([len(s) for s in allscr])
	maxv=np.max(np.concatenate(allscr))
	for s in allscr:
		s.extend([maxv]*(maxl-len(s)))
	allscr=np.array(allscr)
	
	#print(np.min(allscr), np.mean(allscr), np.max(allscr), np.std(allscr))
	allscr=2-allscr
	allscr-=np.min(allscr)
	allscr/=np.max(allscr)
	
	if maxl>1:
		mx=np.max(allscr, axis=1)[:,None]
		allscr=np.exp(allscr*20)
		allscr=allscr/np.sum(allscr, axis=1)[:,None]
		allscr*=mx
	
	try: os.remove(lname)
	except: pass
	lout=LSXFile(lname, False)
	for i, dc in enumerate(dics):
		lc=""
		
		if isinstance(dc["xform.align3d"], list):
			alilist=dc["xform.align3d"]
			scorelist=dc["score"]
		else:
			alilist=[dc["xform.align3d"]]
			scorelist=[dc["score"]]
		for j,xf in enumerate(alilist):
			d=xf.get_params("eman")
			d["score"]=float(allscr[i,j])
			if d["score"]>.05 or j==0:
				lc=lc+str(d)+';'
			
		l=pinfo[i]
		lout.write(-1, l[0], l[1], lc[:-1])

	lout=None

	pb=options.padby
	
	if options.parallel.startswith("mpi") and len(dics)>10000:
		m3dpar="--parallel {}".format(options.parallel)
	else:
		m3dpar=""
		
	if options.seedmap:
		seed="--seedmap "+threedname
	else:
		seed=""
	
	cmd="e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads {trd} --outsize {bx} --apix {apx} --mode gauss_2 --keep {kp} --sym {sm} {seed} {par}".format(
		inp=lname, 
		out=threedout,
		bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=options.sym, trd=options.threads,par=m3dpar, seed=seed)
	
	run(cmd)

	E2end(logid)


def run(cmd):
	print(cmd)
	launch_childprocess(cmd)


def gen_xfs(symc, astep):
	sym=parsesym(symc)
	xfs0=sym.gen_orientations("saff", {"delta":astep,"inc_mirror":1})
	xfs=[x for x in xfs0]
	n0=len(xfs)
	ns=Transform.get_nsym(symc)

	for i in range(1,ns):
		xfs.extend([xf.get_sym(symc, i) for xf in xfs0])

		
	xfmat=np.array([x.get_matrix() for x in xfs])
	xfmat=xfmat.reshape((-1, 3,4 ))[:, :,:3]
	dt=np.tensordot(xfmat, xfmat, axes=(1,1)).transpose(0,2,1,3)
	om=(np.trace(dt, axis1=2, axis2=3)-1)/2.
	om=np.clip(om, -1, 1)
	om=np.arccos(om)*180/np.pi
	dst=om
	dst+=np.eye(len(xfs))*180
	dt=np.min(dst, axis=1)

	mid=np.where(dt[:n0]>astep/.9)[0]
	nid=np.argmin(dst[mid], axis=1)
	xns=[]
	for i, mi in enumerate(mid):
		ni=nid[i]
		dx=(xfs[mi]*xfs[ni].inverse()).get_params('spin')
		dx["omega"]/=-2
		xns.append(xfs[mi]*Transform(dx))

	xfs=xfs[:n0]
	xfs.extend(xns)
	xfsall=[]
	for x in xfs:
		for phi in np.arange(0,360-.1, astep):
			x0=x.get_params("eman")
			x0["phi"]=phi
			xfsall.append(Transform(x0))
	return xfsall

class SptNewTltRefineTask(JSTask):
	
	
	def __init__(self, info, ref, options):
		
		data={"info":info, "ref": ref}
		JSTask.__init__(self,"SptNewTltRefine",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		time0=time.time()
		def test_rot(x, returnxf=False):
			rotcount.append(0)
			fullxf=False
			if isinstance(x, Transform):
				xf=x
				
			else:
				if len(x)<4:
					xf=Transform({"type":"eman", "alt":x[0], "az":x[1], "phi":x[2]})
				else:
					xf=Transform({"type":"eman", "alt":x[0], "az":x[1], "phi":x[2],"tx":x[3], "ty":x[4]})
					fullxf=True			   
			
			pj=refsmall.project('gauss_fft',{"transform":xf, "returnfft":1})
			
			x0=ss//8; x1=int(ss*.4)
			#x0=ss//8; x1=ss//3
			ccf=imgsmall.calc_ccf(pj)
			pos=ccf.calc_max_location_wrap(mxsft, mxsft, 0)
			xf.set_trans(pos)
			pj.translate(xf.get_trans())
		
			fsc=imgsmall.calc_fourier_shell_correlation(pj)
			fsc=np.array(fsc).reshape((3,-1))[:, x0:x1]
			wt=fsc[2]
			if ctfwt:
				wt*=ctfcv[x0:x1]
			
			scr=-np.sum(fsc[1]*wt)/np.sum(wt)
			
			if returnxf:
				return scr, xf
			else:
				return scr
		
		
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		

		ref=EMData(data["ref"],0)
		ny0=ny=ref["ny"]
		
		if options.maxres>0:
			maxrescut=ceil(ny*ref["apix_x"]/options.maxres)
			maxy=good_size(maxrescut*3)
			maxy=int(min(maxy, ny))
		else:
			maxy=ny
			maxrescut=1e5
		
		ref=ref.do_fft()
		ref.process_inplace("xform.phaseorigin.tocenter")
		ref.process_inplace("xform.fourierorigin.tocenter")
		ssrg=2**np.arange(4,12, dtype=int)
		ssrg[:2]=ssrg[:2]*3/2
		ssrg=ssrg.tolist()
		
		for infoi, infos in enumerate(data["info"]):
			rotcount=[]
			ii=infos[0]
			info=infos[1]
			
			img=EMData(info[1],info[0])
			if img["ny"]!=ny: # box size mismatch. simply clip the box
				by=img["ny"]
				img=img.get_clip(Region((by-ny)/2, (by-ny)/2, ny,ny))
				
			img.process_inplace("mask.soft",{"outer_radius":-10,"width":10})
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			img.process_inplace("xform.fourierorigin.tocenter")
			path=[]
			
			npos=32
			istart=0
			initxfs=[]
				
			if not options.fromscratch and isinstance(info[-1], str):
				for xfs in info[-1].split(';'):
					initxf=eval(xfs)
					if "score" in initxf:
						initxf.pop('score')
					initxfs.append(Transform(initxf))
					
			if options.localrefine:	
				newxfs=[]
				istart=1
				npos=npos//2
				for ixf in initxfs:
					ix=ixf.get_params("eman")
					for i in range(npos):
						d={"type":"eman","tx":ix["tx"], "ty":ix["ty"]}
						for ky in ["alt", "az", "phi"]:
							d[ky]=ix[ky]+np.random.randn()*5./np.pi*2
						newxfs.append(Transform(d))
				
			if img.has_attr("ctf"):
				ctf=img["ctf"]
				ds=1./(ny*ref["apix_x"])
				ctf.bfactor=10
				ctfcv=abs(np.array(ctf.compute_1d(ny,ds,Ctf.CtfType.CTF_AMP)))
				ci=np.where(np.diff(ctfcv)<0)[0][0]
				ctfcv[:ci]=ctfcv[ci]
				ctfwt=True
			else:
				ctfwt=False
				
			for si in range(istart, len(ssrg)):
				ss=ssrg[si]
				if ss>=maxy: 
					ss=maxy
					
				refsmall=ref.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
				imgsmall=img.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
				imgsmall.process_inplace("xform.fourierorigin.tocenter")
					
				mxsft=ss//8
				astep=89.999/floor((np.pi/(3*np.arctan(2./ss))))
				sym=Symmetries.get(options.sym)
				score=[]
				
				if options.debug:
					print(ss, npos, astep)
					
				if si==0:
					xfs=sym.gen_orientations("saff",{"delta":astep,"phitoo":astep,"inc_mirror":1})
					newxfs=[]
					for xf in xfs:
						scr, x=test_rot(xf, True)
						score.append(scr)
						newxfs.append(x)
				
				else:
					xfs=newxfs
					if len(initxfs)>0:	
						xfs.append(initxfs[0])
					newxfs=[]
					simplex=np.vstack([[0,0,0], np.eye(3)*astep])
					for xf0 in xfs:
						x=xf0.get_params("eman")
						x0=[x["alt"], x["az"], x["phi"]]
						res=minimize(test_rot, x0, method='Nelder-Mead', options={'ftol': 1e-2, 'disp': False, "maxiter":50, "initial_simplex":simplex+x0})
						scr, x=test_rot(res.x, True)
						score.append(scr)
						newxfs.append(x)
						
				xfs=newxfs
				if options.debug:
					print(ss, "xfs:", len(newxfs))

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
					
				for xf in newxfs:
					xf.set_trans(xf.get_trans()*ny/ss)
					
				if options.verbose>1:
					print("size: {}, xfs: {}".format(ss, len(newxfs)))
					for x in newxfs:
						xf=x.get_params("eman")
						print("\t{:.1f} {:.1f} {:.1f} {:.1f} {:.1f}".format(xf["alt"], xf["az"], xf["phi"], xf['tx'], xf['ty']))
				
				npos=max(1, npos//2)
				lastastep=astep
				path.append(newxfs)
				if ss>=maxy:
					break
				
			#r={"idx":ii, "xform.align3d":newxfs[0], "score":np.min(score)}
			r={"idx":ii, "xform.align3d":newxfs[:options.nkeep], "score":newscore[:options.nkeep], "path":path}
			callback(100*float(infoi/len(self.data["info"])))
			rets.append(r)
			if options.debug:
				print('time',time.time()-time0, len(rotcount))
			
		return rets


class SptTltRefineTask(JSTask):
	
	
	def __init__(self, info, ref, options):
		
		data={"info":info, "ref": ref}
		JSTask.__init__(self,"SptTltRefine",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		time0=time.time()
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		a=EMData(data["ref"],0)
		
		for infoi, infos in enumerate(data["info"]):
			ii=infos[0]
			info=infos[1]
			b=EMData(info[1],info[0])
				
			if options.maxres>0:
				if options.shrink>1:
					b.process_inplace("math.meanshrink",{"n":options.shrink})
				b.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1./options.maxres})
			
			if b["ny"]!=a["ny"]: # box size mismatch. simply clip the box
				if not options.defocus:
					b=b.get_clip(Region((b["nx"]-a["ny"])//2, (b["ny"]-a["ny"])//2, a["ny"],a["ny"]))
				
			if type(info[-1])==str:
				initxf=eval(info[-1])
				if options.shrink>1:
					initxf.set_trans(initxf.get_trans()/float(options.shrink))
				
			if options.transonly:
				if options.smooth:
					m=options.maxshift
					trans=np.indices((m*2+1, m*2+1)).reshape((2,-1)).T-m
					scr=[]
					xf=Transform({"type":"eman","tx":initxf["tx"], "ty":initxf["ty"], "alt":initxf["alt"],"az":initxf["az"],"phi":initxf["phi"]})
					pj=a.project("standard", xf)
					
					for t in trans.tolist():
						pjts=pj.copy()
						pjts.translate(t[0], t[1],0)
						s=b.cmp("frc",pjts, {"minres":30, "maxres":5})
						scr.append(s)
				
				else:
					c=b.align("translational", pj, {"maxshift":options.maxshift})
					trans=c["xform.align2d"].get_trans()
					xf.translate(-trans)
					scr=c.cmp("frc",pj)
					
				r={"idx":ii,"xform.align3d":[xf], "score":scr}
			
			elif options.defocus:
				xf=Transform({"type":"eman","tx":initxf["tx"], "ty":initxf["ty"], "alt":initxf["alt"],"az":initxf["az"],"phi":initxf["phi"]})
				pj=a.project("standard", xf)
				pj=pj.get_clip(Region((pj["nx"]-b["nx"])//2, (pj["ny"]-b["ny"])//2, b["ny"],b["ny"]))
				
				ctf=b["ctf"]
				fft1=b.do_fft()
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				e=fft1.do_ift()
				
				fsc=e.calc_fourier_shell_correlation(pj)
				fsc=np.array(fsc).reshape((3,-1))[1]
				ds=1./(e["apix_x"]*e["ny"])
				zrg=np.arange(-.3,.3,0.01)
				scr=[]
				for dz in zrg:
					ctf1=EMAN2Ctf(ctf)
					ctf1.defocus=ctf1.defocus+dz
					c=ctf1.compute_1d(e["nx"]+2, ds, Ctf.CtfType.CTF_SIGN)
					c=fsc*c
					scr.append(np.sum(c[len(c)//3:len(c)*2//3]))
				
				r={"idx":ii,"xform.align3d":[xf], "score":scr}
				
				
			elif options.fromscratch:
				alignpm={"verbose":options.verbose,"sym":options.sym,"maxres":options.maxres}
				dic=[]
				
				b1=b.do_fft()
				b1.process_inplace("xform.phaseorigin.tocorner")
				c=b1.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)
				dic.append(c[0])
				
				bestmr=int(np.argmin([d["score"] for d in dic]))
				xf=dic[bestmr]["xform.align3d"]
				xf.set_mirror(bestmr)
				r={"idx":ii, "xform.align3d":[xf], "score":[dic[bestmr]["score"]]}
			
			else:
				nxf=options.refinentry
				astep=options.refineastep
				xfs=[]
				
				for i in range(nxf):
					d={"type":"eman","tx":initxf["tx"], "ty":initxf["ty"]}
					for ky in ["alt", "az", "phi"]:
						d[ky]=initxf[ky]+(i>0)*np.random.randn()*astep/np.pi*2
					xfs.append(Transform(d))
						
				alignpm={"verbose":options.verbose,"sym":options.sym,"maxshift":options.maxshift,"initxform":xfs, "maxang":astep*2.}
				#print("lenxfs:", len(xfs))
				if initxf["mirror"]:
					b=b.process("xform.flip", {"axis":'x'})
				b=b.do_fft()
				b.process_inplace("xform.phaseorigin.tocorner")
				c=b.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)

				xf=c[0]["xform.align3d"]
				xf.set_mirror(initxf["mirror"])
				r={"idx":ii,"xform.align3d":[xf], "score":[c[0]["score"]]}
			
			
			callback(100*float(infoi/len(self.data["info"])))
			rets.append(r)
			#print('time', time.time()-time0)
			#print(infoi,r)
		#callback(100)
			
		return rets

	
if __name__ == '__main__':
	main()

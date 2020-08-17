#!/usr/bin/env python
# Muyuan Chen 2018-04
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2jsondb import JSTask
from EMAN2PAR import EMTaskCustomer
#from scipy.optimize import minimize

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
	#parser.add_argument("--scipytest", action="store_true", default=False ,help="test scipy optimizer")
	parser.add_argument("--fromscratch", action="store_true", default=False ,help="align from scratch and ignore previous particle transforms. for spt mostly. will include mirror")
	parser.add_argument("--refineastep", type=float,help="Mean angular variation for refine alignment", default=2.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=4)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=8)
	parser.add_argument("--localrefine", action="store_true", default=False ,help="local refinement")
	parser.add_argument("--defocus", action="store_true", default=False ,help="refine defocus. Still under development")
	parser.add_argument("--seedmap", action="store_true", default=False ,help="seed")
	parser.add_argument("--ctfweight", action="store_true", default=False ,help="weight by ctf")
	parser.add_argument("--skipm3d", action="store_true", default=False ,help="skip make3d. only output aligned list")


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
	etc=EMTaskCustomer(options.parallel, module="e2spt_tiltrefine_oneiter.SptTltRefineTask")
		
	num_cpus = etc.cpu_est()
	
	print("{} total CPUs available".format(num_cpus))
	print("{} jobs".format(nptcl))
	infos=[[] for i in range(num_cpus)]
	for i,info in enumerate(pinfo):
		infos[i%num_cpus].append([i, info])
	
	tids=[]
	for info in infos:
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
		return
	
	maxl=np.max([len(s) for s in allscr])
	maxv=np.max(np.concatenate(allscr))
	for s in allscr:
		s.extend([maxv]*(maxl-len(s)))
	allscr=np.array(allscr)
	
	#print(np.min(allscr), np.mean(allscr), np.max(allscr), np.std(allscr))
	if options.skipm3d:
		pass
	else:
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
	
	
	cmd="e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads {trd} --outsize {bx} --apix {apx} --mode trilinear --keep {kp} --sym {sm} {seed} {par}".format(
		inp=lname, 
		out=threedout,
		bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=options.sym, trd=options.threads,par=m3dpar, seed=seed)
	
	if options.skipm3d:
		print("Skipping 3D reconstruction")
	else:
		run(cmd)

	E2end(logid)


def run(cmd):
	print(cmd)
	launch_childprocess(cmd)



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
				#if not options.defocus:
				b=b.get_clip(Region((b["nx"]-a["ny"])//2, (b["ny"]-a["ny"])//2, a["ny"],a["ny"]))
				
			if type(info[-1])==str:
				initxf=eval(info[-1])
				if options.shrink>1:
					initxf.set_trans(initxf.get_trans()/float(options.shrink))
					
				xf=Transform({"type":"eman","tx":initxf["tx"], "ty":initxf["ty"], "alt":initxf["alt"],"az":initxf["az"],"phi":initxf["phi"]})
				
			if options.transonly:
				#if options.smooth>0:
				m=options.maxshift
				trans=np.indices((m*2+1, m*2+1)).reshape((2,-1)).T-m
				scr=[]
				
				pj=a.project("standard", xf)
				
				for t in trans.tolist():
					pjts=pj.copy()
					pjts.translate(t[0], t[1],0)
					s=b.cmp("frc",pjts, {"minres":options.minres, "maxres":options.maxres})
					scr.append(s)
				
				
				if options.smooth<=0:
					s=np.argmin(scr)
					tx=trans[s].tolist()
					xf.translate(tx)
					scr=[np.min(scr)]
					
				r={"idx":ii,"xform.align3d":[xf], "score":scr}
			
			elif options.defocus:
				pj=a.project("standard", xf)
				#pj=pj.get_clip(Region((pj["nx"]-b["nx"])//2, (pj["ny"]-b["ny"])//2, b["ny"],b["ny"]))
				
				ctf=b["ctf"]
				fft1=b.do_fft()
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				e=fft1.do_ift()
				
				fsc=e.calc_fourier_shell_correlation(pj)
				fsc=np.array(fsc).reshape((3,-1))[1]
				ds=1./(e["apix_x"]*e["ny"])
				zrg=np.arange(-.5,.5,0.001)
				scr=[]
				for dz in zrg:
					ctf1=EMAN2Ctf(ctf)
					ctf1.defocus=ctf1.defocus+dz
					c=ctf1.compute_1d(e["nx"]+2, ds, Ctf.CtfType.CTF_SIGN)
					c=fsc*c
					scr.append(np.sum(c[len(c)//4:len(c)*2//3]))
				
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

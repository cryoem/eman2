#!/usr/bin/env python
# Muyuan Chen 2018-04
from past.utils import old_div
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
	parser.add_argument("--fromscratch", action="store_true", default=False ,help="align from scratch and ignore previous particle transforms. for spt mostly. will include mirror")
	parser.add_argument("--refineastep", type=float,help="Mean angular variation for refine alignment", default=2.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=4)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=8)


	parser.add_argument("--padby", type=float,default=2.0, help="pad by factor. default is 2")
	parser.add_argument("--maxres", type=float,default=-1, help="max resolution for cmp")
	parser.add_argument("--sym", type=str,help="symmetry. will use symmetry from spt refinement by default", default="c1")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

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
	if options.maxres>0:
		options.shrink=max(1, int(options.maxres/apix*.3))
		options.shrink=min(options.shrink, bxsz//48)
		print("Will shrink by {} and filter to {:.0f} A. Box size {}".format(options.shrink, options.maxres, bxsz//options.shrink))
	else:
		options.shrink=1
	#m.process_inplace('normalize.edgemean')
	
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
		tid=etc.send_task(task)
		tids.append(tid)
	
	while 1:
		st_vals = etc.check_task(tids)
		#print(st_vals)
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
	
	allscr=np.array([d["score"] for d in dics])
	print(np.min(allscr), np.mean(allscr), np.max(allscr), np.std(allscr))
	allscr=2-allscr
	s=allscr.copy()
	s-=np.mean(s)
	s/=np.std(s)
	clp=2
	ol=abs(s)>clp
	print("Removing {} outliers from {} particles..".format(np.sum(ol), len(s)))
	s=(s+clp)/clp/2
	s[ol]=0
	allscr=s
	allscr-=np.min(allscr)-1e-5
	allscr/=np.max(allscr)
	
	try: os.remove(lname)
	except: pass
	lout=LSXFile(lname, False)
	for i, dc in enumerate(dics):
		d=dc["xform.align3d"].get_params("eman")
		d["score"]=float(allscr[i])
		l=pinfo[i]
		lout.write(-1, l[0], l[1], str(d))

	lout=None

	pb=options.padby
	
	if options.parallel.startswith("mpi") and len(dics)>10000:
		m3dpar="--parallel {}".format(options.parallel)
	else:
		m3dpar=""
	cmd="e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads {trd} --outsize {bx} --apix {apx} --mode gauss_2 --keep {kp} --sym {sm} {par}".format(
		inp=lname, 
		out=threedout,
		bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=options.sym, trd=options.threads,par=m3dpar)
	
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
		
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		a=EMData(data["ref"],0)
		if options.maxres>0 and options.shrink>1:
			a.process_inplace("math.meanshrink",{"n":options.shrink})
		
		for i, infos in enumerate(data["info"]):
			ii=infos[0]
			info=infos[1]
			
			b=EMData(info[1],info[0])
				
			if options.maxres>0:
				if options.shrink>1:
					b.process_inplace("math.meanshrink",{"n":options.shrink})
				b.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1./options.maxres})
			
			if b["ny"]!=a["ny"]: # box size mismatch. simply clip the box
				b=b.get_clip(Region((b["nx"]-a["ny"])/2, (b["ny"]-a["ny"])/2, a["ny"],a["ny"]))
				
			if type(info[-1])==str:
				initxf=eval(info[-1])
				if options.shrink>1:
					initxf.set_trans(initxf.get_trans()/float(options.shrink))
				
			if options.transonly:
				
				xf=Transform({"type":"eman","tx":initxf["tx"], "ty":initxf["ty"], "alt":initxf["alt"],"az":initxf["az"],"phi":initxf["phi"]})
				pj=a.project("standard", xf)
				c=b.align("translational", pj, {"intonly":0, "maxshift":options.maxshift})
				trans=c["xform.align2d"].get_trans()
				xf.translate(-trans)
				scr=c.cmp("frc",pj)
				r={"idx":ii,"xform.align3d":xf, "score":scr}
				#print(ii, (initxf["tx"], initxf["ty"]), trans)
			
			elif options.fromscratch:
				alignpm={"verbose":0,"sym":options.sym}
				mriter=[False, True]
				dic=[]
				
				for mirror in mriter:
					
					if mirror:
						b1=b.process("xform.flip", {"axis":'x'})
					else:
						b1=b.copy()
						
					b1=b1.do_fft()
					b1.process_inplace("xform.phaseorigin.tocorner")
					c=b1.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)
					dic.append(c[0])
				
				bestmr=int(np.argmin([d["score"] for d in dic]))
				xf=dic[bestmr]["xform.align3d"]
				xf.set_mirror(bestmr)
				r={"idx":ii, "xform.align3d":xf, "score":dic[bestmr]["score"]}
				#print(ii, xf)
	
			else:
				nxf=options.refinentry
				astep=options.refineastep
				xfs=[]
				
				for i in range(nxf):
					d={"type":"eman","tx":initxf["tx"], "ty":initxf["ty"]}
					for ky in ["alt", "az", "phi"]:
						d[ky]=initxf[ky]+(i>0)*np.random.randn()*astep/np.pi*2
					xfs.append(Transform(d))
						
				alignpm={"verbose":0,"sym":options.sym,"maxshift":options.maxshift,"initxform":xfs, "maxang":astep*2.}
				#print("lenxfs:", len(xfs))
			
				b=b.do_fft()
				b.process_inplace("xform.phaseorigin.tocorner")
				c=b.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)

				xf=c[0]["xform.align3d"]
				r={"idx":ii,"xform.align3d":xf, "score":c[0]["score"]}
			
			#print(ii,info, r)
			if options.shrink>1:
				xf=r["xform.align3d"]
				xf.set_trans(options.shrink*xf.get_trans())
				r["xform.align3d"]=xf
			callback(100*float(i/len(self.data["info"])))
			rets.append(r)
		#callback(100)
			
		return rets

	
if __name__ == '__main__':
	main()

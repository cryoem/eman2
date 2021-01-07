#!/usr/bin/env python
# Muyuan Chen 2020-08
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2jsondb import JSTask
from EMAN2PAR import EMTaskCustomer

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ptclin", type=str,help="particle input", default=None)
	parser.add_argument("--output", type=str,help="score output", default=None)

	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12")

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data (threads * 8 particles)")

	parser.add_argument("--maxres", type=float,default=15, help="max resolution for cmp")
	parser.add_argument("--minres", type=float,default=300, help="min resolution for cmp")
	#parser.add_argument("--sym", type=str,help="symmetry. ", default="c1")
	parser.add_argument("--mask", type=str,help="mask. ", default=None)
	#parser.add_argument("--maxshift", type=int,help="max shift.", default=0)
	
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)


	threedname=args
	lst=LSXFile(options.ptclin, True)
	pinfo=[]
	nptcl=lst.n
	for i in range(nptcl):
		pinfo.append(lst.read(i))
	lst=None
	
	e=EMData(options.ptclin,0,True)
	nx=e["nx"]
	apix=e["apix_x"]
	b=good_size(nx*apix/options.maxres*2)
	options.shrink=max(1, nx/b)
	
	print("Initializing parallelism...")
	etc=EMTaskCustomer(options.parallel, module="e2spa_classify.SpaClassifyTask")
		
	num_cpus = etc.cpu_est()
	
	print("{} total CPUs available".format(num_cpus))
	print("{} jobs".format(nptcl))
	infos=[[] for i in range(num_cpus)]
	for i,info in enumerate(pinfo):
		infos[i%num_cpus].append([i, info])
	
	tids=[]
	for info in infos:
		task = SpaClassifyTask(info, threedname, options)
			
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
	
	dics=np.zeros((nptcl, len(threedname)))
	for i in tids:
		ret=etc.get_results(i)[1]
		for r in ret:
			#print(r)
			ii=r.pop("idx")
			dics[ii]=r["score"]
	
	del etc
	np.savetxt(options.output, dics)

	E2end(logid)


def run(cmd):
	print(cmd)
	launch_childprocess(cmd)



class SpaClassifyTask(JSTask):
	
	
	def __init__(self, info, ref, options):
		
		data={"info":info, "ref": ref}
		JSTask.__init__(self,"SpaClassifyTask",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		refnames=data["ref"]
		
		ref=EMData(refnames[0],0,True)
		ny=ref["ny"]
		apix=ref["apix_x"]
	
		pmax=ceil(ny*apix/options.maxres)
		ss=good_size(pmax*3)
		ss=int(min(ss, ny))
	
		pmin=ceil(ny*apix/options.minres)
		shrink=ny/ss
		
			
		refs=[EMData(r) for r in refnames]
		
		if options.mask:
			mask=EMData(options.mask)
			avgr=Averagers.get("mean")
			for r in refs:
				avgr.add_image(r)
			avg=avgr.finish()
		else:
			mask=None
		
		for i,ref in enumerate(refs):
			if mask:
				ref=ref*mask + avg*(1-mask)
			
			ref=ref.do_fft()
			ref.process_inplace("xform.phaseorigin.tocenter")
			ref.process_inplace("xform.fourierorigin.tocenter")
			ref=ref.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
			refs[i]=ref
		
		for infoi, infos in enumerate(data["info"]):
			ii=infos[0]
			info=infos[1]
			img=EMData(info[1],info[0])
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			img.process_inplace("xform.fourierorigin.tocenter")
			img=img.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
			img.process_inplace("xform.fourierorigin.tocenter")
				
			initxf=eval(info[-1])
				
			xf=Transform({"type":"eman","tx":initxf["tx"]/shrink, "ty":initxf["ty"]/shrink, "alt":initxf["alt"],"az":initxf["az"],"phi":initxf["phi"]})
				
			#trans=np.indices((m*2+1, m*2+1)).reshape((2,-1)).T-m
			scr=[]
			for ref in refs:
				pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
				s=img.cmp("frc",pj, {"pmin":pmin, "maxres":pmax})
				scr.append(s)
			
			#print(infoi, scr)
			
			
			#scr=[np.min(scr)]
				
			r={"idx":ii, "score":scr}
			
			
			callback(100*float(infoi/len(data["info"])))
			rets.append(r)
			
		return rets

	
if __name__ == '__main__':
	main()

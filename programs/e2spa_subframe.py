#!/usr/bin/env python
# Muyuan Chen 2021-09
from EMAN2 import *
import numpy as np
import queue
import threading
from EMAN2jsondb import JSTask
from EMAN2PAR import EMTaskCustomer
from scipy.optimize import minimize

def main():
	
	usage="Sub-frame alignment for single particle analysis. Need specific format for motion correction and refinement to work properly. The improvement in resolution seems marginal. Not ready to use yet."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ptclin", type=str,help="particle input", default=None)
	parser.add_argument("--ptclout", type=str,help="particle output", default=None)
	parser.add_argument("--ref", type=str,help="reference input", default=None)
	#parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.8", default=0.8)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12")

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data")
	#parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=-1)
	#parser.add_argument("--localrefine", type=int, default=-1 ,help="local refinement. larger value correspond to smaller local region")
	#parser.add_argument("--goldcontinue", action="store_true", default=False ,help="split even/odd subset and references.")
	#parser.add_argument("--ctfweight", action="store_true", default=False ,help="weight by ctf. not used yet...")
	#parser.add_argument("--slow", action="store_true", default=False ,help="slow but finer search")
	parser.add_argument("--maxres", type=float,default=-1, help="max resolution for cmp")
	parser.add_argument("--minrespx", type=int,default=4, help="skip the first x pixel in fourier space")
	#parser.add_argument("--sym", type=str,help="symmetry. ", default="c1")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	#parser.add_argument("--nkeep", type=int,help="", default=1)
	parser.add_argument("--verbose","-v", type=int,help="Verbose", default=0)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#m=EMData(options.ref)
	#bxsz=m["nx"]
	#apix=m["apix_x"]
	
	options.shrink=1
	pinfo=load_lst_params(options.ptclin)
	nptcl=len(pinfo)
	
	print("Initializing parallelism...")
	etc=EMTaskCustomer(options.parallel, module="e2spa_subframe.SpaAlignTask")	
	num_cpus = etc.cpu_est()
	
	print("{} particles".format(nptcl))
	print("{} total CPUs available".format(num_cpus))
	if options.debug: 
		nptcl=min(4*num_cpus, nptcl)
		print("Debugging mode. running on one thread with 8 particles")
		
	nbatch=min(nptcl//4, num_cpus)
	
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
		
		def test_trans(tx):
			sft=np.linalg.norm(tx)
			if sft>mxsft: return 1
			pjsft=[p.process("xform", {"tx":tx[0], "ty":tx[1]}) for p in projs]

			fsc=[m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmall, pjsft)]
			fsc=np.array(fsc).reshape((len(imgsmall), 3,-1))[:,1]
			scr=np.mean(fsc[:, options.minrespx:maxrescut], axis=1)
			scr=np.sum(scr*wt)/np.sum(wt)
			return -scr
		
		options=self.options
		data=self.data
		callback(0)
		rets=[]
		
		ptclinfo=load_lst_params(options.ptclin)
		
		#### do the even/odd split if required
		reffile=options.ref
		#if options.goldcontinue:
		refnames=[reffile[:-4]+f"_{eo}.hdf" for eo in ["even", "odd"]]
		
		ref=EMData(refnames[0], 0, True)
		ny0=ny=ref["ny"]
		if options.maxres>0:
			maxrescut=ceil(ny*ref["apix_x"]/options.maxres)
			maxy=good_size(maxrescut*3)
			maxy=int(min(maxy, ny))
		else:
			maxy=ny
			maxrescut=ny//2
			
		ss=maxy
		
		refs=[]
		for r in refnames:
			ref=EMData(r)
			ref=ref.do_fft()
			ref.process_inplace("xform.phaseorigin.tocenter")
			ref.process_inplace("xform.fourierorigin.tocenter")
			ref=ref.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
			refs.append(ref)
		
		
		for infoi, infos in enumerate(data["info"]):
			ii=infos[0]
			info=infos[1]
			
			#print(info)
			pinfo=[p for p in ptclinfo if p["src"]==info["src"]]
			pinfo=[p for p in pinfo if p["frameid"]==info["frameid"]]
			if len(pinfo)<5:
				rets.append((ii, info))
				continue
			
			ptcls=[EMData(p["src"], p["idx"], True) for p in pinfo]
			coord=np.array([p["ptcl_source_coord"] for p in ptcls])
			p0=EMData(info["src"],info["idx"], True)
			c=p0["ptcl_source_coord"]
			dst=np.linalg.norm(coord-c, axis=1)
			smoothN=7
			smooth=100
			srtid=np.argsort(dst)[:smoothN]
			dst=dst[srtid]
			wt=np.exp(-dst/smooth)
			
			imgsmall=[]
			initxfs=[]
			projs=[]
			for i in srtid:
				img=EMData(pinfo[i]["src"], pinfo[i]["idx"])
				img.clip_inplace(Region((img["nx"]-ny)//2, (img["ny"]-ny)//2, ny, ny))
				img=img.do_fft()
				img.process_inplace("xform.phaseorigin.tocenter")
				img.process_inplace("xform.fourierorigin.tocenter")
				img=img.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
				img.process_inplace("xform.fourierorigin.tocenter")
				imgsmall.append(img)
				
				xf=Transform(pinfo[i]["xform.projection"])
				xf.set_trans(xf.get_trans()*ss/ny)
				initxfs.append(xf)
				rid=pinfo[i]["class"]
				ref=refs[rid]
				pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
				projs.append(pj)
				
				
			mxsft=ss//8
			r0=test_trans([0,0])
			res=minimize(test_trans, [0,0],
				method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":10})

			newxf=Transform(info["xform.projection"])
			x=res.x*ny/ss
			newxf.translate(x.tolist())
			
			if options.debug:
				print(res.x, r0, res.fun)
				
				
			r={	"src":info["src"], "idx":info["idx"], "frameid":info["frameid"],
				"xform.projection":newxf, "score":float(res.fun),"class":info["class"]}
			
			callback(100*float(infoi/len(self.data["info"])))
			rets.append((ii, r))
			
		return rets




if __name__ == '__main__':
	main()

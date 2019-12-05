#!/usr/bin/env python
from __future__ import division
# Muyuan Chen 2018-07
from EMAN2 import *
import numpy as np
import threading
import Queue

def alifn(jsd,ii, info,a,options):
	#### deal with mirror
	dic=[]
	alignpm={"verbose":0,"sym":options.sym}
	if options.maxshift>=0:
		alignpm["maxshift"]=options.maxshift
	
	#print(ii, info)
	if len(options.refinefrom)>0:
		
		nxf=options.refinentry
		astep=options.refineastep
		xfs=[]
		initxf=eval(info[-1])
		
		for i in range(nxf):
			d={"type":"eman","tx":initxf["tx"], "ty":initxf["ty"]}
			for ky in ["alt", "az", "phi"]:
				d[ky]=initxf[ky]+(i>0)*np.random.randn()*astep
				xfs.append(Transform(d))
				
		alignpm["initxform"]=xfs
		if options.maxshift<0:
			alignpm["maxshift"]=10
	
	if options.incmirror==1:
		mriter=[False, True]
	else:
		mriter=[False]
		
	for mirror in mriter:
		
		b=EMData(info[1],info[0])
		if b["ny"]!=a["ny"]: # box size mismatch. simply clip the box
			b=b.get_clip(Region((b["nx"]-a["ny"])/2, (b["ny"]-a["ny"])/2, a["ny"],a["ny"]))
			
		if mirror:
			b.process_inplace("xform.flip", {"axis":'x'})
			
		b=b.do_fft()
		b.process_inplace("xform.phaseorigin.tocorner")
		#print(alignpm)
		c=b.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)
		dic.append(c[0])
		#print(mirror,c)
		
	bestmr=np.argmin([d["score"] for d in dic])
	xf=dic[bestmr]["xform.align3d"]
	xf.set_mirror(bestmr)
	c={"xform.align3d":xf, "score":dic[bestmr]["score"]}
	
	jsd.put((ii,c))
	
	#if options.verbose>1 : print("{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"]))


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--ref", type=str,help="reference", default="")
	parser.add_argument("--output", type=str,help="output", default="")
	parser.add_argument("--refinefrom", type=str,help="refine alignment using transform from a lst file", default="")
	parser.add_argument("--refineastep", type=float,help="angular step for refine alignment (gauss std)", default=10.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=32)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=-1)
	parser.add_argument("--threads", type=int,help="threads", default=12)
	parser.add_argument("--incmirror", type=int,help="include mirror", default=1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	if options.output=="" or options.ref=="":
		print("No output or reference. exit")
		return
	
	ptcls=args[0]
	ref=options.ref
	
	nptcl=EMUtil.get_image_count(ptcls)
	if ptcls.endswith("lst"):
		#### parse list
		pinfo=[]
		lst=LSXFile(ptcls)
		nptcl=lst.n
		for i in range(nptcl):
			pinfo.append(lst.read(i))
		lst=None
	else:
		pinfo=[(i, ptcls) for i in range(nptcl)]
	
	if len(options.refinefrom)>0:
		lst=LSXFile(options.refinefrom)
		if (lst.n!=nptcl):
			print("particle number in refinefrom file does not match input...")
			return
		
		for i in range(nptcl):
			l=lst.read(i)
			if (pinfo[i][0]==l[0] and pinfo[i][1]==l[1]):
				pinfo[i].append(l[2])
			else:
				print("refinefrom file and particle input mismatch at line {}..".format(i))
				print(pinfo[i])
				print(l)
				return
		lst=None
		
	
	
	e=EMData(ref)
	a=e.do_fft()
	a.process_inplace("xform.phaseorigin.tocorner")
	
	jsd=Queue.Queue(0)
	thrds=[threading.Thread(target=alifn,args=([jsd, i, info, a, options])) for i,info in enumerate(pinfo)]
	#thrds=thrds[:2]
	thrtolaunch=0
	tsleep=threading.active_count()
	
	print("starting threads...")
	ndone=0
	dics=[0]*nptcl
	nthrd=options.threads
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
		if thrtolaunch<len(thrds):
			while (threading.active_count()==nthrd+tsleep ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(.2)


		while not jsd.empty():
			ii, dc=jsd.get()
			dics[ii]=dc
			#print(dc)
			ndone+=1
			#if ndone%1000==0:
			print("\t{}/{} finished.".format(ndone, nptcl), end='\r')
			sys.stdout.flush()

	
	#### weight by score
	allscr=np.array([d["score"] for d in dics])
	print("min\tmax\tmean\tstd")
	print(np.min(allscr), np.max(allscr),  np.mean(allscr), np.std(allscr))
	allscr*=-1
	s=allscr.copy()
	s-=np.mean(s)
	s/=np.std(s)
	clp=2
	ol=abs(s)>clp
	print("Removing {} outliers from {} particles..".format(np.sum(ol), len(s)))
	s=(s+clp)/clp/2
	s[ol]=0
	allscr=s
	
	#allscr-=np.min(allscr)-1e-5
	#allscr/=np.max(allscr)

	lname=options.output
	try: os.remove(lname)
	except: pass
	lout=LSXFile(lname, False)
	for i, dc in enumerate(dics):
		d=dc["xform.align3d"].get_params("eman")
		d["score"]=float(allscr[i])
		l=pinfo[i]
		lout.write(-1, l[0], l[1], str(d))

	lout.close()

	#os.system("rm tmp*.hdf")
	#e=EMData("emd5592_shrink.hdf")
	#xf=Transform()
	#pj=e.project('standard', xf)
	##plt.imshow(pj.numpy())
	
	#ali=pj.align('rotate_translate_2d_to_3d_tree', e, {"verbose":10})
	#print(ali["xform.align3d"])
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

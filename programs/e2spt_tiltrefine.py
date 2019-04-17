#!/usr/bin/env python
# Muyuan Chen 2018-04
from __future__ import print_function
from __future__ import division
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

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path", type=str,help="path", default=None, guitype='strbox',row=0, col=0,rowspan=1, colspan=1)
	parser.add_argument("--iter", type=int,help="start from iteration X", default=-1, guitype='intbox',row=0, col=1,rowspan=1, colspan=1)
	parser.add_argument("--niters", type=int,help="run this many iterations. Default is 4.", default=4, guitype='intbox',row=0, col=2,rowspan=1, colspan=1)

	parser.add_argument("--sym", type=str,help="symmetry. will use symmetry from spt refinement by default", default="", guitype='strbox',row=2, col=0,rowspan=1, colspan=1)
	parser.add_argument("--padby", type=float,help="pad by factor. default is 2", default=2., guitype='floatbox',row=1, col=1,rowspan=1, colspan=1)
	parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.5", default=0.5, guitype='floatbox',row=1, col=2,rowspan=1, colspan=1)

	parser.add_argument("--maxalt", type=float,help="max altitude to insert to volume", default=90.0, guitype='floatbox',row=1, col=0,rowspan=1, colspan=1)	
	parser.add_argument("--nogs", action="store_true", default=False ,help="skip gold standard...", guitype='boolbox',row=2, col=1,rowspan=1, colspan=1)
	#parser.add_argument("--localfilter", type=int, default=-1 ,help="use tophat local. specify 0 or 1 to overwrite the setting in the spt refinement")
	parser.add_argument("--tophat", type=str, default="auto" ,help="auto: same as spt refinement; local; global;")
	parser.add_argument("--mask", type=str, default="None" ,help="Refinement masking. default is the same as the spt refinement. Use Auto for automasking, None for keeping the same masking as spt refinement",guitype='strbox',row=3, col=0,rowspan=1, colspan=2)	

	parser.add_argument("--threads", type=int,help="Number of CPU threads to use. Default is 12.", default=12, guitype='intbox',row=2, col=2,rowspan=1, colspan=1)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12", guitype='strbox',row=4, col=0,rowspan=1, colspan=3)

	parser.add_argument("--refineastep", type=float,help="angular variation for refine alignment (gauss std)", default=8.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=32)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=10)

	parser.add_argument("--buildsetonly", action="store_true", default=False ,help="build sets only")
	parser.add_argument("--output", type=str,help="Write results to this directory. We do not recommend changing this.", default="subtlt")#, guitype='intbox',row=2, col=1,rowspan=1, colspan=1)

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data (threads * 8 particles)")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	itr = options.iter

	oldpath = options.path	
	if not oldpath:
		print("No input path. Exit.")
		return
	

	if options.iter != -1:
		itr = options.iter
	elif "spt" in oldpath:
		for f in sorted(os.listdir(oldpath)):
			if "particle_parms" in f:
				itrstr = f[15:].split(".")[0]
				if os.path.isfile(os.path.join(oldpath,"threed_{}.hdf".format(itrstr))):
					itr = int(itrstr)
	else:
		for f in sorted(os.listdir(oldpath)):
			if re.match("threed_[0-9][0-9].hdf",f):
				itr = int(f[7:].split(".")[0])

	# print(oldpath)
	fromspt=True
	if "0_subtlt_params.json" in os.listdir(oldpath):
		print("Continuing from a subtilt refinement...")
		fromspt=False
		
	
	path = make_path(options.output)
	
	if not os.path.isfile(os.path.join(oldpath,"threed_{:02d}.hdf".format(itr))):
		print("Could not locate {}/threed_{:02d}.hdf".format(oldpath,itr))
		print("Please specify the iteration number (--iter) of a completed subtomogram refinement.")
		sys.exit(1)
	#elif not os.path.isfile("{}/particle_parms_{:02d}.json".format(oldpath,itr)):
		#print("Could not locate {}/particle_parms_{:02d}.json".format(oldpath,itr))
		#print("Please specify the iteration number (--iter) of a completed subtomogram refinement.")
		#sys.exit(1)
	else:
		#copy2("{}/0_spt_params.json".format(oldpath),"{}/0_subtlt_params.json".format(path))

		oldmap = os.path.join(oldpath,"threed_{:02d}.hdf".format(itr))
		oem = os.path.join(oldpath,"threed_{:02d}_even.hdf".format(itr))
		oom = os.path.join(oldpath,"threed_{:02d}_odd.hdf".format(itr))
		oldfsc = os.path.join(oldpath, "fsc_masked_{:02d}.txt".format(itr))

		copy2(oldmap,os.path.join(path,"threed_00.hdf"))
		copy2(oldfsc, os.path.join(path, "fsc_masked_00.txt"))
		copy2(oem,os.path.join(path,"threed_00_even.hdf"))
		copy2(oom,os.path.join(path,"threed_00_odd.hdf"))
		
		if fromspt:
			oldparm = os.path.join(oldpath,"particle_parms_{:02d}.json".format(itr))
			copy2(oldparm,os.path.join(path,"particle_parms_00.json"))
		else:
			for eo in ["even", "odd"]:
				oali = os.path.join(oldpath,"ali_ptcls_{:02d}_{}.lst".format(itr, eo))
				copy2(oali,os.path.join(path,"ali_ptcls_00_{}.lst".format(eo)))


	e=EMData(os.path.join(path,"threed_00.hdf"))
	
	bxsz=e["nx"]
	apix=e["apix_x"]
	jd = js_open_dict(os.path.join(path, "0_subtlt_params.json"))
	jd.update(vars(options))
	jd["cmd"] = " ".join(sys.argv)
	jd["path"] = oldpath
	jd["iter"] = itr
	jd["output"] = path
	
	options.ptclkeep=1.0

	if fromspt:
		sptparms = os.path.join(oldpath,"0_spt_params.json")
	else:
		sptparms = os.path.join(oldpath,"0_subtlt_params.json")
	if os.path.isfile(sptparms):
		oldjd = js_open_dict(sptparms)
		#print(oldjd.keys())
		jd["mass"] = oldjd["mass"]
		jd["setsf"] = oldjd["setsf"]
		jd["sym"] = oldjd["sym"]
		jd["localfilter"]=oldjd["localfilter"]
		jd["mask"]=oldjd["mask"]
		if oldjd.has_key("radref"):
			jd["radref"]=oldjd["radref"]
			
		if fromspt:
			options.ptclkeep=oldjd["pkeep"]
			
		oldjd.close()
	else:
		print("Cannot find {}. exit.".format(sptparms))
	
	if options.mask.lower()!="none":
		print("Overwritting masking")
		if options.mask.lower()=="auto":
			jd["mask"]=""
		else:
			jd["mask"]=options.mask
	
	
	#if options.localfilter==0:
		#jd["localfilter"]=False
	#elif options.localfilter==1:
		#jd["localfilter"]=True
		
	if len(options.sym)>0:
		jd["sym"]=options.sym
		
	jsparams=jd.data
	jd.close()
	jd = jsparams

	if fromspt:
		js=js_open_dict(os.path.join(path,"particle_parms_00.json"))
		k=list(js.keys())[0]
		src=eval(k)[0]
		
		print("loading 3D particles from {}".format(base_name(src)))
		print("box size {}, apix {:.2f}".format(bxsz, apix))

		lname=[os.path.join(path, "ali_ptcls_00_{}.lst".format(eo)) for eo in ["even", "odd"]]
		for l in lname:
			try: os.remove(l)
			except:pass
		
		lst=[LSXFile(m, False) for m in lname]
		n3d=len(list(js.keys()))
		if options.ptclkeep<1.0:
			score=[]
			for k in list(js.keys()):
				score.append(float(js[k]["score"]))
		
			simthr=np.sort(score)[int(len(score)*options.ptclkeep)]
			print("removing bad particles...")
			
		else:
			simthr=10000
			
		for ky in js.keys():
			
			src,ii=eval(ky)
			e=EMData(src, ii, True)
			fname=e["class_ptcl_src"]
			ids=e["class_ptcl_idxs"]
			#ky="('{}', {})".format(src, ii)
			dic=js[ky]
			xali=dic["xform.align3d"]
			scr=float(dic["score"])
			if scr>simthr:
				continue
			
			if "__even" in src:
				eo=0
			elif "__odd" in src:
				eo=1
			else:
				eo=ii%2
			
			#print(src, eo)
			
			for i in ids:
				try:
					m=EMData(fname, i, True)
				except:
					continue
				xf=m["xform.projection"]
				dc=xf.get_params("xyz")
				if abs(dc["ytilt"])>options.maxalt:
					continue
				rot=xf*xali.inverse()
				lst[eo].write(-1, i, fname, str(rot.get_params("eman")))
		for l in lst:
			l.close()
		js=None

	if options.buildsetonly: return

	from EMAN2PAR import EMTaskCustomer
	
	for itr in range(0,options.niters):

		

		for eo in ["even", "odd"]:
			
			if options.nogs:
				threedname=os.path.join(path, "threed_{:02d}.hdf".format(itr))
			else:
				threedname=os.path.join(path, "threed_{:02d}_{}.hdf".format(itr, eo))
			
			lstname=os.path.join(path, "ali_ptcls_{:02d}_{}.lst".format(itr, eo))
			lst=LSXFile(lstname, True)
			m=EMData(threedname)
			
			m.process_inplace('normalize.edgemean')
			
			pinfo=[]
			if options.debug: nptcl=options.threads*8
			else: nptcl=lst.n
			for i in range(nptcl):
				pinfo.append(lst.read(i))
			lst=None
			
			print("Initializing parallelism...")
			etc=EMTaskCustomer(options.parallel, module="e2spt_tiltrefine.SptTltRefineTask")
			num_cpus = etc.cpu_est()
			
			print("{} total CPUs available".format(num_cpus))
			print("{} jobs".format(nptcl))
			
			infos=[[] for i in range(num_cpus)]
			for i,info in enumerate(pinfo):
				infos[i%num_cpus].append([i, info])
			
			tids=[]
			for info in infos:
				task = SptTltRefineTask(info, m, options)
				tid=etc.send_task(task)
				tids.append(tid)
			
			while 1:
				st_vals = etc.check_task(tids)
				#print("{:.1f}/{} finished".format(np.mean(st_vals), 100))
				#print(tids)
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
			allscr*=-1
			s=allscr.copy()
			s-=np.mean(s)
			s/=np.std(s)
			clp=2
			ol=abs(s)>clp
			print("Removing {} outliers from {} particles..".format(np.sum(ol), len(s)))
			s=old_div(old_div((s+clp),clp),2)
			s[ol]=0
			allscr=s
			#allscr-=np.min(allscr)-1e-5
			#allscr/=np.max(allscr)

			lname=os.path.join(path, "ali_ptcls_{:02d}_{}.lst".format(itr+1, eo))
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
			threedout=os.path.join(path, "threed_{:02d}_{}.hdf".format(itr+1, eo))
			
			if options.parallel.startswith("mpi"):
				m3dpar="--parallel {}".format(options.parallel)
			else:
				m3dpar=""
			cmd="e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads {trd} --outsize {bx} --apix {apx} --mode gauss_2 --keep {kp} --sym {sm} {par}".format(
				inp=lname, 
				out=threedout,
				bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=jd["sym"], trd=options.threads,par=m3dpar)
			
			run(cmd)
			run("e2proc3d.py {} {}".format(threedout, os.path.join(path, "threed_raw_{}.hdf".format(eo))))

		s = ""
		
		if jd.has_key("goldstandard"): 
			if jd["goldstandard"] > 0: 
				s += " --align"
		if jd.has_key("setsf"):
			s += " --setsf {}".format(jd['setsf']) #options.setsf)
		
		s+=" --tophat global"
		
		if options.tophat=="auto" and jd.has_key("localfilter"):
			s += " --tophat local"
		elif options.tophat=="local":
			s += " --tophat local"
		elif options.tophat=="global":
			s += " --tophat global"
			
		msk = jd["mask"] #{}/mask_tight.hdf".format(path)
		if len(msk)>0:
			if os.path.isfile(msk):
				msk=" --automask3d mask.fromfile:filename={}".format(msk)
			else:
				msk=" --automask3d {}".format(msk)
				
		s+=msk

		# get target resolution from last iteration map
		ref=os.path.join(path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(path, "fsc_masked_{:02d}.txt".format(itr)))
		rs=1./fsc[fsc[:,1]<0.3, 0][0]
		curres=rs*.5
		
		even=os.path.join(path, "threed_{:02d}_even.hdf".format(itr+1))
		odd=os.path.join(path, "threed_{:02d}_odd.hdf".format(itr+1))
		if jd.has_key("radref") and jd["radref"]!="":
			rfname=str(jd["radref"])
			print("doing radial correction with {}".format(rfname))
			for f in [even, odd]:
				e=EMData(f)
				rf=EMData(rfname).process("normalize")
				#rf=rf.align("translational",e)
				
				e.process_inplace("filter.matchto",{"to":rf})
				e.process_inplace("normalize")

				radrf=rf.calc_radial_dist(rf["nx"]//2, 0.0, 1.0,1)
				rade=e.calc_radial_dist(e["nx"]//2, 0.0, 1.0,1)
				rmlt=np.sqrt(np.array(radrf)/np.array(rade))
				#rmlt[0]=1
				rmlt=from_numpy(rmlt.copy())
				er=e.mult_radial(rmlt)
				er.write_image(f)
		

		#os.system("rm {}/mask*.hdf {}/*unmasked.hdf".format(path, path))
		ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --restarget {} --threads {} --sym {} --mass {} {}".format(even, odd, 
			os.path.join(path, "threed_{:02d}.hdf".format(itr+1)), itr+1, curres, options.threads, jd["sym"], jd["mass"], s)
		run(ppcmd)

		fsc=np.loadtxt(os.path.join(path, "fsc_masked_{:02d}.txt".format(itr+1)))
		rs=1./fsc[fsc[:,1]<0.3, 0][0]
		print("Resolution (FSC<0.3) is ~{:.1f} A".format(rs))
				
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
		for i, infos in enumerate(self.data["info"]):
			ii=infos[0]
			info=infos[1]
			#print(ii, info)
			ptcl=EMData(info[1],info[0])
			a=data["ref"]
			
			nxf=options.refinentry
			astep=options.refineastep
			xfs=[]
			initxf=eval(info[-1])
			for i in range(nxf):
				d={"type":"eman","tx":initxf["tx"], "ty":initxf["ty"]}
				for ky in ["alt", "az", "phi"]:
					d[ky]=initxf[ky]+(i>0)*np.random.randn()*astep
					xfs.append(Transform(d))
					
			alignpm={"verbose":0,"sym":options.sym,"maxshift":options.maxshift,"initxform":xfs}
		
			b=EMData(info[1],info[0])
			if b["ny"]!=a["ny"]: # box size mismatch. simply clip the box
				b=b.get_clip(Region((b["nx"]-a["ny"])/2, (b["ny"]-a["ny"])/2, a["ny"],a["ny"]))
				
				
			b=b.do_fft()
			b.process_inplace("xform.phaseorigin.tocorner")
			c=b.xform_align_nbest("rotate_translate_2d_to_3d_tree",a, alignpm, 1)

			xf=c[0]["xform.align3d"]
			r={"idx":ii,"xform.align3d":xf, "score":c[0]["score"]}
			
			#print(ii,info, r)
			callback(float(i/len(self.data["info"])))
			rets.append(r)
		callback(100)
			
		return rets

	
if __name__ == '__main__':
	main()

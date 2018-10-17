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

def alifn(jsd,ii, info,a,options):
	#### deal with mirror
		
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
	c={"xform.align3d":xf, "score":c[0]["score"]}
	
	jsd.put((ii,c))

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path", type=str,help="path", default=None, guitype='strbox',row=0, col=0,rowspan=1, colspan=1)
	parser.add_argument("--iter", type=int,help="start from iteration X", default=-1, guitype='intbox',row=0, col=1,rowspan=1, colspan=1)
	parser.add_argument("--niters", type=int,help="run this many iterations. Default is 4.", default=4, guitype='intbox',row=0, col=2,rowspan=1, colspan=1)

	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=1, col=0,rowspan=1, colspan=1)
	parser.add_argument("--padby", type=float,help="pad by factor. default is 2", default=2., guitype='floatbox',row=1, col=1,rowspan=1, colspan=1)
	parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.5", default=0.5, guitype='floatbox',row=1, col=2,rowspan=1, colspan=1)

	parser.add_argument("--maxalt", type=float,help="max altitude to insert to volume", default=90.0, guitype='floatbox',row=2, col=0,rowspan=1, colspan=1)	
	parser.add_argument("--nogs", action="store_true", default=False ,help="skip gold standard...", guitype='boolbox',row=2, col=1,rowspan=1, colspan=1)

	parser.add_argument("--buildsetonly", action="store_true", default=False ,help="build sets only")
	parser.add_argument("--output", type=str,help="Write results to this directory. We do not recommend changing this.", default="subtlt")#, guitype='intbox',row=2, col=1,rowspan=1, colspan=1)

	parser.add_argument("--refineastep", type=float,help="angular step for refine alignment (gauss std)", default=10.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=32)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=10)

	parser.add_argument("--threads", type=int,help="Number of CPU threads to use. Default is 12.", default=12, guitype='intbox',row=2, col=2,rowspan=1, colspan=1)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12", guitype='strbox',row=4, col=0,rowspan=1, colspan=3)

	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data (threads * 8 particles)")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	itr = options.iter

	oldpath = options.path	
	if not oldpath:
		print("No input path. Exit.")
		return
	
	path = make_path(options.output)

	if options.iter != -1:
		itr = options.iter
	elif "spt" in oldpath:
		for f in sorted(os.listdir(oldpath)):
			if "particle_parms" in f:
				itrstr = f[15:].split(".")[0]
				if os.path.isfile("{}/threed_{}.hdf".format(oldpath,itrstr)):
					itr = int(itrstr)
	else:
		for f in sorted(os.listdir(oldpath)):
			if re.match("threed_[0-9][0-9].hdf",f):
				itr = int(f[7:].split(".")[0])

	print(oldpath)

	if not os.path.isfile("{}/threed_{:02d}.hdf".format(oldpath,itr)):
		print("Could not locate {}/threed_{:02d}.hdf".format(oldpath,itr))
		print("Please specify the iteration number (--iter) of a completed subtomogram refinement.")
		sys.exit(1)
	elif not os.path.isfile("{}/particle_parms_{:02d}.json".format(oldpath,itr)):
		print("Could not locate {}/particle_parms_{:02d}.json".format(oldpath,itr))
		print("Please specify the iteration number (--iter) of a completed subtomogram refinement.")
		sys.exit(1)
	else:
		copy2("{}/0_spt_params.json".format(oldpath),"{}/0_subtlt_params.json".format(path))

		oldmap = os.path.join(oldpath,"threed_{:02d}.hdf".format(itr))
		oem = os.path.join(oldpath,"threed_{:02d}_even.hdf".format(itr))
		oom = os.path.join(oldpath,"threed_{:02d}_odd.hdf".format(itr))
		oldparm = os.path.join(oldpath,"particle_parms_{:02d}.json".format(itr))
		oldfsc = os.path.join(oldpath, "fsc_masked_{:02d}.txt".format(itr))

		copy2(oldmap,os.path.join(path,"threed_00.hdf"))
		copy2(oldparm,os.path.join(path,"particle_parms_00.json"))
		copy2(oldfsc, os.path.join(path, "fsc_masked_00.txt"))
		copy2(oem,os.path.join(path,"threed_00_even.hdf"))
		copy2(oom,os.path.join(path,"threed_00_odd.hdf"))

		oldmsk = "{}/mask_tight.hdf".format(oldpath)
		if os.path.isfile(oldmsk):
			msk = "{}/mask_tight.hdf".format(path)
			copy2(oldmsk, msk)

	jd = js_open_dict("{}/0_subtlt_params.json".format(path))
	jd["path"] = options.path
	jd["iter"] = itr
	jd["sym"] = options.sym
	jd["niters"] = options.niters
	jd["padby"] = options.padby
	jd["keep"] = options.keep
	jd["maxalt"] = options.maxalt
	jd["refineastep"] = options.refineastep
	jd["refinentry"] = options.refinentry
	jd["maxshift"] = options.maxshift
	jd["nogs"] = options.nogs
	jd["parallel"] = options.parallel
	jd["threads"] = options.threads
	jd["mask"] = msk
	jd.close()

	js=js_open_dict(os.path.join(path,"particle_parms_00.json"))
	k=list(js.keys())[0]
	src=eval(k)[0]
	e=EMData(src, 0,True)
	bxsz=e["nx"]
	apix=e["apix_x"]
	print("loading 3D particles from {}".format(src))
	print("box size {}, apix {:.2f}".format(bxsz, apix))

	lname=[os.path.join(path, "ali_ptcls_00_{}.lst".format(eo)) for eo in ["even", "odd"]]
	for l in lname:
		try: os.remove(l)
		except:pass
	lst=[LSXFile(m, False) for m in lname]
	n3d=len(list(js.keys()))
	for ii in range(n3d):
		e=EMData(src, ii, True)
		fname=e["class_ptcl_src"]
		ids=e["class_ptcl_idxs"]
		ky="('{}', {})".format(src, ii)
		dic=js[ky]
		xali=dic["xform.align3d"]
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
			lst[ii%2].write(-1, i, fname, str(rot.get_params("eman")))
	for l in lst:
		l.close()
	js=None

	if options.buildsetonly: return

	for itr in range(0,options.niters):

		from EMAN2PAR import EMTaskCustomer

		for eo in ["even", "odd"]:
			
			if options.nogs:
				threedname="{}/threed_{:02d}.hdf".format(path, itr)
			else:
				threedname="{}/threed_{:02d}_{}.hdf".format(path, itr, eo)
			
			lstname="{}/ali_ptcls_{:02d}_{}.lst".format(path, itr, eo)
			lst=LSXFile(lstname, True)
			m=EMData(threedname)
			
			m.process_inplace('normalize.edgemean')
			#m.process_inplace("threshold.belowtozero",{"minval":0.5})
			
			pinfo=[]
			if options.debug: nptcl=options.threads*8
			else: nptcl=lst.n
			for i in range(nptcl):
				pinfo.append(lst.read(i))
			lst=None
			
			etc=EMTaskCustomer(options.parallel)
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

			lname="{}/ali_ptcls_{:02d}_{}.lst".format(path, itr+1, eo)
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
			cmd="e2make3dpar.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads {trd} --outsize {bx} --apix {apx} --mode gauss_var --keep {kp} --sym {sm}".format(
				inp=lname, 
				out="{}/threed_{:02d}_{}.hdf".format(path, itr+1, eo),
				bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=options.sym, trd=num_cpus)
			
			run(cmd)

		s = ""
		jd = js_open_dict("{}/0_subtlt_params.json".format(path))
		if jd.has_key("goldstandard"): 
			if jd["goldstandard"] > 0: 
				s += " --align"
		if jd.has_key("setsf"):
			s += " --setsf {}".format(jd['setsf']) #options.setsf)
		if jd.has_key("localfilter"):
			s += " --tophat local"
		if jd.has_key("mass"): 
			s += " --mass {}".format(jd["mass"])
		else: 
			s += " --mass 1000"
		if jd.has_key("sym"):
			if jd["sym"] != "c1": 
				s += " --sym {}".format(jd["sym"])
		msk = jd["mask"] #{}/mask_tight.hdf".format(path)
		if os.path.isfile(msk):
			s += " --automask3d mask.fromfile:filename={}".format(msk)
		else:
			s += " --automask3d {}".format(msk)
		jd.close()

		# get target resolution from last iteration map
		ref=os.path.join(path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(path, "fsc_masked_{:02d}.txt".format(itr)))
		rs=1./fsc[fsc[:,1]<0.3, 0][0]
		curres=rs*.5

		#os.system("rm {}/mask*.hdf {}/*unmasked.hdf".format(path, path))
		ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --restarget {} --threads {} {}".format(
			os.path.join(path, "threed_{:02d}_even.hdf".format(itr+1)), 
			os.path.join(path, "threed_{:02d}_odd.hdf".format(itr+1)), 
			os.path.join(path, "threed_{:02d}.hdf".format(itr+1)), itr+1, curres, options.threads, s)
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

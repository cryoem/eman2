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

def do_copy(a,b):
	try:
		copy2(a,b)
	except:
		print("!!\tCannot find {}".format(a))
	return

def main():
	
	usage="""prog --path <path to previous spt or subtlt refinement> [options]
	This program will run subtilt refinement based on previous subtomogram or subtilt refinement."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_header(name="orblock0", help='Just a visual separation', title="Inputs", row=0, col=1, rowspan=1, colspan=3, mode="model")


	parser.add_argument("--path", type=str,help="Path to the previous spt/subtlt refinement", default=None, guitype='filebox',  browser="EMBrowserWidget(withmodal=True,multiselect=False)",row=1, col=0,rowspan=1, colspan=2)
	parser.add_argument("--iter", type=int,help="Start from iteration X of previous refinement", default=-1, guitype='intbox',row=1, col=2,rowspan=1, colspan=1)
	
	
	#####################
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=3, mode="model")

	
	
	parser.add_argument("--niters", type=int,help="Run this many iterations. Default is 4.", default=4, guitype='intbox',row=3, col=0,rowspan=1, colspan=1)

	parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.8", default=0.8, guitype='floatbox',row=3, col=1,rowspan=1, colspan=1)

	parser.add_argument("--maxalt", type=float,help="max altitude to insert to volume", default=45.0, guitype='floatbox',row=3, col=2,rowspan=1, colspan=1)	
	
	parser.add_argument("--mask", type=str, default="Auto" ,help="Refinement and reprojection masking.",guitype='strbox',row=4, col=0,rowspan=1, colspan=2)	
	
	parser.add_argument("--nogs", action="store_true", default=False ,help="Skip gold standard. This is not a great idea...", guitype='boolbox',row=4, col=2,rowspan=1, colspan=1)
	

	parser.add_argument("--threads", type=int,help="Number of CPU threads to use. Default is 12.", default=12, guitype='intbox',row=5, col=2,rowspan=1, colspan=1)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:12", guitype='strbox',row=5, col=0,rowspan=1, colspan=2)

	parser.add_argument("--buildsetonly", action="store_true", default=False ,help="will only prepare particle set for the refinement but skip the actual refinement process.",guitype='boolbox',row=6, col=0,rowspan=1, colspan=1)
	parser.add_argument("--resume", action="store_true", default=False ,help="continue from previous run",guitype='boolbox',row=6, col=1,rowspan=1, colspan=1)
	
	parser.add_argument("--reproject", action="store_true", default=False ,help="Reproject 3D particles into 2D particles.")
	
	parser.add_argument("--reprj_offset", type=str, default="" ,help="Offset translation before reprojection")
	parser.add_argument("--reprj_clip", type=int, default=-1 ,help="clip after reprojection")

	parser.add_argument("--tophat", type=str, default="auto" ,help="Filter option for refine_postprocess. auto: same as spt refinement; local; global;")
	parser.add_argument("--refineastep", type=float,help="Mean angular variation for refine alignment", default=1.)
	parser.add_argument("--refinentry", type=int,help="number of starting points for refine alignment", default=8)
	parser.add_argument("--maxshift", type=int,help="maximum shift allowed", default=8)


	parser.add_argument("--padby", type=float,default=1.5, help="pad by factor. default is 1.5")
	parser.add_argument("--output", type=str,help="Write results to this directory. We do not recommend changing this.", default="subtlt")#, guitype='intbox',row=2, col=1,rowspan=1, colspan=1)
	parser.add_argument("--debug", action="store_true", default=False ,help="Turn on debug mode. This will only process a small subset of the data (threads * 8 particles)")
	parser.add_argument("--localnorm",action="store_true",help="local normalization. do not use yet....",default=False)
	parser.add_argument("--sym", type=str,help="symmetry. will use symmetry from spt refinement by default", default="")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	parser.add_argument("--transonly", action="store_true", default=False ,help="only refine translation",guitype='boolbox',row=7, col=0,rowspan=1, colspan=1)


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

	
	if not os.path.isfile(os.path.join(oldpath,"threed_{:02d}.hdf".format(itr))):
		print("Could not locate {}/threed_{:02d}.hdf".format(oldpath,itr))
		print("Please specify the iteration number (--iter) of a completed subtomogram refinement.")
		return
	
	
	if "0_subtlt_params.json" in os.listdir(oldpath):
		print("Continuing from a subtilt refinement...")
		fromspt=False
	else:
		print("Start from a spt refinement...")
		fromspt=True
		
	if options.resume:
		if fromspt:
			print("Cannot resume from a spt refinement...")
			return
		path=oldpath
		e=EMData(os.path.join(path,"threed_{:02d}.hdf".format(itr)))
	else:
		path = make_path(options.output)
		print("Writing in {}...".format(path))
	
		oldmap = os.path.join(oldpath,"threed_{:02d}.hdf".format(itr))
		oem = os.path.join(oldpath,"threed_{:02d}_even.hdf".format(itr))
		oom = os.path.join(oldpath,"threed_{:02d}_odd.hdf".format(itr))
		oldfsc = os.path.join(oldpath, "fsc_masked_{:02d}.txt".format(itr))

		do_copy(oldmap,os.path.join(path,"threed_00.hdf"))
		do_copy(oldfsc, os.path.join(path, "fsc_masked_00.txt"))
		do_copy(oem,os.path.join(path,"threed_00_even.hdf"))
		do_copy(oom,os.path.join(path,"threed_00_odd.hdf"))
		
		if fromspt:
			oldparm = os.path.join(oldpath,"particle_parms_{:02d}.json".format(itr))
			do_copy(oldparm,os.path.join(path,"particle_parms_00.json"))
		else:
			for eo in ["even", "odd"]:
				oali = os.path.join(oldpath,"ali_ptcls_{:02d}_{}.lst".format(itr, eo))
				do_copy(oali,os.path.join(path,"ali_ptcls_00_{}.lst".format(eo)))


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
		if "radref" in oldjd:
			jd["radref"]=oldjd["radref"]
			
		if fromspt:
			options.ptclkeep=oldjd["pkeep"]
			
		oldjd.close()
	else:
		print("Cannot find {}. Using default parameters.".format(sptparms))
		jd["mass"] = -1
		jd["sym"] = "c1"
		jd["localfilter"]=False
		jd["mask"]=""
			
		if fromspt:
			options.ptclkeep=0.9
		
	
	if options.mask.lower()!="none":
		#print("Overwritting masking")
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

	if fromspt and not options.reproject:
		js=js_open_dict(os.path.join(path,"particle_parms_00.json"))
		k=list(js.keys())[0]
		src=eval(k)[0]
		
		print("loading {} 3D particles from {}".format(len(js.keys()), base_name(src)))
		print("   box size {}, apix {:.2f}".format(bxsz, apix))
		
		e=EMData(src, 0, True)
		pjname=e["class_ptcl_src"]
		pj=EMData(pjname, 0, True)
		transmult=1.0
		if e["apix_x"]!=pj["apix_x"]:
			print("Apix mismatch between 2D and 3D particles: {:.2f} vs {:.2f}".format(pj["apix_x"], e["apix_x"]))
			transmult=e["apix_x"]/pj["apix_x"]
			print("  Will scale translation by {:.1f} ...".format(transmult))

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
		
			thr=int(len(score)*options.ptclkeep)
			simthr=np.sort(score)[thr]
			print("  Removing {} bad 3D particles...".format(thr))
			
		else:
			simthr=10000
			
		print("Building aligned 2D particle set from spt alignment...")
		for ky in js.keys():
			
			src,ii=eval(ky)
			e=EMData(src, ii, True)
			fname=e["class_ptcl_src"]
			ids=e["class_ptcl_idxs"]
			#ky="('{}', {})".format(src, ii)
			dic=js[ky]
			xali=dic["xform.align3d"]
			xali.set_trans(xali.get_trans()*transmult)
			scr=float(dic["score"])
			if scr>simthr:
				continue
			
			if "__even" in src:
				eo=0
			elif "__odd" in src:
				eo=1
			else:
				eo=ii%2
			
			
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

	if fromspt and options.reproject:
		print("Reprojecting 3D particles...")
		
		if options.reprj_offset!="":
			options.reprj_offset=[int(i) for i in options.reprj_offset.split(',')]
			print("Offset by {} before reprojection".format(options.reprj_offset))
		else:
			options.reprj_offset=None
		
		if options.mask.lower()=="auto":
			mask=EMData(os.path.join(oldpath,"mask.hdf"))
			mask.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.1})
			mask.add(0.5)
			mask.div(1.5)
			mask.write_image(os.path.join(path,"mask_reproj.hdf"))
			
		elif options.mask.lower()=="none":
			mask=None
			
		else:
			mask=EMData(options.mask)
			
		
		ptclfile=[os.path.join(path,"particles_reproj_{}.hdf".format(eo)) for eo in ["even","odd"]]
		for p in ptclfile:
			try: os.remove(p)
			except: pass
		
		lstfile=[os.path.join(path, "ali_ptcls_00_{}.lst".format(eo)) for eo in ["even", "odd"]]
		for p in lstfile:
			try: os.remove(p)
			except: pass


		js=js_open_dict(os.path.join(path,"particle_parms_00.json"))
		jsdata=dict(js.data)
		js.close()
			
		if options.ptclkeep<1.0:
			score=[]
			for k in list(jsdata.keys()):
				score.append(float(jsdata[k]["score"]))
		
			simthr=np.sort(score)[int(len(score)*options.ptclkeep)]
			print("removing bad particles with score > {:.2f}...".format(simthr))
			
		else:
			simthr=10000
			
		keys=list(jsdata.keys())
		nt=options.threads
		if options.debug:
			keys=keys[:nt*32]
		thrds=[threading.Thread(target=do_reproject,
				args=(jsdata, keys[i::nt], ptclfile, options, mask, simthr)) for i in range(nt)]

		print(len(thrds)," threads")
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==options.threads ) : time.sleep(.1)
				#if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(1)
		
		for t in thrds:
			t.join()

		print("Now writing list files...")
		#### write a list file
		for k in [0,1]:
		
			n=EMUtil.get_image_count(ptclfile[k])
			lst=LSXFile(lstfile[k], False)
			for i in range(n):
				e=EMData(ptclfile[k], i, True)
				rot=e["xform.align3d"]
				lst.write(-1, i, ptclfile[k], str(rot.get_params("eman")))
			
			lst=None
			
			
			
	if options.buildsetonly: return

	
	if options.resume:
		starti=itr
	else:
		starti=0
	
	for itr in range(starti,options.niters+starti):

		

		for eo in ["even", "odd"]:
			
			if options.nogs:
				threedname=os.path.join(path, "threed_{:02d}.hdf".format(itr))
			else:
				threedname=os.path.join(path, "threed_{:02d}_{}.hdf".format(itr, eo))
			
			lstname=os.path.join(path, "ali_ptcls_00_{}.lst".format(eo))
			lname=os.path.join(path, "ali_ptcls_{:02d}_{}.lst".format(itr+1, eo))
			threedout=os.path.join(path, "threed_{:02d}_{}.hdf".format(itr+1, eo))
			
			cmd="e2spt_tiltrefine_oneiter.py --ptclin {} --ptclout {} --ref {} --threedout {} --keep {} --threads {} --parallel {} --refineastep {} --refinentry {} --maxshift {} --padby {} --sym {}".format(lstname, lname, threedname, threedout,  options.keep, options.threads, options.parallel, options.refineastep, options.refinentry, options.maxshift, options.padby, jd["sym"])
			if options.debug: 
				cmd+=" --debug"
				
			if options.transonly: 
				cmd+=" --transonly"
				
			run(cmd)
			
			run("e2proc3d.py {} {}".format(threedout, os.path.join(path, "threed_raw_{}.hdf".format(eo))))

		s = ""
		
		if "goldstandard" in jd: 
			if jd["goldstandard"] > 0: 
				s += " --align"
		if "setsf" in jd:
			s += " --setsf {}".format(jd['setsf']) #options.setsf)
		
		
		if options.tophat=="auto" and ("localfilter" in jd) and jd["localfilter"]==True:
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
		
		try:
			rs=1./fsc[fsc[:,1]<0.3, 0][0]
		except:
			rs=10
			
		curres=rs*.5
		
		even=os.path.join(path, "threed_{:02d}_even.hdf".format(itr+1))
		odd=os.path.join(path, "threed_{:02d}_odd.hdf".format(itr+1))


		ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --restarget {} --threads {} --sym {} --mass {} {}".format(even, odd, os.path.join(path, "threed_{:02d}.hdf".format(itr+1)), itr+1, curres, options.threads, jd["sym"], jd["mass"], s)
		run(ppcmd)
		
		if options.localnorm:
			for f in [even, odd]:
				run("e2proc3d.py {} {} --process normalize --process normalize.local:threshold=1:radius=16".format(f,f))
				
			run("e2proc3d.py {} {} --addfile {} --process normalize".format(
				even, os.path.join(path, "threed_{:02d}.hdf".format(itr+1)), odd))

		fsc=np.loadtxt(os.path.join(path, "fsc_masked_{:02d}.txt".format(itr+1)))

		print("Resolution (FSC<0.3) is ~{:.1f} A".format(rs))
				
	E2end(logid)

def do_reproject(js, keys, outfiles, options, mask=None, simthr=10000.):
	
	for ik, ky in enumerate(keys):

		dic=js[ky]
		scr=float(dic["score"])
		if scr>simthr:
			continue
		
		src,ii=eval(ky)
		e=EMData(src, ii)
		fname=e["class_ptcl_src"]
		ids=e["class_ptcl_idxs"]
		#ky="('{}', {})".format(src, ii)
		
		xali=Transform(dic["xform.align3d"])
		if options.reprj_offset:
			of=options.reprj_offset
			xali.translate(-of[0], -of[1], -of[2])
			
		if mask:
			m=mask.copy()
			m.transform(xali.inverse())
			e.mult(m)
		
		
		if "__even" in src:
			eo=0
		elif "__odd" in src:
			eo=1
		else:
			eo=ii%2
			
		outname=outfiles[eo]
		
		for i in ids:
			try:
				m=EMData(fname, i, True)
			except:
				continue
			
			
			xf=m["xform.projection"]
			dc=xf.get_params("eman")
			if abs(dc["alt"])>options.maxalt:
				continue
			
			rot=xf*xali.inverse()
			pj=e.project("standard", xf)
			pj.set_attr_dict(m.get_attr_dict())
			
			ts=rot.get_trans()
			tr=np.round(ts)
			pj.translate(-tr[0], -tr[1],0)
			rot.set_trans((ts-tr).tolist())
			
			pj["xform.projection"]=xf
			pj["xform.align3d"]=rot
			if options.reprj_clip>0:
				clip=options.reprj_clip
				pj=pj.get_clip(Region((pj["nx"]-clip)//2, (pj["ny"]-clip)//2, clip, clip))
			pj.write_image(outname, -1)



def run(cmd):
	print(cmd)
	launch_childprocess(cmd)



	
if __name__ == '__main__':
	main()

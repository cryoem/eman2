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

def split_xf(xf):
	dc=xf.get_params('eman')
	xf3d=Transform({"type":"eman","alt":dc["alt"], "az":dc["az"]})
	
	xf2d=Transform({"type":"2d","alpha":dc["phi"], "tx":dc["tx"], "ty":dc["ty"]})
	return xf3d, xf2d

def merge_xf(xf3d, xf2d):
	dc2=xf2d.get_params('eman')
	dc3=xf3d.get_params('eman')
	xf=Transform({"type":"eman","alt":dc3["alt"], "az":dc3["az"], "phi":dc2["phi"], "tx":dc2["tx"], "ty":dc2["ty"]})
	
	return xf

def refine_ali(ids, pinfo, m, jsd, options):
	sz=m["nx"]
	for ii in ids:
		l=pinfo[ii]
		e=EMData(l[1], l[0])
		b=e["nx"]
		e=e.get_clip(Region(old_div((b-sz),2), old_div((b-sz),2), sz,sz)).process("normalize")
		dc=eval(l[2])
		try:
			dc.pop('score')
		except:
			pass
		
		xf=Transform(dc)
		pj=m.project("standard", xf).process("normalize.edgemean")
		
		#### do snr weighting if ctf info present
		if e.has_attr("ctf"):
			ctf=e["ctf"]
			ctf.bfactor=500 #### use a fixed b factor for now...
			ctf.dsbg=old_div(1.,(e["apix_x"]*e["nx"]))
			s=np.array(ctf.compute_1d(e["nx"], ctf.dsbg, Ctf.CtfType.CTF_INTEN ))
			s[:np.argmax(s)]=np.max(s)
			s=np.maximum(s, 0.001)
			ctf.snr=s.tolist()
			e["ctf"]=ctf
			eali=e.align("refine", pj, {"maxshift":8, "maxiter":50}, "frc", {"snrweight":1, "maxres":options.maxres, "minres":500})
		else:
			#print("missing ctf...")
			eali=e.align("refine", pj, {"maxshift":8, "maxiter":50}, "frc", {"minres":500, "maxres":options.maxres})
		al=eali["xform.align2d"]
		
		x0,x1=split_xf(xf)
		x1=al.inverse()*x1
		xf1=merge_xf(x0, x1)
	


		scr=eali.cmp("frc", pj, {"maxres":options.maxres})
		dc=xf1.get_params("eman")
		dc["score"]=float(scr)
		jsd.put((ii, dc))
		
def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSPTParticleTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=1, col=1, rowspan=1, colspan=1)

	parser.add_argument("--maxres", type=float,help="max resolution for comparison", default=20.0, guitype='floatbox',row=2, col=0,rowspan=1, colspan=1)
	parser.add_argument("--iter", type=int,help="iteration", default=1, guitype='intbox',row=2, col=1,rowspan=1, colspan=1)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=2, col=2,rowspan=1, colspan=1)

	parser.add_argument("--padby", type=float,help="pad by factor. default is 2", default=2., guitype='floatbox',row=3, col=0,rowspan=1, colspan=1)
	parser.add_argument("--keep", type=float,help="propotion of tilts to keep. default is 0.5", default=0.5, guitype='floatbox',row=3, col=1,rowspan=1, colspan=1)
	parser.add_argument("--maxalt", type=float,help="max altitude to insert to volume", default=90.0, guitype='floatbox',row=3, col=2,rowspan=1, colspan=1)

	parser.add_argument("--dopostp", action="store_true", default=False ,help="Do post processing routine", guitype='boolbox',row=4, col=0,rowspan=1, colspan=1)
	parser.add_argument("--nogs", action="store_true", default=False ,help="skip gold standard...", guitype='boolbox',row=4, col=1,rowspan=1, colspan=1)
	parser.add_argument("--fromlast", action="store_true", default=False ,help="continue from a previous tilt refine", guitype='boolbox',row=5, col=0,rowspan=1, colspan=1)

	#parser.add_argument("--unmask", action="store_true", default=False ,help="use unmasked map as references", guitype='boolbox',row=4, col=1,rowspan=1, colspan=1)
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=4, col=2,rowspan=1, colspan=1)

	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	path=options.path
	itr=options.iter
	if not path: 
		print("No input path. Exit.")
		return

	#### start from a spt_align - spt_average
	if options.fromlast:
		
		e=EMData("{}/threed_{:02d}.hdf".format(path, itr-1))
		bxsz=e["nx"]
		apix=e["apix_x"]
		print("box size {}, apix {:.2f}".format(bxsz, apix))
		
	else:
		
		js=js_open_dict(path+"/particle_parms_{:02d}.json".format(itr))
		k=list(js.keys())[0]
		src=eval(k)[0]
		e=EMData(src, 0,True)
		bxsz=e["nx"]
		apix=e["apix_x"]
		print("loading 3D particles from {}".format(src))
		print("box size {}, apix {:.2f}".format(bxsz, apix))
		
		fscs=[os.path.join(path,f) for f in os.listdir(path) if f.startswith("fsc") and f.endswith("{:02d}.txt".format(itr))]
		for f in fscs:
			os.rename(f, f[:-4]+"_raw.txt")
			
		for eo in ["", "_even", "_odd"]:
			os.rename("{}/threed_{:02d}{}.hdf".format(path, itr, eo), 
					"{}/threed_{:02d}_raw{}.hdf".format(path, itr, eo))
		
		lname=[os.path.join(path, "ali_ptcls_{:02d}_{}.lst".format(itr,eo)) for eo in ["even", "odd"]]
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
	
	for eo in ["even", "odd"]:
		
		if options.fromlast:
			#### start from a previous tilt refine
			lstname="{}/ali_ptcls_refine_{:02d}_{}.lst".format(path, itr-1, eo)
			if options.nogs:
				threedname="{}/threed_{:02d}.hdf".format(path, itr-1)
			else:
				threedname="{}/threed_{:02d}_{}.hdf".format(path, itr-1, eo)
			
			#print(lstname, threedname)
			lst=LSXFile(lstname, True)
			m=EMData(threedname)
		
		else:
		
			lst=LSXFile("{}/ali_ptcls_{:02d}_{}.lst".format(path, itr, eo), True)
			
			
			#if options.unmask:
				#m=EMData("{}/threed_{}_unmasked.hdf".format(path, eo))
				#m.process_inplace("mask.soft",{"outer_radius":-4})
			#else:
			if options.nogs:
				m=EMData("{}/threed_{:02d}_raw.hdf".format(path, itr))
			else:
				m=EMData("{}/threed_{:02d}_raw_{}.hdf".format(path, itr, eo))
			
		m.process_inplace('normalize.edgemean')
		#m.process_inplace("threshold.belowtozero",{"minval":0.5})
		
		pinfo=[]
		nptcl=lst.n
		for i in range(nptcl):
			pinfo.append(lst.read(i))
		lst=None

		jsd=queue.Queue(0)
		jobs=[]
		print("Refining {} set with {} 2D particles..".format(eo, nptcl))
		batchsz=100
		for tid in range(0,nptcl,batchsz):
			ids=list(range(tid, min(tid+batchsz, nptcl)))
			jobs.append([ids, pinfo, m, jsd, options])

		thrds=[threading.Thread(target=refine_ali,args=(i)) for i in jobs]
		thrtolaunch=0
		tsleep=threading.active_count()

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
				ndone+=1
				if ndone%2000==0:
					print("\t{}/{} finished.".format(ndone, nptcl))

		for t in thrds: t.join()
		#np.savetxt("tmpout1.txt", dics, fmt='%s')
			

		allscr=np.array([d["score"] for d in dics])
		print(np.min(allscr), np.mean(allscr), np.max(allscr), np.std(allscr))
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

		lname="{}/ali_ptcls_refine_{:02d}_{}.lst".format(path, itr, eo)
		try: os.remove(lname)
		except: pass
		lout=LSXFile(lname, False)
		for i, d in enumerate(dics):
			d["score"]=float(allscr[i])
			l=pinfo[i]
			lout.write(-1, l[0], l[1], str(d))

		lout=None

		pb=options.padby
		cmd="make3dpar_rawptcls.py --input {inp} --output {out} --pad {pd} --padvol {pdv} --threads 12 --outsize {bx} --apix {apx} --mode gauss_2 --keep {kp} --sym {sm}".format(
			inp=lname, 
			out="{}/threed_{:02d}_{}.hdf".format(path, itr, eo),
			bx=bxsz, pd=int(bxsz*pb), pdv=int(bxsz*pb), apx=apix, kp=options.keep, sm=options.sym)
		
		run(cmd)
		
	if options.dopostp:
		sfx="{}/threed_{:02d}".format(path, itr )
		if len(options.setsf)>0:
			sf=" --setsf {}".format(options.setsf)
		else:
			sf=""
		cmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {} --mass 1000.0 --restarget 10.0 --sym {}  --align {}".format(
			sfx+"_even.hdf", sfx+"_odd.hdf", sfx+".hdf", itr, options.sym, sf)
		
		run(cmd)
	
		fscs=[os.path.join(path,f) for f in os.listdir(path) if f.startswith("fsc") and f.endswith("{:02d}.txt".format(itr))]
		for f in fscs:
			os.rename(f, f[:-4]+"_ali.txt")
			
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()


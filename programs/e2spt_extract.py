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
from EMAN2_utils import *



def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="Specify tomograms from which you wish to extract boxed particles.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="extract")
	parser.add_argument("--boxsz_unbin", type=int,help="box size in unbinned tomogram", default=-1, guitype='intbox',row=2, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--label", type=str,help="Only extract particle with this name. Leave blank to extract all particles.", default=None, guitype='strbox',row=2, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--newlabel", type=str,help="Label of output particles. Same as original particle label by default.", default="", guitype='strbox',row=6, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=3, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=4, col=1,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--maxtilt", type=int,help="max tilt", default=100, guitype='intbox',row=4, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--padtwod", type=float,help="padding factor", default=2.0, guitype='floatbox',row=5, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--noctf", action="store_true", default=False ,help="skip ctf correction..", guitype='boolbox',row=5, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--wiener", action="store_true", default=False ,help="wiener filter the particles using ctf information..", guitype='boolbox',row=6, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--alltomograms", action="store_true", default=False ,help="use all tomograms.", guitype='boolbox',row=1, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--dotest", action="store_true", default=False ,help="only make 1 batch of subtomograms for testing")
	parser.add_argument("--curves", action="store_true", default=False ,help="extract particles from saved curves")
	parser.add_argument("--curves_overlap", type=float, help="fraction of overlap when generating particle along curves. default is 0.5",default=0.5)

	parser.add_argument("--shrink", type=int, help="Shrinking factor for output particles. Default is 1 (no shrink)",default=1, guitype='intbox',row=8, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--tltkeep", type=float,help="keep a fraction of tilt images with good score determined from tomogram reconstruction", default=1.0, guitype='floatbox',row=8, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--rmbeadthr", type=float,help="remove 2d particles with high contrast object beyond N sigma at 100A. Note that this may result in generating fewer particles than selected. Default is -1 (include all particles). 5 might be a good choice for removing gold beads but may need some testing...", default=-1, guitype='floatbox',row=9, col=0, rowspan=1, colspan=1, mode="extract")
	
	parser.add_header(name="orblock1", help='Just a visual separation', title="Re-extraction from spt", row=10, col=0, rowspan=1, colspan=1, mode="extract")

	parser.add_argument("--jsonali", type=str,help="re-extract particles using a particle_param_xx json file from a spt alignment", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=11, col=0,rowspan=1, colspan=2, mode="extract")
	#parser.add_argument("--alioffset", type=str,help="coordinate offset when re-extract particles. (x,y,z)", default="0,0,0", guitype='strbox', row=12, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--postxf", type=str,help="a file listing post transforms", default="")


	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	
	if len(options.jsonali)>0:
		print("re-extracting particles based on previous alignment. Ignoring particle/tomogram input...")
		#try:
			#offset=list(float(i) for i in options.alioffset.split(','))
		#except:
			#offset=(0,0,0)
		
		#print("extracting {} sub-particles per particle, with offset".format(len(offset)//3))
		#print(np.array(offset).reshape((-1,3)))
		allxfs=parse_json(options.jsonali, options.postxf)
		for fname in allxfs.keys():
			#### it seems options is changed inplace somewhere...
			(options, args) = parser.parse_args()
			do_extraction(fname, options, allxfs[fname])
		
		return
	
	

	if options.alltomograms:
		fld="tomograms/"
		args=[fld+f for f in os.listdir(fld) if f.endswith(".hdf")]
		### do not count a tomogram twice when multiple versions exist
		uq, uid=np.unique([base_name(a) for a in args], return_index=True)
		args=[args[i] for i in uid]
		#print(args)
	
		
	if len(args)==0:
		print("No input. Exit.")
		return

	
	if len(args)==1:
		
		do_extraction(args[0], options)
	else:
		print("Processing {} files in sequence..".format(len(args)))
		
		#cmd=sys.argv
		#opt=' '.join([s for s in cmd[1:] if s not in args])
		#opt=opt.replace("--alltomograms","")
		for a in args:
			do_extraction(a, options)
			
	print("Finished.")
	
	
	E2end(logid)
	return
	
	
	
def do_extraction(pfile, options, xfs=[]):
	#pfile=args[0]
	print("Reading from {}...".format(pfile))
	options.boxsz=options.boxsz_unbin // int(options.shrink)
	
	#### reading alignment info...
	js=js_open_dict(info_name(pfile))
	ttparams=np.array(js["tlt_params"])
	ttparams[:,:2]/=options.shrink
	
	if options.noctf==False and "defocus" in js:
		#### read ctf info when exist
		defocus=np.array(js["defocus"])
		phase=np.array(js["phase"])
		voltage=float(js["voltage"])
		cs=float(js["cs"])
		print("CTF information exists. Will do phase flipping..")
	else:
		defocus=[]
		phase=[]
		
	tltkeep=[]
	if options.tltkeep<1.0:
		if js.has_key("ali_loss"):
			aliloss=js["ali_loss"]
			tltkeep=np.argsort(aliloss)[:int(len(ttparams)*options.tltkeep)]
		else:
			print("warning: --tltkeep specified, but cannot find ali_loss in info. setting tltkeep to 1.")
	
			
	tfile=str(js["tlt_file"])
	options.tltfile=tfile
	ptcl=EMData(pfile, 0, True)
	img=EMData(tfile,0, True)
	apix_ptcl=ptcl["apix_x"]
	apix_tlt=img["apix_x"] * float(options.shrink)
	
	
	try:
		zshift=ptcl["zshift"]#/2
	except:
		zshift=0
	#yt=ptcl["ytilt"]
	yt=0
	#options.ytilt=ttparams[np.argmin(abs(ttparams[:,3]-yt)),3]
	options.ytilt=0
	scale=apix_ptcl/apix_tlt

	print("Scaling factor: {:.1f}, y-tilt: {:.1f}, z-shift: {:d}.".format(scale, options.ytilt, int(zshift)))
	
	
	#### reading particle location from json file corresponding to the tomogram input
	ptclpos=[]
	nptcl=EMUtil.get_image_count(pfile)
	e=EMData(pfile, 0, True)
	if len(xfs)>0:
		if options.newlabel=="":
			lab=pfile[pfile.rfind('__')+2:pfile.rfind('.')]
			lab+="_reextract"
		else:
			lab=options.newlabel
			
		outname=str(base_name(pfile)+"__"+lab+".hdf")
		sz=int(options.boxsz//2)
		if scale!=1:
			for t in xfs:
				t.set_trans(t.get_trans()*scale)
		towrite=[[xfs, outname, sz]]
		
		
	elif e["nz"]!=e["nx"] and nptcl==1:
		print("Reading particle location from a tomogram...")
		js=js_open_dict(info_name(pfile))
		towrite=[]
		if options.curves:
			print("Generating particles along curves...")
			overlap=options.curves_overlap
			if overlap>=1 or overlap<0:
				print("Overlap has to be in [0,1)")
				return
			
			if js.has_key("curves") and len(js["curves"])>0:
				pts=np.array(js["curves"]).copy()
				js.close()
				
				if "apix_unbin" in js:
					pts[:,:3]/=options.shrink
				else:
					pts[:,:3]-=[e["nx"]//2, e["ny"]//2, e["nz"]//2]
					pts[:,:3]*=scale
					pts[:,2]-=zshift
				
				if options.newlabel=="":
					lab="curve"
				else:
					lab=options.newlabel
					
				if options.shrink>1:
					lab+="_bin{:d}".format(options.shrink)
				outname=str(base_name(pfile)+"__"+lab+".hdf")
				
				sz=int(options.boxsz//2)
				
				
				bxs=[]
				drs=[]
				for li in np.unique(pts[:,3]):
					pt=pts[pts[:,3]==li][:,:3].copy()
					pt=pt[np.append(True, np.linalg.norm(np.diff(pt, axis=0), axis=1)>0.1)]
					ln=np.linalg.norm(pt[-1]-pt[0])
					#     print np.round(ln)//2
					if len(pt)<2: continue
					ipt=interp_points(pt, npt=int(np.round(ln/options.boxsz/(1-overlap))))
					
					if len(ipt)<4: continue
					bxs.append(ipt[1:-1])
					drs.append(ipt[2:]-ipt[:-2])
				
				bxs=np.vstack(bxs)
				drs=np.vstack(drs)
				bxs=np.hstack([bxs,drs])
				
				towrite.append((bxs, outname, sz))
		
		elif "class_list" in js and "boxes_3d" in js:
			clslst=js["class_list"]
			boxes=js["boxes_3d"]
			for ky in list(clslst.keys()):
				val=clslst[ky]
				if options.label:
					if str(val["name"])!=options.label:
						continue
						
				bxs=np.array([[b[0], b[1], b[2]] for b in boxes if b[5]==int(ky)], dtype=float)
				if len(bxs)==0: continue
				
				if "apix_unbin" in js:
					bxs/=options.shrink
					if options.boxsz<0:
						sz=int(val["boxsize"])//2
					else:
						sz=int(options.boxsz//2)
				else:
					bxs-=[e["nx"]//2, e["ny"]//2, e["nz"]//2]
					bxs*=scale
					bxs[:,2]-=zshift
					if options.boxsz<0:
						sz=int(val["boxsize"])*scale//2
					else:
						sz=int(options.boxsz//2)
						
					
				if options.newlabel=="":
					lab=val["name"]
				else:
					lab=options.newlabel
					
				if options.shrink>1:
					lab+="_bin{:d}".format(options.shrink)
				outname=str(base_name(pfile)+"__"+lab+".hdf")
				
				towrite.append((bxs, outname, sz))
				print("{} : {} boxes, unbinned box size {}".format(val["name"], len(bxs), int(sz*2))) 
		
		if len(towrite)==0:
			print("No particles. exit..")
			return
	
	else:
		print("Reading particle location from image stack..")
		for i in range(nptcl):
			ptcl=EMData(pfile, i, True)
			ptclpos.append(ptcl["ptcl_source_coord"])
		
	
		sfx=pfile[pfile.find("__"):]
		for i in range(2,10): sfx=sfx.replace("_bin{:d}".format(i),"")
		if options.shrink>1:
			sfx+="_bin{:d}".format(options.shrink)
		outname=base_name(pfile)+sfx
		
		
		ptclpos=np.array(ptclpos, dtype=float)
		
		ptclpos*=scale
		ptclpos[:,2]-=zshift#-2
		
		print("Reading {} particles".format(len(ptclpos)))
		
		if options.boxsz<0:
			boxsz=int(ptcl["nx"]//2*scale)
		else:
			boxsz=int(options.boxsz//2*scale)
		
		
		
		towrite=[(ptclpos,outname, boxsz)]
	
	
	print("Reading tilt series file: {}".format(tfile))
	#### make sure this works on image stack or mrc volume
	img=EMData(tfile,0)
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	for m in imgs: 
		if options.shrink>1:
			m.process_inplace("math.fft.resample",{"n":options.shrink})
		m.process_inplace("normalize")
	ntlt=len(imgs)
		
	try: os.mkdir("particles3d")
	except: pass
	try: os.mkdir("particles")
	except: pass

	for pinfo in towrite:
		ptclpos, outname, boxsz=pinfo
		if options.dotest:
			nptcl=options.threads
		else:
			nptcl=len(ptclpos)
	
		options.output=os.path.join("particles3d", outname)
		options.output2d=os.path.join("particles", outname)
		options.pad=pad=good_size(boxsz*2*options.padtwod)
		
		print("Writing {} particles to {}".format(nptcl, outname))
		print("Box size {}, padded to {}".format(int(boxsz*2), pad))
		
		try: os.remove(options.output)
		except: pass
		try: os.remove(options.output2d)
		except: pass

		jsd=queue.Queue(0)
		jobs=[]
		
		batchsz=4
		if len(defocus)>0:
			ctf=[defocus, phase, voltage, cs]
		else:
			ctf=[]
			
		for tid in range(0,nptcl,batchsz):
			ids=list(range(tid, min(tid+batchsz, nptcl)))
			jobs.append([jsd, ids, imgs, ttparams, pinfo, options, ctf, tltkeep])
		
		
		thrds=[threading.Thread(target=make3d,args=(i)) for i in jobs]

		thrtolaunch=0
		tsleep=threading.active_count()
		ndone=0
		while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
			if thrtolaunch<len(thrds):
				while (threading.active_count()==options.threads+tsleep ) : 
					#print threading.active_count(), options.threads, tsleep, thrtolaunch, len(thrds)
					time.sleep(.1)
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(.2)
			
			
			while not jsd.empty():
				pid, threed, projs=jsd.get()
				#pjids=[pid*ntlt+i for i in range(ntlt)]
				try: pji=EMUtil.get_image_count(options.output2d)
				except: pji=0
				pjids=[]
				for i,pj in enumerate(projs):
					pj.write_image(options.output2d, pji)
					pjids.append(pji)
					pji+=1
				threed["class_ptcl_src"]=options.output2d
				threed["class_ptcl_idxs"]=pjids
				#if options.rmbeadthr>0:
				#### give up maintaining the order since there are empty particles...
				#pid=-1
				threed.write_image(options.output, -1)
				ndone+=1
				#if ndone%10==0:
				sys.stdout.write("\r{}/{} finished.".format(ndone, nptcl))
				sys.stdout.flush()

		for t in thrds: t.join()
			
		print("Particles written to {}".format(options.output))
	
	

def make3d(jsd, ids, imgs, ttparams, pinfo, options, ctfinfo=[], tltkeep=[]):
	ppos, outname, boxsz=pinfo
	
	bx=boxsz*2
	pad=options.pad
	p3d=good_size(int(pad*1.5))
	apix=imgs[0]["apix_x"]
	if len(ctfinfo)>0:
		defocus, phase, voltage, cs=ctfinfo
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":voltage, "bfactor":50., "cs":cs,"ampcont":0, "apix":apix})
	
	for pid in ids:
		
		pos=ppos[pid]
		
		if type(pos)==type(Transform()):
			tf_dir=Transform(pos.get_rotation()).inverse()
			pos=np.array(pos.get_trans())
			
			
		elif len(pos)>3:
			drs=pos[3:].copy()
			pos=pos[:3].copy()
			drs/=np.linalg.norm(drs)
			drs=drs.tolist()
			tf_dir=Transform()
			tf_dir.set_rotation(drs)
			tf_dir.invert()

		else:
			tf_dir=None

		recon=Reconstructors.get("fourier", {"sym":'c1', "size":[p3d, p3d, p3d], "mode":"gauss_2"})
		recon.setup()
		projs=[]
		for nid, m in enumerate(imgs):
			if len(tltkeep)>0:
				if nid not in tltkeep:
					continue
			
			tpm=ttparams[nid]
			
			yt=tpm[3]-options.ytilt
			if abs(yt)>options.maxtilt: continue

			pxf=get_xf_pos(ttparams[nid], pos)

			tx=m["nx"]//2 +pxf[0]
			ty=m["ny"]//2 +pxf[1]

			txint=int(tx)
			tyint=int(ty)

			txdf=tx-txint
			tydf=ty-tyint

			e=m.get_clip(Region(txint-pad//2, tyint-pad//2, pad, pad), fill=0)
			
			e.process_inplace("mask.zeroedgefill",{"nonzero":1})		# This tries to deal with particles that were boxed off the edge of the micrograph
			if e.has_attr("hadzeroedge") and e["hadzeroedge"]!=0:
				continue

			e.mult(-1)
			e.process_inplace("normalize.edgemean")
			e.process_inplace("mask.soft",{"outer_radius":-1})
			if e["sigma"]==0:
				continue
			
			if options.rmbeadthr>0:
				er=e.process("filter.lowpass.gauss",{"cutoff_freq":0.01})
				er.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
				er.process_inplace("threshold.binary",{"value":options.rmbeadthr})
				er.process_inplace("mask.addshells.gauss",{"val1":e["nx"]//40,"val2":e["nx"]//25})
				
				if er["mean"]>0.1: 
					### too many beads. better just remove this subtilt
					continue
				
				e.mult(1-er)
				e.mult(1-er)
				noise=e.copy()
				noise.to_zero()
				noise.process_inplace("math.addnoise",{"noise":1})
				noise.process_inplace("normalize")
				e.add(noise*er)
				
			rot=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
			p1=rot.transform(pos.tolist())
			pz=p1[2]
			dz=pz*apix/10000.
			
			if len(ctfinfo)>0:
				df=defocus[nid]-dz
				ctf.set_phase(phase[nid]*np.pi/180.)
				ctf.defocus=df
				e["ctf"]=ctf
				fft1=e.do_fft()
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				if options.wiener:
					ctf.dsbg=1./(e["apix_x"]*e["nx"])
					s=np.array(ctf.compute_1d(e["nx"], ctf.dsbg, Ctf.CtfType.CTF_INTEN ))
					s[:np.argmax(s)]=np.max(s)
					s=np.maximum(s, 0.0001)
					ctf.snr=s.tolist()
					flt=fft1.copy_head()
					ctf.compute_2d_complex(flt,Ctf.CtfType.CTF_WIENER_FILTER)
					fft1.mult(flt)
				
				e=fft1.do_ift()
				e["ctf"]=ctf
				#if options.dotest:
					#print(tpm[3], pz, dz, defocus[nid]-dz)
			

			xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4], "ztilt":tpm[2], "tx":txdf, "ty":tydf})
			e["xform.projection"]=xform
			e["ptcl_source_coord"]=[float(txint), float(tyint), float(dz)]
			e["ptcl_source_src"]=options.tltfile
			e["model_id"]=pid
			e["tilt_id"]=nid
			e["file_threed"]=options.output
			e["ptcl_source_coord_3d"]=pos.tolist()
			e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.45})
			projs.append(e)
			
			sz=e["nx"]
			e0=e.get_clip(Region((sz-p3d)//2,(sz-p3d)//2,p3d,p3d))
			trans=Transform({"type":"2d", "tx":-txdf, "ty":-tydf})
			e1=recon.preprocess_slice(e0, trans)
			recon.insert_slice(e1,xform,1)

		if len(projs)<len(imgs)/5:
			#### too many bad 2D particles
			continue
		else:
			threed=recon.finish(True)
			threed.process_inplace("math.gausskernelfix",{"gauss_width":4.0})
			threed=threed.get_clip(Region((p3d-bx)//2,(p3d-bx)//2,(p3d-bx)//2,bx,bx,bx))
			
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix
		threed["ptcl_source_coord"]=pos.tolist()
		threed["file_twod"]=options.output2d
		threed["model_id"]=pid
		
		if tf_dir:
			threed["xform.align3d"]=tf_dir
		#print(ids,projs)
		jsd.put((pid, threed, projs))

	return
	


#### get 2D position on a tilt given 3D location
def get_xf_pos(tpm, pk):
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0], pk[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)
	
	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0], p1[1]])
	

	return [p2[0], p2[1]]


#### parse a json file from spt_align to re-extract particles
def parse_json(jsfile, xffile=""):
	js=js_open_dict(jsfile)
	coord=[]
	allxfs={}
	
	postxfs=[]
	if xffile!="":
		f=open(xffile,'r')
		lines=f.readlines()
		f.close()
		for l in lines:
			if len(l)>3:
				postxfs.append(Transform(eval(l)))
	else:
		postxfs.append(Transform())
	
	
	for ky in js.keys():
		src, ii = eval(ky)
		e=EMData(src, ii, True)
		tomo=e["class_ptcl_src"]
		dic=js[ky]
		dxf=dic["xform.align3d"]
		c=e["ptcl_source_coord"]
		for xf in postxfs:
			ali=Transform(dxf)
			ali=xf.inverse()*ali
			a=ali.inverse()
			a.translate(c[0], c[1], c[2])
			#print(a)
			
			if allxfs.has_key(tomo):
				allxfs[tomo].append(a)
			else:
				allxfs[tomo]=[a]
			
	js.close()
	#print(len(allxfs))
	return allxfs
	

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()


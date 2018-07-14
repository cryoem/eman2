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


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="Specify tomograms from which you wish to extract boxed particles.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="extract")
	parser.add_argument("--boxsz", type=int,help="box size in binned tomogram", default=-1, guitype='intbox',row=2, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--label", type=str,help="Only extract particle with this name. Leave blank to extract all particles.", default=None, guitype='strbox',row=2, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=3, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=4, col=1,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--maxtilt", type=int,help="max tilt", default=100, guitype='intbox',row=4, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--padtwod", type=float,help="padding factor", default=2.0, guitype='floatbox',row=5, col=0, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--noctf", action="store_true", default=False ,help="skip ctf correction..", guitype='boolbox',row=5, col=1, rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--dotest", action="store_true", default=False ,help="only make 1 batch of subtomograms for testing")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if len(args)==1:
		print("Reading from {}...".format(args[0]))
	else:
		print("Processing {} files in sequence..".format(len(args)))
		cmd=sys.argv
		opt=' '.join([s for s in cmd if s.startswith("-")])
		for a in args:
			run("{} {} {}".format(cmd[0], a, opt))
		return
	
	pfile=args[0]
	
	#### reading alignment info...
	js=js_open_dict(info_name(pfile))
	ttparams=np.array(js["tlt_params"])
	
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
		
	tfile=str(js["tlt_file"])
	
	ptcl=EMData(pfile, 0, True)
	img=EMData(tfile,0, True)
	apix_ptcl=ptcl["apix_x"]
	apix_tlt=img["apix_x"]
	zshift=ptcl["zshift"]#/2
	#yt=ptcl["ytilt"]
	yt=0
	#options.ytilt=ttparams[np.argmin(abs(ttparams[:,3]-yt)),3]
	options.ytilt=0
	scale=old_div(apix_ptcl,apix_tlt)

	print("Scaling factor: {:.1f}, y-tilt: {:.1f}, z-shift: {:d}.".format(scale, options.ytilt, int(zshift)))
	
	
	#### reading particle location from json file corresponding to the tomogram input
	ptclpos=[]
	nptcl=EMUtil.get_image_count(pfile)
	e=EMData(pfile, 0, True)
	if e["nz"]!=e["nx"] and nptcl==1:
		print("Reading particle location from a tomogram...")
		js=js_open_dict(info_name(pfile))
		towrite=[]
		if "class_list" in js and "boxes_3d" in js:
			clslst=js["class_list"]
			boxes=js["boxes_3d"]
			for ky in list(clslst.keys()):
				val=clslst[ky]
				if options.label:
					if str(val["name"])!=options.label:
						continue
						
				bxs=np.array([[b[0], b[1], b[2]] for b in boxes if b[5]==int(ky)], dtype=float)
				bxs-=[old_div(e["nx"],2), old_div(e["ny"],2), old_div(e["nz"],2)]
				bxs*=scale
				bxs[:,2]-=zshift
				outname=str(base_name(pfile)+"__"+val["name"]+".hdf")
				if options.boxsz<0:
					sz=int(val["boxsize"])*scale/2
				else:
					sz=int(np.round(options.boxsz/2.*scale))
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
		outname=base_name(pfile)+sfx
		
		
		ptclpos=np.array(ptclpos, dtype=float)
		
		ptclpos*=scale
		ptclpos[:,2]-=zshift#-2
		
		print("Reading {} particles".format(len(ptclpos)))
		
		if options.boxsz<0:
			boxsz=int(ptcl["nx"]/2*scale)
		else:
			boxsz=int(options.boxsz/2*scale)
		
		
		
		towrite=[(ptclpos,outname, boxsz)]
	
	
	print("Reading tilt series file: {}".format(tfile))
	#### make sure this works on image stack or mrc volume
	img=EMData(tfile,0)
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	for m in imgs: m.process_inplace("normalize")
	ntlt=len(imgs)
		
	try: os.mkdir("particles3d")
	except: pass
	try: os.mkdir("particles")
	except: pass

	for pinfo in towrite:
		ptclpos, outname,options.boxsz=pinfo
		if options.dotest:
			nptcl=1
		else:
			nptcl=len(ptclpos)
	
		options.output=os.path.join("particles3d", outname)
		options.output2d=os.path.join("particles", outname)
		options.pad=pad=good_size(options.boxsz*2*options.padtwod)
		
		print("Writing {} particles to {}".format(nptcl, outname))
		print("Box size {}, padded to {}".format(options.boxsz*2, pad))
		
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
			jobs.append([jsd, ids, imgs, ttparams, ptclpos, options, ctf])
		
		
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
				threed.write_image(options.output, pid)
				ndone+=1
				#if ndone%10==0:
				sys.stdout.write("\r{}/{} finished.".format(ndone, nptcl))
				sys.stdout.flush()

		for t in thrds: t.join()
			
		print("Particles written to {}".format(options.output))
	
	E2end(logid)
	

def make3d(jsd, ids, imgs, ttparams, ppos, options, ctfinfo=[]):
	
	
	bx=options.boxsz*2
	pad=options.pad
	p3d=good_size(int(pad*1.5))
	apix=imgs[0]["apix_x"]
	if len(ctfinfo)>0:
		defocus, phase, voltage, cs=ctfinfo
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":voltage, "bfactor":0., "cs":cs,"ampcont":0, "apix":apix})
	
	for pid in ids:
		
		pos=ppos[pid]
		recon=Reconstructors.get("fourier", {"sym":'c1', "size":[p3d, p3d, p3d]})
		recon.setup()
		projs=[]
		for nid, m in enumerate(imgs):
			tpm=ttparams[nid]
			
			yt=tpm[3]-options.ytilt
			if abs(yt)>options.maxtilt: continue

			pxf=get_xf_pos(ttparams[nid], pos)

			tx=old_div(m["nx"],2) +pxf[0]
			ty=old_div(m["ny"],2) +pxf[1]

			txint=int(tx)
			tyint=int(ty)

			txdf=tx-txint
			tydf=ty-tyint

			e=m.get_clip(Region(txint-old_div(pad,2), tyint-old_div(pad,2), pad, pad), fill=0)

			e.mult(-1)
			e.process_inplace("normalize.edgemean")
			e.process_inplace("mask.soft",{"outer_radius":-1})
			if e["sigma"]==0:
				continue
			
			if len(ctfinfo)>0:
				rot=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
				p1=rot.transform(pos.tolist())
				pz=p1[2]
				dz=pz*apix/10000.
				df=defocus[nid]-dz
				ctf.set_phase(phase[nid]*np.pi/180.)
				ctf.defocus=df
				e["ctf"]=ctf
				fft1=e.do_fft()
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				e=fft1.do_ift()
				
				if options.dotest:
					print(tpm[3], pz, dz, defocus[nid]-dz)
			

			xform=Transform({"type":"xyz","ytilt":tpm[3],"xtilt":tpm[4], "ztilt":tpm[2], "tx":txdf, "ty":tydf})
			e["xform.projection"]=xform
			e["model_id"]=pid
			e["tilt_id"]=nid
			e["file_threed"]=options.output
			projs.append(e)
			
			sz=e["nx"]
			e0=e.get_clip(Region(old_div((sz-p3d),2),old_div((sz-p3d),2),p3d,p3d))
			trans=Transform({"type":"2d", "tx":-txdf, "ty":-tydf})
			e1=recon.preprocess_slice(e0, trans)
			recon.insert_slice(e1,xform,1)

		threed=recon.finish(True)
		threed=threed.get_clip(Region(old_div((p3d-bx),2),old_div((p3d-bx),2),old_div((p3d-bx),2),bx,bx,bx))
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix
		threed["ptcl_source_coord"]=pos.tolist()
		threed["file_twod"]=options.output2d
		
		jsd.put((pid, threed, projs))

	return threed
	


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


def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()


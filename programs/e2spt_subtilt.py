#!/usr/bin/env python
# Muyuan Chen 2018-04
from EMAN2 import *
import numpy as np
import Queue
import threading


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particles",help="Specify particles on which you want to perform sub-tilt refinement.", default="", guitype='filebox', browser="EMSPTParticleTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="model")
	parser.add_header(name="orblock1", help='Just a visual separation', title="** Subtilt Refinement Options **", row=2, col=0, rowspan=1, colspan=1, mode="model")
	parser.add_argument("--output", type=str,help="output file name", default=None, guitype='strbox',row=4, col=0,rowspan=1, colspan=2, mode="model")
	parser.add_argument("--boxsz", type=int,help="box size in binned tomogram", default=1, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--threads", type=int,help="threads", default=1, guitype='intbox',row=6, col=1,rowspan=1, colspan=1, mode="model")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	pfile=args[0]
	
	ptclpos=[]
	nptcl=EMUtil.get_image_count(pfile)
	for i in range(nptcl):
		ptcl=EMData(pfile, i, True)
		ptclpos.append(ptcl["ptcl_source_coord"])
	ptclpos=np.array(ptclpos, dtype=float)
	
	js=js_open_dict(info_name(pfile))
	ttparams=np.array(js["tlt_params"])
	ttparams[:,:2]*=2
	tfile=str(js["tlt_file"])
	print("Reading tilt series file: {}".format(tfile))
	   
	img=EMData(tfile,0)
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	apix_ptcl=ptcl["apix_x"]
	apix_tlt=imgs[0]["apix_x"]
	zshift=ptcl["zshift"]
	scale=apix_ptcl/apix_tlt
	ptclpos*=scale
	ptclpos[:,2]-=zshift
	print("Scaling factor: {:.1f}, z-shift: {:d}.".format(scale, int(zshift)))
	
	if options.boxsz<0:
		options.boxsz=int(ptcl["nx"]/2*scale)
	else:
		options.boxsz=int(options.boxsz/2*scale)
	
	k=0
	
	if options.output==None:
		sfx=pfile[pfile.find("__"):]
		for i in range(2,10): sfx=sfx.replace("_bin{:d}".format(i),"")
		options.output=os.path.join("particles3d", base_name(pfile)+sfx)
		try: os.mkdir("particles3d")
		except: pass
	try: os.remove(options.output)
	except: pass

	jsd=Queue.Queue(0)
	jobs=[]
	
	batchsz=10
	for tid in range(0,nptcl,batchsz):
		ids=range(tid, min(tid+batchsz, nptcl))
		jobs.append([jsd, ids, imgs, ttparams, ptclpos, options])
		
	thrds=[threading.Thread(target=make3d,args=(i)) for i in jobs]
	thrtolaunch=0
	tsleep=threading.active_count()
	
	ndone=0
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads+1 ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
		
		
		while not jsd.empty():
			pid, threed=jsd.get()
			threed.write_image(options.output, pid)
			ndone+=1
			if ndone%10==0:
				print("{}/{} finished.".format(ndone, nptcl))

	for t in thrds: t.join()
	
		
	print("Particles written to {}".format(options.output))
	E2end(logid)
	

def make3d(jsd, ids, imgs, ttparams, ppos, options):
	
	
	bx=options.boxsz
	pad=good_size(bx*4)
	for pid in ids:
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad]})
		recon.setup()

		for nid in range(len(imgs)):

			tpm=ttparams[nid]

			pxf=get_xf_pos(ttparams[nid], ppos[pid])

			pxf[0]+=imgs[nid]["nx"]/2
			pxf[1]+=imgs[nid]["ny"]/2

			e=imgs[nid].get_clip(Region(pxf[0]-pad/2,pxf[1]-pad/2, pad, pad))
			e.process_inplace("normalize")
			e.process_inplace("threshold.clampminmax.nsigma",{"nsigma":5})
			rot=Transform({"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4]})
			e1=recon.preprocess_slice(e, rot)
			recon.insert_slice(e1,rot,1)

		threed=recon.finish(True)
		threed=threed.get_clip(Region((pad-bx*2)/2,(pad-bx*2)/2,(pad-bx*2)/2,bx*2,bx*2,bx*2))
		threed.process_inplace("normalize.edgemean")
		#threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
		threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=imgs[0]["apix_x"]
		threed.process_inplace("mask.soft",{"outer_radius":-1})
		threed.mult(-1)
		jsd.put((pid, threed))
		#print(pid)
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


def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
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
import scipy.spatial.distance as scidist
import queue
import threading

def main():
	
	usage="A simple template matching script. run [prog] <tomogram> <reference> to extract particles from tomogram. Results will be saved in the corresponding info files and can be visualized via spt_boxer"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="Specify tomograms containing reference-like particles to be exctracted.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="boxing")
	
	parser.add_argument("--reference",help="Specify a 3D reference volume.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=2, mode="boxing")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=0, rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--label", type=str,help="Assign unique label to particles resembling specified reference. This allows specific particles to be extracted in the next step and aids project organization with easily interpreted filenames.\nIf --label is not specified, this set of particles will be labeled according to the file name of the reference without file extension.", default=None, guitype='strbox',row=3, col=0, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--nptcl", type=int,help="maximum number of particles", default=500, guitype='intbox', row=3, col=1,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--dthr", type=float,help="distance threshold", default=-1, guitype='floatbox', row=4, col=0,rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--vthr", type=float,help="value threshold (n sigma)", default=2.0, guitype='floatbox', row=4, col=1,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--delta", type=float,help="delta angle", default=30.0, guitype='floatbox', row=5, col=0,rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox', row=5, col=1,rowspan=1, colspan=1, mode="boxing")
	
	parser.add_argument("--rmedge", action="store_true",help="Remove particles on the edge.", default=False, guitype='boolbox', row=6, col=0,rowspan=1, colspan=1, mode="boxing[True]")
	parser.add_argument("--rmgold", action="store_true",help="Remove particles near gold fiducial.", default=False, guitype='boolbox', row=6, col=1,rowspan=1, colspan=1, mode="boxing[True]")

	parser.add_argument("--boxsz", type=int,help="Overwrite box size", default=-1, guitype='intbox', row=7, col=0,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--ppid", type=int,help="ppid", default=-2)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	time0=time.time()

	tmpname=options.reference #args[1]

	sym=parsesym("c1")
	dt=options.delta
	oris=sym.gen_orientations("eman",{"delta":dt, "phitoo":dt})
	print("Try {} orientations.".format(len(oris)))
	
	

	for filenum,imgname in enumerate(args):
		
		print("Locating reference-like particles in {} (File {}/{})".format(imgname,filenum+1,len(args)))
		img=EMData(imgname)
		nbin=int(img["nx"]//500)
		print("Will shrink tomogram by {}".format(nbin))
		img.process_inplace("math.meanshrink",{'n':nbin})
		tomo=img.copy()
		img.mult(-1)
		img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.4})
		img.process_inplace('normalize')
		img.process_inplace('threshold.clampminmax.nsigma', {"nsigma":1})
		
		m=EMData(tmpname)
		mbin=img["apix_x"]/m["apix_x"]
		print("Will shrink reference by {:.1f}".format(mbin))
		m.process_inplace("math.fft.resample",{'n':mbin})
		m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
		m.process_inplace("mask.soft",{"outer_radius":-2})
		m.process_inplace('normalize')
		sz=m["nx"]
		if options.dthr<0:
			options.dthr=sz/np.sqrt(2)

		hdr=m.get_attr_dict()
		ccc=img.copy()*0-65535

		jsd=queue.Queue(0)
		thrds=[threading.Thread(target=do_match,args=(jsd, m,o, img)) for o in oris]
		thrtolaunch=0
		tsleep=threading.active_count()

		ndone=0
		while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==13 ) : time.sleep(.1)
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(1)

			while not jsd.empty():
				cf=jsd.get()
				ccc.process_inplace("math.max", {"with":cf})
				ndone+=1
				#if ndone%10==0:
				sys.stdout.write("\r{}/{} finished.".format(ndone, len(oris)))
				sys.stdout.flush()
		print("")

		cbin=ccc.process("math.maxshrink", {"n":2})
		msk=cbin.copy()
		msk.to_one()
		
		if options.rmedge:
			try:
				js=js_open_dict(info_name(fname))
				tpm=np.array(js["tlt_params"])
				js.close()
				rt=np.mean(tpm[:,2])
			except:
				rt=0

			eg=16
			msk.process_inplace("mask.zeroedge3d",{"x0":eg,"x1":eg,"y0":eg,"y1":eg,"z0":eg,"z1":eg})
			msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			msk.rotate(0,0,-rt)
			
			cbin.mult(msk)
		
		if options.rmgold:
			tomo.process_inplace("math.meanshrink",{"n":2})
			tomo.process_inplace("normalize")
			tomo.mult(-1)
			tomo.process_inplace("threshold.binary",{"value":10})
			tomo.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
			tomo.process_inplace("normalize.edgemean")
			tomo=1-tomo
			
			cbin.mult(tomo)
		
		cbin.process_inplace("threshold.belowtozero")
		cbin.process_inplace("normalize.edgemean")
		cbin.write_image("ccc.hdf")
		#cbin.write_image("tmp0.hdf")
		#msk.write_image("tmp1.hdf")
		cc=cbin.numpy().copy()
		cshp=cc.shape
		ccf=cc.flatten()
		asrt= np.argsort(-ccf)
		pts=[]
		vthr=np.mean(ccf)+np.std(ccf)*options.vthr
		
		dthr=options.dthr/4.
		scr=[]
		#print vthr,cc.shape
		for i in range(len(asrt)):
			aid=asrt[i]
			pt=np.unravel_index(aid, cshp)
			if len(pts)>0:
				dst=scidist.cdist(pts, [pt])
				if np.min(dst)<dthr:
					continue

			pts.append(pt)
			scr.append(float(ccf[aid]))
			if cc[pt]<vthr or len(pts)>=options.nptcl:
				break
				
		pts=np.array(pts)
		print("Found {} particles".format(len(pts)))
		js=js_open_dict(info_name(imgname))
		n=min(options.nptcl, len(pts))
		
		if "class_list" in js:
			clst=js['class_list']
			try: kid=max([int(k) for k in list(clst.keys())])+1
			except: kid=0 # In case someone manually edited the info file. Unlikely.
		else:
			clst={}
			kid=0
		
		if "boxes_3d" in js:
			bxs=js["boxes_3d"]
		else:
			bxs=[]
			
		e=img
		if "apix_unbin" in js:
			apix_unbin=js["apix_unbin"]
			apix=e["apix_x"]
			shp=np.array([e["nz"], e["ny"], e["nx"]])
			#print(sz,nbin,apix, apix_unbin)
			if options.boxsz<0:
				boxsz=int(np.round(sz*apix/apix_unbin))
			else:
				boxsz=options.boxsz
			
			box=(pts*2-shp/2)*apix/apix_unbin
			bxs.extend([[p[2], p[1],p[0], 'tm', scr[i] ,kid] for i,p in enumerate(box[:n])])
			
		else:
			bxs.extend([[p[2], p[1],p[0], 'tm', scr[i] ,kid] for i,p in enumerate(pts[:n]*4)])
			if options.boxsz<0:
				boxsz=sz*nbin
			else:
				boxsz=options.boxsz
				
			
		js['boxes_3d']=bxs
		
		if options.label:
			clst[str(kid)]={"boxsize":boxsz, "name":options.label}
		else:
			clst[str(kid)]={"boxsize":boxsz, "name":base_name(tmpname)}
			
		js["class_list"]=clst
		js.close()

	E2end(logid)
	
def do_match(jsd, m, o, img):

	for i in [0,1]:
		e=m.copy()
		e.transform(o)
		cf=img.calc_ccf(e)
		cf.process_inplace("xform.phaseorigin.tocenter")
		jsd.put(cf)
		o.rotate(Transform({"type":"eman","alt":180}))
		
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()


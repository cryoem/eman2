#!/usr/bin/env python
# Muyuan Chen 2018-04

from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
import scipy.spatial.distance as scidist
import queue
import threading
from scipy import ndimage
from scipy.spatial import KDTree

def main():
	
	usage="A simple template matching script. run [prog] <tomogram> <reference> to extract particles from tomogram. Results will be saved in the corresponding info files and can be visualized via spt_boxer"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="Specify tomograms containing reference-like particles to be exctracted. These should be 'dark contrast'", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="boxing")
	
	parser.add_argument("--reference",help="Specify a 3D reference volume. This should be 'light contrast', ie - positive isosurface.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=2, mode="boxing")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=0, rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--label", type=str,help="Assign unique label to particles resembling specified reference. This allows specific particles to be extracted in the next step and aids project organization with easily interpreted filenames.\nIf --label is not specified, this set of particles will be labeled according to the file name of the reference without file extension.", default=None, guitype='strbox',row=3, col=0, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--nptcl", type=int,help="maximum number of particles per tomogram", default=500, guitype='intbox', row=3, col=1,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--dthr", type=float,help="The program will remove particles closer than this distance threshold. By default, this will be 1/2 of box size or the reference. Otherwise, specify the distance in Angstrom.", default=-1, guitype='floatbox', row=4, col=0,rowspan=1, colspan=1, mode="boxing")
	
	parser.add_argument("--vthr", type=float,help="template matching value threshold (n sigma). Particles with score lower than this will be removed.", default=10, guitype='floatbox', row=4, col=1,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--delta", type=float,help="Anglular sampling to rotate the reference.", default=30.0, guitype='floatbox', row=5, col=0,rowspan=1, colspan=1, mode="boxing")
	
	parser.add_argument("--sym", type=str,help="Symmetry of reference.", default="c1", guitype='strbox', row=5, col=1,rowspan=1, colspan=1, mode="boxing")
	
	parser.add_argument("--rmedge", action="store_true",help="Remove particles on the edge.", default=False, guitype='boolbox', row=6, col=0,rowspan=1, colspan=1, mode="boxing[True]")
	parser.add_argument("--rmgold", action="store_true",help="Remove particles near gold fiducial.", default=False, guitype='boolbox', row=6, col=1,rowspan=1, colspan=1, mode="boxing[True]")

	parser.add_argument("--boxsz", type=int,help="Overwrite box size of the reference. This should be the box size of unbinned micrographs if specified.", default=-1, guitype='intbox', row=7, col=0,rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--threads", type=int,help="number of threads to use", default=12, guitype='intbox', row=8, col=0,rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--shrink", type=int,help="binning factor. Default (-1) will downsample the tomograms to ~500px for template matching", default=-1)
	parser.add_argument("--ppid", type=int,help="ppid", default=-2)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	time0=time.time()

	tmpname=options.reference #args[1]

	sym=parsesym(options.sym)
	dt=options.delta
	oris=sym.gen_orientations("eman",{"delta":dt, "phitoo":dt,"inc_mirror":1})
	print("Testing {} orientations...".format(len(oris)))
	
	def do_match(jsd, xfs):
		for xf in xfs:
			r=ref.process("xform", {"transform":xf})
			cf=tomo.calc_ccf(r)
			jsd.put(cf)

	for filenum,imgname in enumerate(args):
		
		print("Locating reference-like particles in {} (File {}/{})".format(imgname,filenum+1,len(args)))
		tomo=EMData(imgname)
		if options.shrink>0:
			nbin=options.shrink
		else:
			nbin=int(tomo["nx"]//450)
		if nbin>1:
			print("Will shrink tomogram by {}".format(nbin))
			tomo.process_inplace("math.meanshrink",{'n':nbin})
			
			
		ref=EMData(tmpname)
		mbin=tomo["apix_x"]/ref["apix_x"]
		print("Will shrink reference by {:.1f}".format(mbin))
		ref.process_inplace("math.fft.resample",{'n':mbin})
		ref.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.4})
		ref.process_inplace("filter.highpass.gauss",{"cutoff_pixels":4})
		ref.process_inplace('normalize.edgemean')
		ref.mult(-1)
		
		boxsz=ref["ny"]//2
		tomo.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.4})
		tomo.process_inplace("filter.highpass.gauss",{"cutoff_pixels":boxsz})
		tomo.process_inplace('normalize.edgemean')
		tomo0=tomo.copy()
		tomo.process_inplace('threshold.clampminmax.nsigma', {"nsigma":3})
		
		
		nthrds=options.threads
		rg=np.arange(len(oris))
		tasks=[oris[i::nthrds] for i in range(nthrds)]

		ccc=tomo.copy()*0-65535

		jsd=queue.Queue(0)
		thrds=[threading.Thread(target=do_match,args=(jsd, o)) for o in tasks]
		thrtolaunch=0
		tsleep=threading.active_count()

		ndone=0
		for t in thrds:
			t.start()
		
		while threading.active_count()>tsleep or not jsd.empty():
		
			while not jsd.empty():
				cf=jsd.get()
				ccc.process_inplace("math.max", {"with":cf})
				ndone+=1
				sys.stdout.write("\r{}/{} finished.".format(ndone, len(oris)))
				sys.stdout.flush()
			
			time.sleep(.5)
		print("")
		
		ccc.process_inplace("xform.phaseorigin.tocenter")
		ccc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		ccc.process_inplace("filter.highpass.gauss",{"cutoff_pixels":24})
		ccc.process_inplace("normalize.edgemean")
		ccc.process_inplace("threshold.belowtozero")

		
		if options.rmedge:
			try:
				js=js_open_dict(info_name(imgname))
				tpm=np.array(js["tlt_params"])
				js.close()
				rt=np.mean(tpm[:,2])
				print("removing edge: rotation {:.1f}".format(rt))
			except:
				rt=0

			msk=ccc.copy()
			msk.to_one()
			eg=32
			msk.process_inplace("mask.zeroedge3d",{"x0":eg,"x1":eg,"y0":eg,"y1":eg})
			msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
			msk.process_inplace("xform", {"phi":90-rt})
			ccc.mult(msk)
		
		if options.rmgold:
			
			img=tomo0.copy()
			img.mult(-1)
			img.process_inplace("threshold.binary",{"value":8})
			img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
			img.process_inplace("normalize.edgemean")
			img=1-img
			ccc.mult(img)
		
		
		ccc.write_image("tmp_ccc.hdf")
		img=ccc.numpy().copy()
		img[img<options.vthr]=0
		lb, nlb=ndimage.measurements.label(img)
		pks=np.array(ndimage.maximum_position(img,lb,list(range(1,nlb))))
		#pks=np.array(ndimage.center_of_mass(img,lb,list(range(1,nlb))))
		pksize=np.array(ndimage.measurements.sum(img,lb,list(range(1,nlb))))
		n=len(pks)
		print(len(pks))
		
		#### filter out small peaks
		if options.boxsz>0:
			boxsz=options.boxsz
			
		kpid=pksize>5
		pks=pks[kpid]
		pkscore=pksize[kpid]
		
		srt=np.argsort(-pkscore)
		pks=pks[srt]
		pkscore=pkscore[srt]
		pkscore/=np.max(pkscore)
		
		
		tree=KDTree(pks)

		tokeep=np.ones(len(pks), dtype=bool)
		if options.dthr>0:
			dthr=options.dthr/tomo["apix_x"]
		else:
			dthr=boxsz
			
		print(np.min(pks, axis=0), np.max(pks, axis=0), boxsz, dthr)
		for i in range(len(pks)):
			if tokeep[i]:
				k=tree.query_ball_point(pks[i], dthr)
				tokeep[k]=False
				tokeep[i]=True
			
		pts=pks[tokeep]
		scr=pkscore[tokeep]
		
		if len(pts)>options.nptcl:
			print("Found {} particles. Keep the best {}.".format(len(pts), options.nptcl))
			pts=pts[:options.nptcl]
		
		else:
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
			
		if "apix_unbin" in js:
			apix_unbin=js["apix_unbin"]
			apix=tomo["apix_x"]
			shp=np.array([tomo["nz"], tomo["ny"], tomo["nx"]])
			#print(sz,nbin,apix, apix_unbin)
			if options.boxsz<0:
				boxsz=int(np.round(boxsz*apix/apix_unbin*2))
			else:
				boxsz=options.boxsz
			
			box=(pts-shp/2)*apix/apix_unbin
			bxs.extend([[p[2], p[1],p[0], 'tm', scr[i] ,kid] for i,p in enumerate(box[:n])])
			
		else:
			bxs.extend([[p[2], p[1],p[0], 'tm', scr[i] ,kid] for i,p in enumerate(pts[:n]*nbin)])
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
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()


#!/usr/bin/env python
# Muyuan Chen 2016-09
from past.utils import old_div
from builtins import range
from EMAN2 import *
import numpy as np
from scipy import ndimage

def main():
	
	usage="This program is designed to extract subtomograms from segmentation results from the Tomoseg workflow in EMAN2. Run e2spt_boxer22 to visualize the results, and follow the SPT protocol for following steps."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomogram",help="Specify a tomogram from which you want to extract particles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='tomograms')", row=0, col=0,rowspan=1, colspan=2, mode="extract")
	parser.add_pos_argument(name="segmentation",help="The corresponding segmentation of the tomogram.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='segmentations')", row=1, col=0,rowspan=1, colspan=2, mode="extract")

	parser.add_argument("--thresh", type=float,help="Threshold of density value for particle extraction.", default=1.0, guitype='floatbox', row=2, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--featurename", type=str,help="name of the current feature to extract", default="", guitype='strbox', row=2, col=1,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--boxsz", type=int,help="Box size", default=32, guitype='intbox', row=3, col=0,rowspan=1, colspan=1, mode="extract")
	
	parser.add_argument("--random", type=int,help="Specifying N will randomly seed N particles on density above threshold. default is -1, means only choosing peaks. Useful for non-globular particles", default=-1, guitype='intbox', row=4, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--massthresh", type=float,help="Threshold of total mass of each continous object to be considered a particle. ", default=20., guitype='floatbox', row=4, col=1,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--edge", type=int,help="mininmum distance to the edge of tomogram", default=4, guitype='intbox', row=5, col=0,rowspan=1, colspan=1, mode="extract")
	parser.add_argument("--sort", action="store_true",help="Sort particles by segmentation density.", default=False, guitype='boolbox', row=5, col=1,rowspan=1, colspan=1, mode="extract[True]")
	parser.add_argument("--ppid", type=int,help="ppid", default=-1)
	
	
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	tomoname=args[0]
	segname=args[1]
	
	e=EMData(segname)
	img=e.numpy()
	img[img<options.thresh]=0
	
	#if options.shrink==0:
		#tm=EMData(tomoname,0,True)
		#shrinkz=old_div(float(tm["nz"]),e["nz"])
		#shrinkxy=old_div(tm["nx"],e["nx"])
		#print("Shrink by {} in x-y plane, and shrink {} in z axis".format(shrinkxy, shrinkz))
	#else:
		#shrinkz=shrinkxy=options.shrink
	
	if options.random<=0:
		
		lb, nlb=ndimage.measurements.label(img)
		pks=np.array(ndimage.maximum_position(img,lb,list(range(1,nlb))))
		#pks=np.array(ndimage.center_of_mass(img,lb,list(range(1,nlb))))
		pksize=np.array(ndimage.measurements.sum(img,lb,list(range(1,nlb))))
		n=len(pks)
		
		#### filter out small peaks
		kpid=pksize>options.massthresh
		pks=pks[kpid]
		pksize=pksize[kpid]
		
		pkscore=[]
		pk_new=[[-100,-100,-100]]
		for ip,p in enumerate(pks):
			nb=np.sum(np.sum(np.array(pk_new-p)**2,axis=1)<(options.boxsz/4)**2)
			
			if nb<1:
				pk_new.append(p)
				pkscore.append(pksize[ip])
		pks=np.array(pk_new)[1:]
		pkscore=np.array(pkscore)
	
	else:
		pts=np.array(np.where(img>0)).T
		ip=np.arange(len(pts))
		np.random.shuffle(ip)
		pts=pts[ip]
		pks=pts[:options.random].copy()
		pkscore=np.array([1.0 for p in pks])
		
	
	n=len(pks)
	if n==0:
		print("No particles found. Exit...")
		return
	allbox=[]
	allbox3d=[]
	if options.sort:
		srt=np.argsort(-pkscore)
		pks=pks[srt]
		pkscore=pkscore[srt]
		
	pkscore/=np.max(pkscore)

	jsname=info_name(tomoname)
	js=js_open_dict(jsname)
	if "class_list" in js:
		clst=[int(k) for k in list(js["class_list"].keys())]
		for ii in range(100):
			if ii not in clst:
				mytag=ii
				break
			
	else:
		mytag=0
	
	if options.featurename=="":
		tomotag=tomoname[tomoname.rfind('__')+2:-4]
		segtag=segname[segname.rfind('__')+2:-4]
		segtag=[s for s in segtag.split('_') if s not in tomotag.split('_')]
		#print(segtag)
		options.featurename='_'.join(segtag)
	
	if "apix_unbin" in js:
		apix_unbin=js["apix_unbin"]
		apix=e["apix_x"]
		
	shp=np.array([e["nz"], e["ny"], e["nx"]])
	
	for j in range(n):
		box=pks[j].copy()
		if min(box[2],box[1],e["nx"]-box[2],e["ny"]-box[1])<options.edge+options.boxsz//2:
			continue
		
		if "apix_unbin" in js:
			box=(box-shp/2)*apix/apix_unbin
		allbox3d.append([box[2], box[1],int(box[0]), "tomoseg", float(np.round(pkscore[j],3)), mytag])
	
	if "boxes_3d" in js:
		b3d=js["boxes_3d"]
		b3d.extend(allbox3d)
	else:
		b3d=allbox3d
	
	js["boxes_3d"]=b3d
	if "class_list" in js:
		clst=js["class_list"]
	else:
		clst={}
		
	if "apix_unbin" in js:
		options.boxsz=options.boxsz*apix/apix_unbin
	clst[int(mytag)]={"name":options.featurename, "boxsize":options.boxsz}
	js["class_list"]=clst
	#js["boxes"]=allbox
	js.close()
	print("{} boxes found..".format(len(allbox3d)))

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	

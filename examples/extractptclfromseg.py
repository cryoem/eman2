#!/usr/bin/env python
# Muyuan Chen 2016-09
from EMAN2 import *
import numpy as np
havescipy=True
try: 
	from scipy import ndimage
except:
	havescipy=False
def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--thresh", type=float,help="Threshold of density value for particle extraction.", default=1.0)
	parser.add_argument("--massthresh", type=float,help="Threshold of total mass of each continous object to be considered a particle. ", default=20.)
	parser.add_argument("--edge", type=int,help="min distance to edge", default=4)
	parser.add_argument("--featurename", type=str,help="name of the current feature to extract", default=None)
	parser.add_argument("--sort", action="store_true",help="sort by density", default=False)
	parser.add_argument("--random", type=int,help="randomly seed particles on density above threshold", default=-1)
	parser.add_argument("--genptcls", type=str,help="generate particles", default=None)
	parser.add_argument("--genmask", type=str,help="generate mask", default=None)
	parser.add_argument("--boxsz", type=int,help="box size", default=48)
	parser.add_argument("--shrink", type=float,help="shrink factor", default=0)
	parser.add_argument("--zthick", type=int,help="make projection of this thickness", default=-1)
	parser.add_argument("--apix", type=float,help="apix", default=-1.)
	
	
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.genptcls or options.genmask:
		rawname=args[0]
		tomoname=args[1]
		raw=EMData(rawname,0,True)
		tomo=EMData(tomoname,0,True)
		if options.shrink>0:
			shrink=options.shrink
		else:
			shrink=raw["nx"]/tomo["nx"]
		
		if options.apix<=0:
			options.apix=tomo["apix_x"]
			print "Reading apix from data: {}".format(options.apix)
			
			
		nn=EMUtil.get_image_count(tomoname)
		if nn>1:
			print "Particle stack? getting coordinates from header.."
			pks=[]
			for i in range(nn):
				e=EMData(tomoname, i, True)
				pks.append(e["box"])
			pks=np.array(pks)*shrink
			
		else:
			jsname=info_name(tomoname)
			js=js_open_dict(jsname)
			pks=(np.array([[j[0],j[1],j[3]] for j in js["boxes"]])*shrink).astype(int)
			js=None
		
		bxsz=int(options.boxsz*shrink)
		b2=bxsz/2
		
		if options.zthick>0:
			print "Making projection of {} pixel thickness".format(options.zthick)
			zthick=options.zthick
		else:
			zthick=bxsz
		
		if options.genptcls:
			
			pname=options.genptcls
			
			try: os.remove(pname)
			except: pass
			print len(pks), " particles, unbin by ", shrink

			for p in pks:
				
				pj=EMData(rawname, 0, False,Region(p[0]-b2,p[1]-b2,p[2]-zthick/2, bxsz, bxsz, zthick))
				
				
				pj.process_inplace("normalize")
				pj.mult(-1)
				pj["apix_x"]=pj["apix_y"]=pj["apix_z"]=options.apix/shrink
				if options.zthick>0:
					pj=pj.project("standard", Transform())
				
				pj["box"]=[p[0], p[1], p[2]]
				pj["src"]=rawname
				pj.write_image(pname, -1)
		else:
			pname=options.genmask
			print len(pks), " particles, unbin by ", shrink
			try: os.remove(pname)
			except: pass
			
			e=EMData(raw["nx"], raw["ny"], raw["nz"])
			a=EMData(bxsz,bxsz,bxsz)
			a.to_one()
			a.process_inplace("mask.soft",{"outer_radius":-5})
			e.to_zero()
			
			for ii,p in enumerate(pks):
				if ii%100==0: print ii
				print ii,p
				e.insert_scaled_sum(a, p.tolist())
			
			e.process_inplace("threshold.clampminmax",{"maxval":1, "minval":0})
			e.write_image(pname)
			
	
	else:
		segname=args[0]
		tomoname=args[1]
		
		e=EMData(segname)
		img=e.numpy()
		img[img<options.thresh]=0
		
		if options.shrink==0:
			tm=EMData(tomoname,0,True)
			shrinkz=tm["nz"]/e["nz"]
			shrinkxy=tm["nx"]/e["nx"]
			print "Shrink by {} in x-y plane, and shrink {} in z axis".format(shrinkxy, shrinkz)
		else:
			shrinkz=shrinkxy=options.shrink
		
		if options.random<=0:
			if havescipy:
				lb, nlb=ndimage.measurements.label(img)
				pks=np.array(ndimage.maximum_position(img,lb,range(1,nlb)))
				pksize=np.array(ndimage.measurements.sum(img,lb,range(1,nlb)))
				n=len(pks)
			
			else:
				e.process_inplace("mask.onlypeaks")
				
				#print np.sum(img>0)
				pks= np.array(np.where(img>0)).T
				pksize=np.zeros(len(pks))+options.massthresh+1 
			
			#### filter out small peaks
			kpid=pksize>options.massthresh
			pks=pks[kpid]
			pksize=pksize[kpid]
			
			pkscore=[]
			pk_new=[[-100,-100,-100]]
			for ip,p in enumerate(pks):
				nb=np.sum(np.sum(np.array(pk_new-p)**2,axis=1)<(options.boxsz/4)**2)
				#print p, nb
				if nb<1:
					pk_new.append(p)
					pkscore.append(pksize[ip])
			pks=np.array(pk_new)[1:]
			pkscore=np.array(pkscore)
		
		else:
			pts=np.array(np.where(img>0)).T
			ip=np.arange(len(pts))
			#print len(pts)
			np.random.shuffle(ip)
			pts=pts[ip]
			pks=pts[:options.random].copy()
			pkscore=[0 for p in pks]
			#print pks
		
		n=len(pks)
		#e.write_image("tmp2.hdf")
		allbox=[]
		allbox3d=[]
		if options.sort:
			srt=np.argsort(-pkscore)
			pks=pks[srt]
			pkscore=pkscore[srt]
			#denmap=EMData(tomoname)
			#dmp=denmap.numpy()
			#den=np.array(ndimage.measurements.maximum(np.abs(dmp),lb,range(1,nlb)))
			#pks=pks[]
			#

		jsname=info_name(tomoname)
		js=js_open_dict(jsname)
		if js.has_key("class_list"):
			clst=[int(k) for k in js["class_list"].keys()]
			for ii in range(100):
				if ii not in clst:
					mytag=ii
					break
				
		else:
			mytag=0
		
		if options.featurename==None:
			options.featurename="feature_{:02d}".format(mytag)
		
		for j in range(n):
			box=pks[j]
			if min(box[2],box[1],e["nx"]-box[2],e["ny"]-box[1])<options.edge:
				continue
			#box*=2
			box[2]*=shrinkxy
			box[1]*=shrinkxy
			box[0]*=shrinkz
			allbox3d.append([box[2], box[1],int(box[0]), "tomoseg", pkscore[j], mytag])
			allbox.append([box[2], box[1],"manual", int(box[0])])
		
		if js.has_key("boxes_3d"):
			b3d=js["boxes_3d"]
			b3d.extend(allbox3d)
		else:
			b3d=allbox3d
		
		js["boxes_3d"]=b3d
		if js.has_key("class_list"):
			clst=js["class_list"]
		else:
			clst={}
		clst[int(mytag)]={"name":options.featurename, "boxsize":options.boxsz}
		js["class_list"]=clst
		js["boxes"]=allbox
		js.close()
		print "{} boxes found..".format(len(allbox))

	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
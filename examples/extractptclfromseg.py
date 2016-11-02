#!/usr/bin/env python
# Muyuan Chen 2016-09
from EMAN2 import *
import numpy as np
#from scipy import ndimage
def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--thresh", type=float,help="threshold", default=.0)
	parser.add_argument("--edge", type=int,help="min distance to edge", default=0)
	parser.add_argument("--sort", action="store_true",help="sort by density", default=False)
	parser.add_argument("--genptcls", type=str,help="generate particles", default=None)
	parser.add_argument("--genmask", type=str,help="generate mask", default=None)
	parser.add_argument("--boxsz", type=int,help="box size", default=48)
	parser.add_argument("--shrink", type=float,help="shrink", default=0)
	parser.add_argument("--apix", type=float,help="apix", default=1.0)
	
	
	
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
		
		jsname=info_name(tomoname)
		js=js_open_dict(jsname)
		pks=(np.array([[j[0],j[1],j[3]] for j in js["boxes"]])*shrink).astype(int)
		js=None
		
		bxsz=options.boxsz*shrink
		b2=bxsz/2
		
		if options.genptcls:
			
			pname=options.genptcls
			
			try: os.remove(pname)
			except: pass
			print len(pks), " particles, unbin by ", shrink

			for p in pks:
				
				pj=EMData(rawname, 0, False,Region(p[0]-b2,p[1]-b2,p[2]-b2, bxsz, bxsz,bxsz))
				
				pj["box"]=[p[0], p[1], p[2]]
				pj["src"]=rawname
				pj.process_inplace("normalize")
				pj.mult(-1)
				pj["apix_x"]=pj["apix_y"]=pj["apix_z"]=options.apix/shrink
				pj.write_image(pname, -1)
		else:
			pname=options.genmask
			
			try: os.remove(pname)
			except: pass
			
			e=EMData(raw["nx"], raw["ny"], raw["nz"])
			a=EMData(bxsz,bxsz,bxsz)
			a.to_one()
			a.process_inplace("mask.soft",{"outer_radius":-5})
			e.to_zero()
			
			for ii,p in enumerate(pks):
				if ii%100==0: print ii
				e.insert_scaled_sum(a, p.tolist())
			
			e.process_inplace("threshold.clampminmax",{"maxval":1, "minval":0})
			e.write_image(pname)
			
	
	else:
		segname=args[0]
		tomoname=args[1]
		
		e=EMData(segname)
		img=e.numpy()
		img[img<options.thresh]=0
		
		#lb, nlb=ndimage.measurements.label(img)
		#pks=np.array(ndimage.maximum_position(img,lb,range(1,nlb)))
		#n=len(pks)
		#print n
		
		e.process_inplace("mask.onlypeaks")
		#print np.sum(img>0)
		pks= np.array(np.where(img>0)).T
		
		pk_new=[[-100,-100,-100]]
		for p in pks:
			
			nb=np.sum(np.sum(np.array(pk_new-p)**2,axis=1)<(options.boxsz/4)**2)
			#print p, nb
			if nb<1:
				pk_new.append(p)
		pks=np.array(pk_new)
		n=len(pks)
		#e.write_image("tmp2.hdf")
		print "{} boxes found..".format(n)
		allbox=[]
		#if options.sort:
			#denmap=EMData(tomoname)
			#dmp=denmap.numpy()
			#den=np.array(ndimage.measurements.maximum(np.abs(dmp),lb,range(1,nlb)))
			#pks=pks[np.argsort(den)]
		for j in range(n):
			box=pks[j]
			if min(box[2],box[1],e["nx"]-box[2],e["ny"]-box[1])<options.edge:
				continue
			#box*=2
			if options.shrink>0:
				box*=int(options.shrink)
			allbox.append([box[2], box[1],"manual", int(box[0])])

		jsname=info_name(tomoname)
		
		js=js_open_dict(jsname)
		js["boxes"]=allbox
		js=None

	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
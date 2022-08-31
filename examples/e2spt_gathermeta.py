#!/usr/bin/env python
# Muyuan Chen 2022-08
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="", default=None)
	parser.add_argument("--lstin", type=str,help="", default=None)
	parser.add_argument("--path", type=str,help="", default=None)
	#parser.add_argument("--sym", type=str,help="", default="c1")
	
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	lst0=load_lst_params(options.lstin)
	key0=["{}_{:04d}".format(l["src"], l["idx"]) for l in lst0]
	dic0={k:i for i,k in enumerate(key0)}
	lst1=load_lst_params(options.ptcls)
	
	lout3d=[]
	lout2d=[]
	info3d=[]
	info2d=[]
	ndxf=0
	for li,l in enumerate(lst1):
		e=EMData(l["src"], l["idx"], True)
		ky="{}_{:04d}".format(e["orig_ptcl"], e["orig_idx"])
		ii=dic0[ky]
		l0=lst0[ii]
		xf3d=e["xform.align3d"]

		dc3={"src":l["src"], "idx":l["idx"], 
			"score":l0["score"], "xform.align3d":xf3d,
			"orig_idx":ii}
		lout3d.append(dc3)

		imgsrc=e["class_ptcl_src"]
		imgidx=e["class_ptcl_idxs"]
		coord=e["ptcl_source_coord"]
		boxsz=e["nx"]
		
		rhdrs=EMData.read_images(imgsrc,imgidx,IMAGE_UNKNOWN,True)
		idx2d=[]
		for k,i in enumerate(imgidx): 
			h=rhdrs[k]
			xf=Transform(h["xform.projection"])
			xf=xf*(xf3d.inverse())
			dc2={"src":imgsrc,"idx":i, "class":ii%2,
				"ptcl3d_id":li, "xform.projection":xf,
				"tilt_id":h["tilt_id"]}
			if h.has_attr("orig_dxf"):
				dc2["dxf"]=h["orig_dxf"]
				ndxf+=1
				
			dc={"src":imgsrc,"idx":i,
				"idx3d":li,
				"xform.projection":Transform(h["xform.projection"]), "tilt_id":h["tilt_id"]}
			
			idx2d.append(len(lout2d))
			lout2d.append(dc2)
			info2d.append(dc)
				
		dc={"src":l["src"], "idx":l["idx"],
			"coord":coord, "idx2d":idx2d, "xform.align3d":xf3d}
		info3d.append(dc)
		
		sys.stdout.write("\r {}/{}".format(li+1, len(lst1)))
		sys.stdout.flush()
		
	print("Total {} 3d particles and {} 2d particles".format(len(lout3d), len(lout2d)))
	print("{} particles ( {:.1f}% ) has 2d transform".format(ndxf, 100*float(ndxf)/len(lout2d)))

	save_lst_params(lout3d, "{}/aliptcls3d_01.lst".format(options.path))
	save_lst_params(lout2d, "{}/aliptcls2d_01.lst".format(options.path))
	save_lst_params(info2d, "{}/particle_info_2d.lst".format(options.path))
	save_lst_params(info3d, "{}/particle_info_3d.lst".format(options.path))
	
	#for eo in ["even", "odd"]:
		#run(f"e2spa_make3d.py --input {options.path}/aliptcls2d_00.lst --output {options.path}/threed_00_{eo}.hdf --clsid {eo} --outsize {boxsz} --sym {options.sym} --parallel thread:32")
		
	
	E2end(logid)

	
if __name__ == '__main__':
	main()
	

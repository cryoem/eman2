#!/usr/bin/env python
# Muyuan Chen 2022-06
from EMAN2 import *
import numpy as np

def main():
	
	usage="map a lst file of aligned particles to box coordinates in tomograms. "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--label", type=str,help="particle label", default="fromlst")
	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	
	n=EMUtil.get_image_count(args[0])
	hdrs=[EMData(args[0], i, True) for i in range(n)]
	bxsz=hdrs[0]["nx"]
	apix=hdrs[0]["apix_x"]

	coord=[p["ptcl_source_coord"] for p in hdrs]
	iname=[info_name(p["class_ptcl_src"]) for p in hdrs]

	alipm=load_lst_params(args[0])
	alipm=np.array([a["xform.align3d"].inverse().get_trans() for a in alipm])
	coord=np.array(coord)+alipm
	coord=coord.tolist()
	dic={b:[] for b in np.unique(iname)}
	
	for i in range(n):
		b=coord[i]+["fromlst", 0, 0]
		dic[iname[i]].append(b)
		
	for ky in dic.keys():
		boxes=dic[ky]
		
		js=dict(js_open_dict(ky)).copy()
		tfile=js["tlt_file"]
		e=EMData(tfile, 0, True)
		apix_unbin=e["apix_x"]
		scale=apix/apix_unbin
		for b in boxes:
			for i in range(3): b[i]=round(b[i]*scale)
		
		if "class_list" in js:
			mx=list(js["class_list"].keys())
			mx=1+np.max([int(x) for x in mx])
			for d in boxes:
				d[-1]+=mx
		else:
			js["class_list"]={}
			mx=0
		
		js["class_list"][str(mx)]={'boxsize': good_size_small(bxsz*scale), 'name': options.label}
		if "boxes_3d" not in js: js["boxes_3d"]=[]
		js["boxes_3d"].extend(boxes)
		print("{} -> {} particles".format(base_name(ky), len(boxes)))
		f=js_open_dict(ky)
		f.update(js)
		f.close()
    
	
	
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	

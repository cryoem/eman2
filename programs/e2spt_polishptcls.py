#!/usr/bin/env python
# Muyuan Chen 2022-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ali2d", type=str,help="", default="")
	parser.add_argument("--ali3d", type=str,help="", default="")
	parser.add_argument("--outpath", type=str,help="", default="ptcl3d")
	parser.add_argument("--outsize", type=int,help="", default=-1)
	parser.add_argument("--threads", type=int,help="", default=12)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#ali2d=load_lst_params(options.ali2d)
	#a2d=[a for a in ali2d if a["tilt_id"]==30]
	#xfs=[a["xform.projection"] for a in a2d]
	#n=len(xfs)
	#print(len(ali2d), n)

	ali3d=load_lst_params(options.ali3d)
	xf3=[a["xform.align3d"] for a in ali3d]
	print(len(xf3))
	
	try: os.mkdir(options.outpath)
	except: pass
	
	e2=EMData(options.ali2d, 0, True)
	pad=e2["nx"]
	e3=EMData(options.ali3d, 0, True)
	if options.outsize<=0:
		sz=e3["nx"]
	else:
		sz=options.outsize
	
	
	for i in range(len(ali3d)):
		om=f"{options.outpath}/ptcl_{i:03d}_raw.hdf"
		cmd=f"e2spa_make3d.py --input {options.ali2d} --output {om} --p3did {i} --keep 1,.95,1 --outsize {sz} --pad {pad} --sym c1 --threads {options.threads}"
		print(cmd)
		launch_childprocess(cmd)
		hdr=EMData(options.ali3d, i, True)
		nm=base_name(ali3d[i]['src'])
		ni=ali3d[i]['idx']
		fname=f"{options.outpath}/{nm}_{ni:03d}__ptcl.hdf"
		
		e=EMData(om)
		x=Transform(xf3[i]).inverse()
		print(fname, x)
		
		x.set_trans(0,0,0)
		e.process_inplace("xform",{"transform":x})
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./25})
		e.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
		e.mult(-1)
		e.process_inplace("normalize.edgemean")
		e.write_compressed(fname,0, 12, nooutliers=True)
		os.remove(om)

	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	

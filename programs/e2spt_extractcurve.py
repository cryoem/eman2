#!/usr/bin/env python
# Muyuan Chen 2023-02
from EMAN2 import *
from EMAN2_utils import *
import numpy as np


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--topn", type=int,help="save top N longest curves", default=-1)
	parser.add_argument("--curveid", type=int,help="curve id", default=0)
	parser.add_argument("--boxsize", type=int,help="box size. default 64", default=64)
	parser.add_argument("--writetmp", action="store_true", default=False ,help="write tmp files")

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	try:os.mkdir("curves")
	except: pass

	for fname in args:
		print("extract curves from",fname)
		
		info=dict(js_open_dict(info_name(fname))).copy()
		if "curves" not in info:
			print("no curve in ",fname)
			continue
		
		curves=np.array(info["curves"])
		cvs=curves[curves[:,-1]==options.curveid,:-1].copy()
		cvs=[cvs[cvs[:,-1]==i,:-1] for i in np.unique(cvs[:,-1])]
		cvl=[np.sum(np.linalg.norm(np.diff(c, axis=0), axis=1)) for c in cvs]
		print("{} curves, average length {:.1f} pixels".format(len(cvs), np.mean(cvl)))
		
		tomo=EMData(fname)
		sz=np.array([tomo["nx"],tomo["ny"],tomo["nz"]])
		
		ci=np.argsort(cvl)[::-1]
		if options.topn>0:
			ci=ci[:options.topn]
		
		for ii in ci:
			cv=cvs[ii].copy()
			cv=np.unique(cv, axis=0)
			bx=options.boxsize
			cv=cv/4+sz/2
			ln=np.sum(np.linalg.norm(np.diff(cv, axis=0), axis=1))
			cv=interp_points(cv, npt=int(np.round(ln/4)))
			df=cv[1:]-cv[:-1]
			cv=(cv[1:]+cv[:-1])/2.

			ln=np.sum(np.linalg.norm(np.diff(cv, axis=0), axis=1))
			ln=good_size(ln)
			fib=EMData(bx, ln, bx)
			fib.set_attr_dict(tomo.get_attr_dict())

			fib.to_zero()
			wt=fib.copy()
			msk=EMData(bx,bx,bx)
			msk.to_one()
			msk.process_inplace("mask.zeroedge3d",{"y0":20,"y1":20})
			msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.025})
			
			if options.writetmp:
				pname="curves/ptcls_tmp.hdf"
				if os.path.isfile(pname):
					os.remove(pname)

			pos=ln#bx/2
			for i,c in enumerate(cv):
				xf=Transform()

				d=df[i].copy()
				dl=np.linalg.norm(d)
				d/=dl
				xf.set_rotation(d.tolist())
				xf=Transform({"type":"xyz","ztilt":90})*xf
				xf=xf*Transform({"type":"xyz","ytilt":90})

				xf.set_trans(c.tolist())
				e=tomo.get_rotated_clip(xf, (bx,bx,bx))
				e.mult(msk)

				if options.writetmp:
					e.write_image(pname, -1)
					
				fib.insert_scaled_sum(e, (bx/2,pos,bx/2))
				wt.insert_scaled_sum(msk, (bx/2,pos,bx/2))
				pos-=dl

			wt.process_inplace("math.reciprocal")
			fib.mult(wt)
			# fib.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./30})
			fib.process_inplace("filter.highpass.gauss",{"cutoff_freq":1./400})
			fib.write_image(f"curves/{base_name(fname)}_{ii:02d}.hdf")
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	


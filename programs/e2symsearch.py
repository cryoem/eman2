#!/usr/bin/env python
# Muyuan Chen 2022-05
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize
#import queue
#import threading

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sym", type=str,help="symetry", default="c2")
	parser.add_argument("--boxsz", type=int,help="", default=64)
	#parser.add_argument("--ntry", type=int,help="number of tries", default=20)
	parser.add_argument("--applysym", action="store_true", default=False ,help="apply symmetry after alignment")
	parser.add_argument("--refineonly", action="store_true", default=False ,help="start near correct solution")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	fname=args[0]
	sz=options.boxsz
	e=EMData(fname)
	sym=options.sym
	
	def test_rot_refine(x):
		
		xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2],"tx":x[3],"ty":x[4],"tz":x[5]})
		ref0=ref.process("xform", {"transform":xf})
		ref_sym=ref0.process("xform.applysym",{"sym":sym})
		fsc=ref0.calc_fourier_shell_correlation(ref_sym)
		fsc=np.array(fsc).reshape((3,-1))[1]
			
		return -np.mean(fsc)
	
	if options.refineonly:
		x0=[0,0,0,0,0,0]
		ref=e.copy()
		res=minimize(test_rot_refine, x0, method='Powell', 
					options={'ftol': 1e-3, 'disp': False, "maxiter":20})
		
		x=res.x
		xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2],"tx":x[3],"ty":x[4],"tz":x[5]})
		a=ref.process("xform", {"transform":xf})
		a["xform.align3d"]=xf
		print(xf)
		if options.applysym:
			a.process_inplace("xform.applysym", {"averager":"mean.tomo", "sym":sym})
		outname=fname[:-4]+"_sym.hdf"
		a.write_image(outname)
		print("Done. Output written to {}".format(outname))
		
		E2end(logid)
		return
	
	scale=e["nx"]/sz
	e.process_inplace("math.fft.resample",{"n":scale})
	e.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.45})
	c=e.calc_center_of_mass(0)
	c=np.array(c)-sz/2
	e.process_inplace("xform", {"tx":-c[0], "ty":-c[1], "tz":-c[2]})
	print(c)
	#e.process_inplace("xform.centerofmass")
	sym=options.sym
	
	ref=e.copy()
	ref=ref.do_fft()
	ref.process_inplace("xform.phaseorigin.tocenter")
	ref.process_inplace("xform.fourierorigin.tocenter")
	
	xf=Transform()
	astep=7.4
	symx=Symmetries.get('c1')
	xfcrs=symx.gen_orientations("eman",{"delta":astep,"inc_mirror":1})
	score=[]
	nsym=xf.get_nsym(sym)
	t0=Transform()
	for xf in xfcrs:
		
		pj0=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
		scr=[]
		for i in range(1,nsym):
			xs=t0.get_sym(sym,i)*xf
			pj1=ref.project('gauss_fft',{"transform":xs, "returnfft":1})
			fsc=pj0.calc_fourier_shell_correlation(pj1)
			fsc=np.array(fsc).reshape((3,-1))[1, 1:20]
			scr.append(np.mean(fsc))
			
		score.append(np.mean(scr))
	
	score=np.array(score)
	idx=np.argsort(-score)[:8]
	
	
	def test_rot(x):
		
		xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2]})
		xt=Transform({"type":"eman","tx":x[3],"ty":x[4],"tz":x[5]})
		# xt=Transform()

		pj0=ref.project('gauss_fft',{"transform":xf*xt, "returnfft":1})
		scr=[]
		for i in range(1,nsym):
			xs=t0.get_sym(sym,i)*xf
			pj1=ref.project('gauss_fft',{"transform":xs*xt, "returnfft":1})
			fsc=pj0.calc_fourier_shell_correlation(pj1)
			fsc=np.array(fsc).reshape((3,-1))[1, 1:64]
			scr.append(np.mean(fsc))
			
		return -np.mean(scr)
	
	
	xfnew=[xfcrs[i] for i in idx]
	score=[]
	xfrefine=[]
	for xx in xfnew:
		x=xx.get_params('xyz')
		x0=[x["xtilt"], x["ytilt"], x["ztilt"],0,0,0]
		
		res=minimize(test_rot, x0, method='Powell', 
					options={'ftol': 1e-3, 'disp': False, "maxiter":20})
		score.append(res.fun)
		xfrefine.append(res.x)
		print(test_rot(x0), test_rot(res.x))
		
	
	a=EMData(fname)
	c=c*a["nx"]/sz
	a.process_inplace("xform", {"tx":-c[0], "ty":-c[1], "tz":-c[2]})

	ref=a.do_fft()
	ref.process_inplace("xform.phaseorigin.tocenter")
	ref.process_inplace("xform.fourierorigin.tocenter")

	x0=xfrefine[np.argmin(score)].copy()
	x0[3:]*=scale
	xf=Transform({"type":"xyz", "xtilt":x0[0], "ytilt":x0[1], "ztilt":x0[2]})
	xf.set_trans(x0[3],x0[4],x0[5])
	
	res=minimize(test_rot, x0, method='Powell', 
				options={'ftol': 1e-3, 'disp': False, "maxiter":40})
	score.append(res.fun)
	x=res.x
	xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2]})
	xf.set_trans(x[3],x[4],x[5])
	print(test_rot(x0), test_rot(x))
		
	a=EMData(fname)
	t=Transform({"type":"eman", "tx":-c[0], "ty":-c[1], "tz":-c[2]})
	xf=xf*t
	a.process_inplace("xform",{"transform":xf})
	print(xf)
	a["xform.align3d"]=xf
	#outname=fname[:-4]+"_sym.hdf"
	#a.write_image(outname)
	if options.applysym:
		a.process_inplace("xform.applysym", {"averager":"mean.tomo", "sym":sym})
	outname=fname[:-4]+"_sym.hdf"
	a.write_image(outname)
	print("Done. Output written to {}".format(outname))
	
	E2end(logid)
	
if __name__ == '__main__':
	main()
	

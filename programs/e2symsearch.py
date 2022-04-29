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
	#parser.add_argument("--threads", type=int,help="threads", default=12)
	#parser.add_argument("--ntry", type=int,help="number of tries", default=20)
	parser.add_argument("--applysym", action="store_true", default=False ,help="apply symmetry after alignment")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	fname=args[0]
	e=EMData(fname)
	scale=e["nx"]/48
	e.process_inplace("math.fft.resample",{"n":scale})
	e.process_inplace("xform.centerofmass")
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
			fsc=np.array(fsc).reshape((3,-1))[1, 1:-4]
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
	a.process_inplace("xform.centerofmass")
	ref=a.do_fft()
	ref.process_inplace("xform.phaseorigin.tocenter")
	ref.process_inplace("xform.fourierorigin.tocenter")

	x0=xfrefine[np.argmin(score)].copy()
	x0[3:]*=scale

	res=minimize(test_rot, x0, method='Powell', 
				options={'ftol': 1e-3, 'disp': False, "maxiter":20})
	score.append(res.fun)
	x=res.x
	xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2]})
	xf.set_trans(x[3],x[4],x[5])
	print(test_rot(x0), test_rot(x))
	print(xf)
		
	a.process_inplace("xform",{"transform":xf})
	a["xform.align3d"]=xf
	if options.applysym:
		a.process_inplace("xform.applysym", {"averager":"mean.tomo", "sym":sym})
	outname=fname[:-4]+"_sym.hdf"
	a.write_image(outname)
	print("Done. Output written to {}".format(outname))
	
	E2end(logid)
	
if __name__ == '__main__':
	main()
	

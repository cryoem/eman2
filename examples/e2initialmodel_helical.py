#!/usr/bin/env python
# Muyuan Chen 2021-04
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize

def main():
	#####################
	def test_rot(x, returnxf=False):
		if isinstance(x, Transform):
			xf=x
		else:
			x=list(x)
			x.extend(curxf[len(x):])
			xf=Transform({"type":"xyz", "xtilt":x[0], "ytilt":x[1], "ztilt":x[2],"tx":x[3], "ty":x[4]})

		pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
		#### only use the first 50 fourier pixels
		x0=1; x1=50

		ccf=e.calc_ccf(pj)
		mxsft=40
		pos=ccf.calc_max_location_wrap(mxsft, mxsft, 0)
		xf.translate(pos[0], pos[1],0)
		pj.process_inplace("xform", {"tx":pos[0], "ty":pos[1]})

		fsc=e.calc_fourier_shell_correlation(pj)
		fsc=np.array(fsc).reshape((3,-1))[:, x0:x1]

		wt=np.ones_like(fsc[1])

		scr=-np.sum(fsc[1]*wt)/np.sum(wt)

		if returnxf:
			return scr, xf
		else:
			return scr
	#####################
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="ptcls", default=None)
	parser.add_argument("--ref", type=str,help="ref", default=None)
	parser.add_argument("--output", type=str,help="output", default="output.hdf")
	parser.add_argument("--niter", type=int,help="number of iters", default=3)
	parser.add_argument("--sym", type=str,help="sym: number of copy:1:rotation:rise in pixel. (h15:1:-1:4.53)", default="h15:1:-1:4.53")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	ref=EMData(options.ref)
	sz=ref["ny"]
	pad=good_size(sz*1.4)
	for itr in range(options.niter):
		ref=ref.do_fft()
		ref.process_inplace("xform.phaseorigin.tocenter")
		ref.process_inplace("xform.fourierorigin.tocenter")
		
		n=EMUtil.get_image_count(options.ptcls)
		idx=range(n)
		output=[]
		for ii in idx:
			e=EMData(options.ptcls,ii)
			e=e.do_fft()
			e.process_inplace("xform.phaseorigin.tocenter")
			rg=np.arange(0,180,5)
			score=[]
			poss=[]
			for az in rg:
				xf=Transform({"type":"eman", "alt":90, "az":float(az)})
				pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
				ccf=e.calc_ccf(pj)
				pos=ccf.calc_max_location_wrap(50,50,0)
				pj.process_inplace("xform",{"tx":pos[0], "ty":pos[1]})

				fsc=pj.calc_fourier_shell_correlation(e)
				fsc=np.array(fsc).reshape((3,-1))

				score.append(np.mean(fsc[1,:40]))
				poss.append(pos)

			ix=np.argmax(score)
			output.append([ii, rg[ix], poss[ix]])
			print(output[-1], np.max(score))
				
			
		output2=[]
		for o in output:
			e=EMData(options.ptcls, o[0])
			xf0=Transform({"type":"eman", "alt":90, "az":float(o[1])})
			e=e.do_fft()
			e.process_inplace("xform.phaseorigin.tocenter")
			
			x=xf0.get_params("xyz")
			curxf=[x["xtilt"], x["ytilt"], x["ztilt"],x["tx"],x["ty"]]
			x0=curxf[:3]
			res=minimize(test_rot, x0, method='Powell', options={'ftol': 1e-3, 'disp': False, "maxiter":20})
			scr, x=test_rot(res.x, True)
			print(o,scr, x)
			output2.append([o[0], x])
			
		t=Transform()
		recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"trilinear"})
		recon.setup()
		for o in output2:
			e=EMData(options.ptcls, o[0])
			xf=o[1]
			nsym=t.get_nsym(options.sym)
			for i in range(nsym):
				x=t.get_sym(options.sym, i)
				mp=recon.preprocess_slice(e, x*xf)
				recon.insert_slice(mp,xf,1)
			
		threed=recon.finish(True)
		threed.clip_inplace(Region((pad-sz)//2, (pad-sz)//2, (pad-sz)//2, sz, sz, sz))
		threed.process_inplace("normalize.edgemean")
		outname=options.output[:-4]+"_{:02d}.hdf".format(itr)
		threed.write_image(outname)
		
		s=[float(i) for i in options.sym[1:].split(":")]
		sym2="h5:1:{:.2f}:{:.2f}".format(s[2]*10,s[3]*10)
		print(sym2)
		tsym=threed.process("xform.applysym",{'sym':sym2})
		tsym.write_image(outname[:-4]+"_sym.hdf")
		tsym.process_inplace("normalize.edgemean")
		ref=tsym.copy()

	E2end(logid)
	
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
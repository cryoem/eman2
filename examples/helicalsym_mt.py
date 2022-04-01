#!/usr/bin/env python
# Muyuan Chen 2018-02

from EMAN2 import *
import numpy as np

def main():
	
	usage="apply microtubule helical symmetry "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sym", type=str,help="sym", default="")
	parser.add_argument("--p0", type=float,help="starting angle of seam line", default=-138)
	parser.add_argument("--w", type=float,help="width  of seam line", default=10)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	sym=options.sym
	e=EMData(args[0])
	if sym=="":
		rise=40.7*3/13/e["apix_x"]
		sym="h13:1:27.6923:{:.4f}".format(rise)
		print("using symmetry", sym)
		
	e.process_inplace("filter.highpass.gauss", {"cutoff_pixels":4})
	
	xf=Transform()
	xfs=[xf.get_sym(sym, i) for i in range(26)]
	p0=options.p0
	avg=e.copy()
	msk=e.copy()
	msk.to_one()
	nrm=msk.copy()

	for x in xfs[1:]:
		a=e.process("xform", {"transform":x})
		phi=x.get_params("eman")["phi"]
		tz=x.get_params("eman")["tz"]
	#	 print(tz, phi)
		if tz>0:
			phi=360-phi
			pcen=p0+phi/2
		else:
			pcen=p0-phi/2
			
		m=msk.process("mask.cylinder",{"outer_radius":e["nx"]//2, "phicen": pcen, "phirange":phi/2-options.w})
		m.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.02})
		a.mult(m)
		avg.add(a)
		nrm.add(m)
		a.write_image("threed_xfs.hdf", -1)
		
	nrm.process_inplace("math.reciprocal", {"zero_to":1})
	avg.mult(nrm)
	
	
	if len(args)>1:
		out=args[1]
	else:
		out=args[0][:-4]+"_sym.hdf"
	print("Writing to {}".format(out))
	if os.path.isfile(out):
		os.remove(out)
	avg.write_image(out)


	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
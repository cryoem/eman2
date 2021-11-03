#!/usr/bin/env python
# Muyuan Chen 2018-05
from EMAN2 import *
import numpy as np
import queue
import threading

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sym", type=str,help="symetry", default="c2")
	#parser.add_argument("--threads", type=int,help="threads", default=12)
	parser.add_argument("--ntry", type=int,help="number of tries", default=20)
	parser.add_argument("--applysym", action="store_true", default=False ,help="apply symmetry after alignment")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	fname=args[0]
	e=EMData(fname)
	e.process_inplace("math.fft.resample",{"n":4})
	sym=options.sym
	ntry=options.ntry
	
	s=parsesym(sym)
	oris=s.gen_orientations("rand",{"n":ntry, "phitoo":True})
	oris[0]=Transform()
	oris=[Transform({"type":"eman", "phi":i*18}) for i in range(20)]
	jsd=queue.Queue(0)
	thrds=[threading.Thread(target=do_ali,args=([e, o, sym, jsd])) for o in oris]
	for t in thrds: t.start()
	for t in thrds: t.join()
	
	alis=[]
	while not jsd.empty():
		alis.append(jsd.get())
	#print alis
	scr=[a["score"] for a in alis]
	im=np.argmin(scr)
	a=alis[im].copy()
	xf=a["xform.align3d"]
	
	e=EMData(fname)
	xf.set_trans(xf.get_trans()*4)
	a=e.align("symalignquat", e, {"sym":sym, "xform.align3d":xf}, "ccc")
	#a=e.process("xform",{"transform":xf})
	
	#a.process_inplace("xform.centeracf")
	if sym.startswith('c'):
		ac=a.copy()
		ac.process_inplace("xform.applysym", {"sym":sym})
		ac.process_inplace("mask.cylinder",{"outer_radius":a["nx"]*3/8, "phirange":360})
		ac.process_inplace("xform.centerofmass")
		ct=ac["xform.align3d"].get_trans()
		xf.translate(0,0,ct[2])
		a=e.copy()
		a.transform(xf)
	
	print(xf)
	a["xform.align3d"]=xf
	if options.applysym:
		a.process_inplace("xform.applysym", {"averager":"mean.tomo", "sym":sym})
	outname=fname[:-4]+"_sym.hdf"
	a.write_image(outname)
	print("Done. Output written to {}".format(outname))
	
	E2end(logid)
	
def do_ali(e,o,sym,jsd):
	a=e.align("symalignquat", e, {"sym":sym, "xform.align3d":o}, "ccc")
	jsd.put(a)
	
if __name__ == '__main__':
	main()
	
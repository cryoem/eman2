#!/usr/bin/env python
# Muyuan Chen 2017-10
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from EMAN2 import *
import numpy as np

def main():
	
	usage="Test EMAN2 functionalities.. "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="output file", default="testout")
	parser.add_argument("--binaryout", action="store_true",help="binary output", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	t0=time.time()
	e=test_image()
	tosave=[]
	
	#### basic processors

	e.process_inplace("math.fft.resample",{"n":2})
	e.process_inplace("normalize")
	e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.22})
	erot=e.copy()
	erot.rotate(0,0,30)
	erot.process_inplace("math.addnoise",{"noise":3, "seed":123})
	erot.process_inplace("normalize")

	#### aligner

	eali=erot.align("rotate_translate_flip", e, {}, "ccc")
	tosave.append(eali.numpy().flatten().copy())
	
	#### reconstruction

	sz=e["nx"]
	pad=good_boxsize(int(sz*1.5))
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"gauss_2"})
	recon.setup()
	epad=eali.get_clip(Region(old_div(sz,2)-old_div(pad,2), old_div(sz,2)-old_div(pad,2), pad, pad))
	xfs=[{}, {"type":"eman", "alt":90}, {"type":"eman", "alt":90, "az":90}]
	for xf in xfs:
		recon.preprocess_slice(epad, Transform(xf))
		recon.insert_slice(epad, Transform(xf), 1)
	threed=recon.finish(True)
	threed.clip_inplace(Region(old_div((pad-sz),2), old_div((pad-sz),2), old_div((pad-sz),2), sz,sz,sz))
	tdnpy=threed.numpy().copy()
	tosave.append(tdnpy.flatten())
	
	#### projection
	sym=parsesym("d2")
	xfs=sym.gen_orientations("eman",{"n":10})
	epj=threed.project("standard", xfs[5])
	pjnp=epj.numpy().copy()
	tosave.append(pjnp.flatten())
	
	tosave=np.hstack(tosave)
	
	if options.binaryout:
		np.save(options.output, tosave)
	else:
		np.savetxt(options.output, tosave)


	
	print("Done. Time elapse: {:.3f}s".format(float(time.time()-t0)))
	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	

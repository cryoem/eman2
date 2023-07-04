#!/usr/bin/env python
# Muyuan Chen 2023-02
from EMAN2 import *
import numpy as np

def main():
	
	usage="""
	Merge results from multiple focused refinement runs with different masks. For example, when you have a global refinement in gmm_00, and focused refinement in gmm_01 and gmm_02 using differnet masks, run
	e2gmm_merge_patch.py gmm_01 gmm_02 --base gmm_00
	
	The result will be labeled iteraton 99 in the base folder.
	This is also called automatically by the patch-by-patch refinement.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--base", type=str, help="path of the global refinement. required for merging multiple folders",default=None)
	parser.add_argument("--sym", type=str, help="",default="c1")
	parser.add_argument("--skippp", action="store_true", default=False ,help="skip post process")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	if len(args)==1:
		## read from gmm_refine_patch in a single folder
		path=args[0]
		print("Reading from patch-wise refinement in",path)
		iprg=[]
		for ip in range(100):
			mfile=f"{path}/mask_patch_{ip:02d}.hdf"
			if not os.path.isfile(mfile):
				break
			
			iprg.append(ip)
			
			
		for eo in ["even","odd"]:
			avg=EMData(f"{path}/threed_00_{eo}.hdf")
			avg.to_zero()
			wt=avg.copy()

			for ci in iprg:
				m=EMData(f"{path}/mask_patch_{ci:02d}.hdf")
				e=EMData(f"{path}/threed_patch_{ci:02d}_raw_{eo}.hdf")
				e.mult(m)
				avg.add(e)
				wt.add(m)
				
			wt.process_inplace("math.reciprocal", {"zero_to":0})
			avg.mult(wt)
			avg.write_image(f"{path}/threed_99_{eo}.hdf")
		
		
		if not options.skippp:
			run(f"e2refine_postprocess.py --even {path}/threed_99_even.hdf --res 5 --tophat localwiener --thread 32 --setsf sf.txt --align")
		
	elif len(args)>1:
		print("Merging refinement results from multiple folders..")
		if not options.base:
			print("A global refinement folder (--base) required")
			return
		
		for eo in ["even","odd"]:
			avg=EMData(f"{options.base}/threed_raw_{eo}.hdf")
			avg.to_zero()
			wt=avg.copy()
			
			m0=EMData(f"{options.base}/mask.hdf")
			for pt in args:
				if pt==options.base: continue
				js=js_open_dict(f"{pt}/0_gmm_params.json")
				mfile=js["mask"]
				m=EMData(mfile)
				print(f"Reading from {pt}, using mask file {mfile}")
				m.process_inplace("threshold.binary",{"value":.5})
				m.process_inplace("mask.addshells.gauss",{"val1":2,"val2":4})
				e=EMData(f"{pt}/threed_raw_{eo}.hdf")
				e.mult(m)
				if options.sym!="c1":
					nsym=Transform.get_nsym(options.sym)
					e.process_inplace("xform.applysym",{"sym":options.sym})
					m.process_inplace("xform.applysym",{"sym":options.sym})
					e.mult(nsym)
					m.mult(nsym)
				avg.add(e)
				wt.add(m)
				m0.add(m*-1)

			e=EMData(f"{options.base}/threed_raw_{eo}.hdf")
			# m=EMData(f"mask_zone0.hdf")
			m=m0.copy()
			m.process_inplace("threshold.binary",{"value":.2})
			m.process_inplace("mask.addshells.gauss",{"val1":5,"val2":10})
			e.mult(m)
			
			avg.add(e)
			wt.add(m)
			# m=EMData("gmm_00/mask.hdf")
			
			wt.process_inplace("threshold.clampminmax", {"minval":1e-3, "maxval":20})
			wt.process_inplace("math.reciprocal", {"zero_to":0})
			avg.mult(wt)
			avg.write_image(f"{options.base}/threed_99_{eo}.hdf")
			
		if not options.skippp:
			run(f"e2refine_postprocess.py --even {options.base}/threed_99_even.hdf --res 5 --tophat localwiener --thread 32 --setsf sf.txt --align")
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	



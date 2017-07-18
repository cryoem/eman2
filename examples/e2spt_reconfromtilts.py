#!/usr/bin/env python
# Muyuan Chen 2017-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--json", type=str,help="e2spt_align.py output json file. Usually particle_parms_*.json", default=None)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--m3dkeep", type=float,help="The fraction of slices to keep, based on quality scores", default=1)
	parser.add_argument("--threads", type=int, help="Number of threads", default=10)
	parser.add_argument("--first", type=int, help="Use only first N particles to test", default=-1)
	parser.add_argument("--mass", type=int, help="Mass", default=1000)
	parser.add_argument("--clip", type=int, help="Clip the volume at the end", default=-1)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	sptpath=os.path.dirname(options.json)
	js=js_open_dict(options.json)
	print "Read alignment of {:d} 3D particles.".format(len(js.keys()))
	ks=[eval(str(s)) for s in sorted(js.keys())]
	
	### lists in dictionaries in dictionaries.....
	ptclnames={}
	for keys in ks:
		e=EMData(keys[0],keys[1], True)
		try:
			srcfile = e["source_2dptcl"]
		except:
			print "No source 2D particles found for the 3D particles. The program will not work. Please make sure the 3D particles are generated from 2D particles through sptmake3d.py."

		mid = e["data_n"]
		if srcfile not in ptclnames:
			ptclnames[srcfile]={}
			
		ptclnames[srcfile][mid]=[keys[1], js[keys]['xform.align3d']]
		
	
	mid=-1
	lname=[os.path.join(sptpath, "ali_ptcls_{}.lst".format(eo)) for eo in ["even", "odd"]]
	for l in lname:
		try: os.remove(l)
		except:pass
	
	lst=[LSXFile(m, False) for m in lname]
	
	pnames=ptclnames.keys()
	for pname in pnames:
		modelids=sorted(ptclnames[pname].keys())
		num=EMUtil.get_image_count(pname)
		print "Processing {}, including {:d} particles of {:d} models".format(pname, num, len(modelids))
		for i in range(num):
			p=EMData(pname, i, True)
			mi=p["model_id"]
			if mi not in modelids:
				print "Cannot find particle {}:{:d} in the alignment result. The program will proceed but something wrong may be going on. Please check your protocol...".format(pname, i)
				continue
			xf=ptclnames[pname][mi][1]
			xi=ptclnames[pname][mi][0]
			xx=xf.inverse()
			pj=p["xform.projection"]
			pj=pj*xx
			rr=pj.get_params("eman")
			
			lst[xi%2].write(-1, i, pname, str(rr))

	for l in lst: l=None			

	print "Aligned particle lists written to {} and {}.".format(lname[0], lname[1])
	
	#################
	
	print "Now make 3D volume. This may take a while...."
	
	e=EMData(pname, 0, True)
	boxsize=e["nx"]
	pad=good_size(boxsize*3/2)
	apix=e["apix_x"]
	
	oname=[os.path.join(sptpath, "map_from_ali_{}.hdf".format(eo)) for eo in ["even", "odd"]]
	for i in [0,1]:
		run("make3dpar_rawptcls.py --input {inp} --sym {sym} --output {out} --keep {kp:.2f} --pad {pd:d} --mode gauss_2 --threads {thd:d} --apix {ap:f}".format(inp=lname[i], out=oname[i], pd=pad, thd=options.threads, ap=apix, sym=options.sym, kp=options.m3dkeep ))
	
	combfile=os.path.join(sptpath, "map_from_ali_avg.hdf")
	
	if options.clip>0:
		for i in [0,1]:
			tname=oname[i][:-4]+"_clip.hdf"
			run("e2proc3d.py {} {} --clip {:d} --process normalize.edgemean".format(oname[i], tname, options.clip))
			oname[i]=tname
	
	run("e2refine_postprocess.py --even {} --odd {} --output {} --mass {}  --sym {}".format(oname[0], oname[1], combfile, options.mass, options.sym))
	
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
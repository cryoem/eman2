#!/usr/bin/env python
# Muyuan Chen 2015-03
from EMAN2 import *
import numpy as np

def main():
	
	usage="""
	Iterative protocol for single particle tomogram refinement. Runs refinement in fully gold-standard way using projections of sub-tomograms.
	prog --path [name of refinement path] --model [initial model] --ptcls [2d particle stack from refinetomo_buildptcls.py]
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path of refinement", default="refine3d2d_00")
	parser.add_argument("--model", type=str,help="initial model", default=None)
	parser.add_argument("--ptcl", type=str,help="particle file generated from refinetomo_buildptcls.py", default=None)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--niter", type=int,help="number of iterations", default=3)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	tarres=50
	e=EMData(options.ptcl,0,True)
	apix=e["apix_x"]
	sz=e["nx"]
	
	try: os.mkdir(options.path)
	except: pass
	
	#### split the data
	lstfile=options.ptcl
	num=EMUtil.get_image_count(lstfile)
	lsteo=[lstfile[:-4]+"_even.lst", lstfile[:-4]+"_odd.lst"]
	for l in lsteo:
		try:
			os.remove(l)
		except: 
			pass
	lsteo=[LSXFile(l,False) for l in lsteo]
	mid=0
	for i in range(num):
		e=EMData(lstfile, i, True)
		lsteo[e["model_id"]%2].write(-1,  e["data_n"],e["data_source"])

	lsteo=None
	
	for eo in ["even","odd"]:
		ptclname=options.ptcl[:-4]+"_{}.lst".format(eo)
		run("e2proc3d.py {model} {path}/threed_00_{eo}.hdf --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix}".format(model=options.model,path=options.path,freq=1.0/(tarres),apix=apix, eo=eo))
		
		run("e2project3d.py {path}/threed_00_{eo}.hdf --outfile {path}/projections_01_{eo}.hdf -f --projector standard --orientgen eman:delta=8:inc_mirror=0:perturb=0 --sym {sym} --parallel thread:12 ".format(path=options.path, eo=eo, sym=options.sym))
		
		run("e2simmx2stage.py {path}/projections_01_{eo}.hdf {ptcl} {path}/simmx_01_{eo}.hdf {path}/proj_simmx_01_{eo}.hdf {path}/proj_stg1_01_{eo}.hdf {path}/simmx_stg1_01_{eo}.hdf --saveali --cmp frc:maxres=30.0 --align rotate_translate_flip --aligncmp ccc --ralign refine --raligncmp ccc --shrinks1 2 --parallel thread:12".format(path=options.path, eo=eo, ptcl=ptclname))
		
		run("refine_tomo.py --ptcl {ptcl} --mapfile {path}/threed_00_{eo}.hdf --lstout {path}/aliptcl_01_{eo}.lst --simmx {path}/simmx_01_{eo}.hdf --pjfile {path}/projections_01_{eo}.hdf".format(eo=eo, path=options.path, ptcl=ptclname))
		
		run("make3dpar_rawptcls.py --input {path}/aliptcl_01_{eo}.lst --sym {sym} --output {path}/threed_01_{eo}.hdf --keep 0.7 --apix {apix} --pad {pad} --mode gauss_5 --threads 12".format(eo=eo, path=options.path, sym=options.sym,apix=apix, pad=int(sz*1.5)))
	
	for it in range(1,options.niter):
		for eo in ["even", "odd"]:
			ptclname=options.ptcl[:-4]+"_{}.lst".format(eo)
			cmd="refine_tomo.py --ptcl {ptcl} --mapfile {path}/threed_{it0:02d}_{eo}.hdf --lstout {path}/aliptcl_{it1:02d}_{eo}.lst --lstin {path}/aliptcl_{it0:02d}_{eo}.lst".format(eo=eo, path=options.path, it0=it, it1=it+1, ptcl=ptclname)
			run(cmd)
			
			cmd="make3dpar_rawptcls.py --input {path}/aliptcl_{it1:02d}_{eo}.lst --sym {sym} --output {path}/threed_{it1:02d}_{eo}.hdf --keep 0.7 --apix {apix} --pad {pad} --mode gauss_5 --threads 12".format(eo=eo, path=options.path, it1=it+1, sym=options.sym,apix=apix, pad=int(sz*1.5))
			run(cmd)
			
		cmd="e2refine_postprocess.py --even {path}/threed_{it1:02d}_even.hdf --odd {path}/threed_{it1:02d}_odd.hdf --output {path}/threed_{it1:02d}.hdf --automaskexpand -1 --align --mass 1200.0 --iter {it1} --sym={sym} --restarget=30.0 --underfilter --ampcorrect=strucfac".format(path=options.path, it1=it+1, sym=options.sym)
		run(cmd)
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
#!/usr/bin/env python

from EMAN2 import *
import os
import sys
import time
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2compare.py [options] < EM map >"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--pdb", default=None,required=True, type=str, help="Image stack containing phase-flipped particles used for alignment")
	parser.add_argument("--sym", default="c1", type=str, help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos.")
	#parser.add_argument("--flip", action="store_true", default=False,help = "Use the traditional e2make3d program instead of the new e2make3dpar program")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	############ Parse options

	struct = args[0]
	base,ext = os.path.splitext(struct)
	dirpath = base_name(struct)

	try: 
		hdr = EMData(struct,0,True)
	except:
		print("Could not load structure header.")
		exit(1)

	try:
		apix = hdr["apix_x"]
		res=2*apix
	except:
		print("Could not determine the apix of the input structure from the header.")
		exit(1)

	try: 
		box = hdr["nx"]
	except:
		print("Could not determine the box size from the header.")
		exit(1)
	
	pdb = options.pdb
	pdb_base = os.path.basename(pdb)
	mrc = "{}.mrc".format(pdb_base)
	ali = "{}_ali.mrc".format(pdb_base)

	sym = options.sym

	fsc = "{}_vs_{}_fsc.txt".format(dirpath,pdb_base)

	############### convert

	print("Converting {} to {}".format(pdb,mrc))
	print("Filtering to {} angstroms.".format(res))
	cmd = "e2pdb2mrc.py {} {} --res={} --box={} --apix={}".format(pdb,mrc,res,box,apix)
	run(cmd)

	############### flip

	#if options.flip:
	#	print("Flipping the handedness of {}".format(mrc))
	#	flp = "{}_flip.mrc".format(pdb_base)
	#	cmd = "e2proc3d.py {} {} --process xform.flip:axis=z".format(mrc,flp)
	#	run(cmd)
	#	mrc = flp

	############### align 

	print("Aligning {} to {}".format(mrc,struct))
	cmd = "e2proc3d.py {} {} --align=rotate_translate_3d_tree --align=refine_3d --alignref={} --sym {}".format(mrc,ali,struct,sym)
	if "d" in "sym": cmd += " --alignctod"
	run(cmd)

	############### compare 

	print("Computing FSC between {} and {}".format(ali,struct))
	cmd = "e2proc3d.py {} {} --calcfsc {}".format(ali,fsc,struct)
	run(cmd)
	print("FSC saved as: {}".format(fsc))


def run(command):
	print "{}: {}".format(time.ctime(time.time()),command)
	append_html("<p>{}: {}</p>".format(time.ctime(time.time()),command),True)
	ret=launch_childprocess(command)
	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)
	return

if "__name__" == "__main__":
	main()

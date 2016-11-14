#!/usr/bin/env python
# Muyuan Chen 2015-10
from EMAN2 import *
import numpy as np
import random
import os

def main():
	
	usage="""e2pathwalker_auto.py <density map> [options]
	Run the whole pathwalker protocol"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="Output folder", default=None)
	parser.add_argument("--natoms", type=int,help="Number of C-alpha atoms", default=500)
	parser.add_argument("--denthr", type=float,help="Density threshold for segmentation", default=-1)
	parser.add_argument("--mapweight", type=int,help="Weight for the density map in pathwalker.", default=200)
	parser.add_argument("--minhlxlen", type=int,help="Minimum length of an alpha helix.", default=4)
	parser.add_argument("--minshtlen", type=int,help="Minimum length of a beta sheet.", default=3)
	parser.add_argument("--nsht", type=int,help="Maximum number of beta sheet strains.", default=30)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	## make folder 
	
	if options.output == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:9]=="pathwalk_" and len(i)==11 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.output = "pathwalk_{:02d}".format(max(fls)+1)
	try: os.makedirs(options.output)
	except: pass
	
	
	mapin=os.path.join(options.output,"map.mrc")
	path=os.path.join(options.output,"path.pdb")
	segout=os.path.join(options.output,"pseudoatoms.pdb")
	hlx=os.path.join(options.output,"helix.pdb")
	sht=os.path.join(options.output,"sheet.pdb")
	
	cmd="e2proc3d.py {} {} --process threshold.belowtozero".format(args[0], mapin)
	run(cmd)
	
	if options.denthr<0:
		e=EMData(mapin)
		print "Density threshold not provided, using mean+sigma..."
		options.denthr=e["mean"]+e["sigma"]
	
	cmd="e2segment3d.py {} --pdbout={} --process=segment.kmeans:ampweight=1:nseg={}:verbose=1:minsegsep=1:pseudoatom=1:thr={}".format(mapin,segout,options.natoms,options.denthr)
	run(cmd)
	
	cmd="e2pathwalker.py {} --mapfile={} --output={} --solver=lkh --overwrite --dmin=1 --dmax=10 --mapthresh={} --mapweight={}".format(segout, mapin, path, options.denthr, options.mapweight)
	run(cmd)
	
	cmd="e2pwhelixfit.py --mapin {} --pdbin {} --output {} --denthr {}  --mapwohelix map_nohlx.mrc --minlen {} --lenthr 10".format(mapin, path, hlx, options.denthr, options.minhlxlen)
	run(cmd)
	
	cmd="e2pwsheetfit.py --pdbin {} --output {} --nsht {} --minlen {}".format(hlx, sht, options.minshtlen, options.nsht)
	run(cmd)
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
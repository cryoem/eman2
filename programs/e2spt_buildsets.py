#!/usr/bin/env python

# Muyuan Chen 2018-07
from __future__ import print_function
from __future__ import division
from EMAN2 import *
import numpy as np


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particle_stacks",help="Specify particles input.", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="sets")
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=3, col=0, rowspan=1, colspan=1, mode="sets")
	parser.add_argument("--label", type=str,help="label of particles for sets", default="", guitype='strbox',row=2, col=0,rowspan=1, colspan=1, mode="sets")
	parser.add_argument("--allparticles", action="store_true", default=False ,help="make sets for all particles",guitype='boolbox',row=3, col=1, rowspan=1, colspan=1,mode="sets")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	plist=args
	if options.allparticles:
		fld="particles3d/"
		plist=[fld+f for f in os.listdir(fld) if f.endswith(".hdf")]
	
	tags=[f[f.find("__")+2:-4] for f in plist]
	tags=np.unique(tags)
	print("Found particles with following labels: ", ", ".join(tags))
	
	if len(options.label)>0:
		if options.label in tags:
			print("Only build sets with particles with {} label".format(options.label))
			tags=[options.label]
		else:
			print("Cannot find any particles with specified label. Exit.")
			return
	
	for tag in tags:
		fulltg="__{}.hdf".format(tag)
		pts=[f for f in plist if f.endswith(fulltg)]
		out="sets/{}.lst".format(tag)
		try: os.remove(out)
		except: pass
		cmd="e2proclst.py {} --create {}".format(' '.join(pts), out)
		run(cmd)
	
	print("Done")
	
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

if __name__ == '__main__':
	main()

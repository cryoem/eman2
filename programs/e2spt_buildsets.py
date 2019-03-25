#!/usr/bin/env python

# Muyuan Chen 2018-07
from __future__ import print_function
from __future__ import division
from EMAN2 import *
import numpy as np


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particle_stacks",help="Specify particles input.", default="", guitype='filebox', browser="EMSPTParticleTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=2, mode="sets")
	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=3, col=0, rowspan=1, colspan=1, mode="sets")
	parser.add_argument("--label", type=str,help="label of particles for sets", default="", guitype='strbox',row=2, col=0,rowspan=1, colspan=1, mode="sets")
	parser.add_argument("--allparticles", action="store_true", default=False ,help="make sets for all particles",guitype='boolbox',row=3, col=1, rowspan=1, colspan=1,mode="sets")
	parser.add_argument("--spliteo", action="store_true", default=False ,help="split even/odd set so there is no overlap of particles from the two sets.",guitype='boolbox',row=4, col=1, rowspan=1, colspan=1,mode="sets")
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
		if options.spliteo:
			print("Splitting even/odd sets...")
			lsts=[]
			for eo in ["even","odd"]:
				out="sets/{}__{}.lst".format(tag, eo)
				try: os.remove(out)
				except: pass
				lsts.append(LSXFile(out, False))
			
			### need to check location of each particle_stacks
			for fm in pts:
				n=EMUtil.get_image_count(fm)
				p=[]
				imgs=EMData.read_images(fm, (), True)
				p=np.array([e["ptcl_source_coord"] for e in imgs])
				#for i in range(n):
					#e=EMData(fm, i, True)
					#p.append(e["ptcl_source_coord"])
				#p=np.array(p)
				m=np.median(p[:,0])+np.random.randn()*0.01
				for i in range(n):
					lsts[p[i,0]<m].write(-1, i, fm)
					
				print("{} : center: {:.1f}, even {}, odd {}".format(base_name(fm), m, np.sum(p[:,0]<m), np.sum(p[:,0]>m)))
			
			for l in lsts:
				l.close()
	
	print("Done")
	
	E2end(logid)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

if __name__ == '__main__':
	main()

#!/usr/bin/env python

# Author: Steven Ludtke, 06/16/14 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

import os
import sys
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
Will take a movie stackfile with a name corresponding to existing info/*json files, and extract
each particle from each frame, putting all particles into a single stack. If there are N frames
in the movie, each particle will be in the output file N times sequentially, with header data
indicating its position in the movie. 

- 'movie' folder must contain movie frames, each exposure sequentially in one file
- info_name of movie must correspond to info/*json file
- movieparticles folder will be created with corresponding named files
- movie stacks must have at least 3 frames

	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--boxsize", type=int, help="Box size for particle extraction, may be larger than the standard size for the project to compensate for motion",default=-1)
	parser.add_argument("--useset", type=str, help="If a set is specified as input, only the specific particles in that set will be extracted. Warning, the resulting movie-particle stacks will have missing images.",default=None)
	parser.add_argument("--invert",action="store_true",help="Extracted images are multiplied by -1 ",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.boxsize<10 : 
		print "Please specify a box size"
		sys.exit(1)

	box=options.boxsize

	try: os.mkdir("movieparticles")
	except: pass

	logid=E2init(sys.argv)


	uniq={}
	
	if options.useset == None :
		for i in os.listdir("particles"):
			if ".hdf" in i: uniq[base_name(i,nodir=True)]=[]
	else :
		lsx=LSXFile(options.useset,True)
		for i in xrange(len(lsx)): 
			ln=lsx.read(i)
			try: uniq[base_name(ln[1],nodir=True)].append(int(ln[0]))
			except: uniq[base_name(ln[1],nodir=True)]=[int(ln[0])]

	for u in sorted(uniq.keys()):
		if os.path.exists("movie/{}_raw_proc_align.hdf".format(u)): m="{}_raw_proc_align.hdf".format(u)
		elif os.path.exists("movie/{}.hdf".format(u.replace("_aliavg","_align"))): m="{}.hdf".format(u.replace("_aliavg","_align"))
		elif os.path.exists("movie/{}.hdf".format(u)): m="{}_proc_align.hdf".format(u)
		elif os.path.exists("movie/{}.mrcs".format(u)): m="{}.mrcs".format(u)
		elif os.path.exists("movie/"+ u.replace("aliavg","align")+".hdf"):m=u.replace("aliavg","align")+".hdf"
		elif os.path.exists("movie/"+ u.replace("_aliavg","").replace("_proc","")+".mrcs"):m=u.replace("_aliavg","").replace("_proc","")+".mrcs"
		elif os.path.exists("movie/"+ u.split(".")[0]+"_raw_proc_corr.hdf") : m= u.split(".")[0]+"_raw_proc_corr.hdf"
		else :
			print "Couldn't find movie for ",u
			print "Tried: ","{}_raw_proc_align.hdf".format(u), "{}.hdf".format(u.replace("_aliavg","_align")), "{}_proc_align.hdf".format(u), "{}.mrcs".format(u), u.replace("aliavg","align")+".hdf",u.replace("_aliavg","").replace("_proc","")+".mrcs","movie/"+ u.split(".")[0]+"_raw_proc_corr.hdf"
			continue
		print "Movie found {} -> {}".format(u,m)


		fsp=os.path.join("movie",m)
		outfsp=os.path.join("movieparticles",u+"_ptcls.hdf")

		try: 
			try: n=EMUtil.get_image_count(fsp)
			except:
				m=m.replace("_proc_align","")
				fsp=os.path.join("movie",m)
				n=EMUtil.get_image_count(fsp)
			if n<=2 : 
				print "skipping ",m," (too few images)"
				continue	
		except:
			print "skipping ",m," (read error)"
			continue

		#if not os.path.exists("info/{}_info.json".format(u)) :
			#print "No info file for {} (info/{}.json)".format(m,u)
			#continue
		#try:
			#db=js_open_dict(info_name(u))
			#boxes=db["boxes"]
		#except:
			#print "No box locations for {} ({})".format(m,info_name(u))
			#continue

		ptclfile="particles/{}_ptcls.hdf".format(u)
		try: nptcl=EMUtil.get_image_count(ptclfile)
		except: 
			ptclfile="particles/{}.hdf".format(u)
			nptcl=EMUtil.get_image_count(ptclfile)

		if nptcl==0 : 
			print "no particles found :",u
			continue

		if len(uniq[u])==0 : uniq[u]=xrange(nptcl)
		
		for i in xrange(n):
			fm=EMData(fsp,i)
			fm.process_inplace("normalize.edgemean")
			for ib in uniq[u]:
				b=EMData(ptclfile,ib,True)["ptcl_source_coord"]
				ptcl=fm.get_clip(Region(b[0]-box/2,b[1]-box/2,box,box))
				ptcl.process_inplace("mask.zeroedgefill")
				if options.invert : ptcl.mult(-1.0)
				ptcl["movie_frames"]=n
				ptcl["movie_n"]=i
				ptcl.write_image(outfsp,i+ib*n)
	
	E2end(logid)



if __name__ == "__main__":
	main()

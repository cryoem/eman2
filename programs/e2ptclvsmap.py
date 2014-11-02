#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/31/14 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
from EMAN2 import *
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2ptclvsmap.py [options]

This program is used to visually, and somewhat quantitatively, assess how well a set of particles actually agree with a given 3-D map.

Provide a 3-D map, and a set of particles. The program will classify and align all of the particles, matching each with it's best-matching 
projection, very much like the first step in 3-D refinement. Once classified, each particle is compared to its associated projection and
the set of per-projection particles is sorted in order of similarity. An output file is then created containing each projection
followed by a sorted list of associated particles for each.

The reason this script does not just extract similar information from a refine_xx folder is to avoid the map perturbations which occur in
gold-standard refinement."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	#parser.add_header(name="multirefineheader", help='Options below this label are specific to e2refinemulti', title="### e2refinemulti options ###", row=1, col=0, rowspan=1, colspan=3, mode="refinement")
	#parser.add_header(name="multimodelheader", help='Options below this label are specific to e2refinemulti Model', title="### e2refinemulti model options ###", row=4, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image file containing the particle data", guitype='filebox', browser='EMSetsTable(withmodal=True,multiselect=False)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--model", dest="model", type=str,default=None, help="The map to use as a starting point for refinement", guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', filecheck=False, row=3, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--filterout", action="store_true", help="Filters output particles to match projections")
 
	parser.add_argument("--angstep",type=float,default=9.0,help="Angular separation of projections. Default 9.0 degrees.")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos.", guitype='strbox', row=10, col=1, rowspan=1, colspan=1, mode="refinement")

	parser.add_argument("--simalign",type=str,help="Default=auto. The name of an 'aligner' to use prior to comparing the images", default="rotate_translate_flip")
	parser.add_argument("--simaligncmp",type=str,help="Default=auto. Name of the aligner along with its construction arguments",default="ccc")
	parser.add_argument("--simralign",type=str,help="Default=auto. The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_argument("--simraligncmp",type=str,help="Default=auto. The name and parameters of the comparitor used by the second stage aligner.",default="ccc")
	parser.add_argument("--simcmp",type=str,help="Default=auto. The name of a 'cmp' to be used in comparing the aligned images", default="optsub:maxres=12")

	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>=<value>. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel",default=None, guitype='strbox', row=24, col=0, rowspan=1, colspan=2, mode="refinement[thread:4]")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--path", default=None, type=str,help="The name of a directory where results are placed. Default = create new ptclmap_xx")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.threads<=1 :
		if options.parallel!=None and options.parallel[:6]=="thread" : 
			options.threads=int(options.parallel[7:])
			print "Note: automatically setting --threads:{}".format(options.threads)
		else: print "WARNING: specifying --threads=<N> (where N is the number of cores to use on a single processor) is strongly recommended, even if already specifying --parallel"

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:8]=="ptclmap_" and len(i)==10 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.path = "ptclmap_{:02d}".format(max(fls)+1)
	try: os.makedirs(options.path)
	except: pass

	# make sure the box sizes match
	if options.input!=None :
		xsize3d=EMData(options.model,0,True)["nx"]
		xsize=EMData(options.input,0,True)["nx"]
		if ( xsize3d != xsize ) :
			print "WARNING: the dimensions of the particles (%d) do not match the dimensions of the starting model (%d). I will attempt to adjust the model appropriately."%(xsize,xsize3d)
			img1 = EMData(options.input,0,True)
			img3 = EMData(options.model,0,True)
			try:
				scale=img3["apix_x"]/img1["apix_x"]
				print "Reference is {box3} x {box3} x {box3} at {apix3:1.2f} A/pix, particles are {box2} x {box2} at {apix2:1.2f} A/pix. Scaling by {scale:1.3f}".format(box3=img3["nx"],box2=img1["nx"],apix3=img3["apix_x"],apix2=img1["apix_x"],scale=scale)
			except:
				print "A/pix unknown, assuming scale same as relative box size"
				scale=float(xsize)/xsize3d
			if scale>1 : cmd="e2proc3d.py %s %s/scaled_model.hdf --clip=%d,%d,%d --scale=%1.5f"%(options.model,options.path,xsize,xsize,xsize,scale)
			else :       cmd="e2proc3d.py %s %s/scaled_model.hdf --scale=%1.5f --clip=%d,%d,%d"%(options.model,options.path,scale,xsize,xsize,xsize)
			run(cmd)

			options.model="%s/scaled_model.hdf"%options.path

	if options.verbose==0 : verbose=""
	else: verbose="--verbose {}".format(options.verbose)

	if options.parallel!=None : parallel="--parallel {}".format(options.parallel)
	elif options.threads>1: parallel="--parallel thread:{}".format(options.threads)
	else: parallel=""

	# Ok do the classification
	orientgen="eman:delta={:1.5f}:inc_mirror=0:perturb=0".format(options.angstep)
	cmd = "e2project3d.py {map} --outfile {path}/projections.hdf -f --projector standard --orientgen {orient} --sym {sym} --parallel thread:{threads} {verbose}".format(
		path=options.path,map=options.model,orient=orientgen,sym=options.sym,threads=options.threads,verbose=verbose)
	run(cmd)

	if options.simralign==None : simralign=""
	else: simralign="--ralign {} --raligncmp {}".format(options.simralign,options.simraligncmp)

	cmd = "e2simmx2stage.py {path}/projections.hdf {inputfile} {path}/simmx.hdf {path}/proj_simmx.hdf {path}/proj_stg1_.hdf {path}/simmx_stg1.hdf --saveali --cmp {simcmp} \
--align {simalign} --aligncmp {simaligncmp} {simralign} {verbose} {parallel}".format(
		path=options.path,inputfile=options.input,simcmp=options.simcmp,simalign=options.simalign,simaligncmp=options.simaligncmp,simralign=simralign,
		verbose=verbose,parallel=parallel)
	run(cmd)

	cmd = "e2classify.py {path}/simmx.hdf {path}/classmx.hdf -f {verbose}".format(
		path=options.path,verbose=verbose)
	run(cmd)

	# Read in the final classification info and process the particles
	pathmx="{path}/classmx.hdf".format(path=options.path)
	classmx=EMData(pathmx,0)
	cmxtx=EMData(pathmx,2)
	cmxty=EMData(pathmx,3)
	cmxalpha=EMData(pathmx,4)
	cmxmirror=EMData(pathmx,5)

	nref=int(classmx["maximum"])+1
	nptcl=classmx["ny"]

	for iref in xrange(nref):
		if options.verbose==1 : print "Class ",iref
		outname="{}/class_{:04d}.hdf".format(options.path,iref)
		ref=EMData("{path}/projections.hdf".format(path=options.path),iref)
		ref.process_inplace("normalize.edgemean")
		ref.write_image(outname,0)

		allptcl=[]		
		for iptcl in xrange(nptcl):
			if classmx[0,iptcl]!=iref : continue		# only proceed if the particle is in this class

			ptcl=EMData(options.input,iptcl)
			ptcl.process_inplace("normalize.edgemean")

			# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
			ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,iptcl],"mirror":int(cmxmirror[0,iptcl]),"tx":cmxtx[0,iptcl],"ty":cmxty[0,iptcl]})
			ptclx=ptcl.process("xform",{"transform":ptclxf})
			ptclx.process_inplace("mask.sharp",{"outer_radius":ptcl["ny"]/2.5})
			c=ptclx.cmp("optsub",ref,{"maxres":18.0})
			if options.filterout : ptclx=ref.process("math.sub.optimal",{"return_subim":1,"ref":ptclx})

			allptcl.append((c,ptclx))
		
		for q,p in sorted(allptcl):
			p["qual"]=q
			p.write_image(outname,-1)
			




def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print "{}: {}".format(time.ctime(time.time()),command)

	ret=launch_childprocess(command)

	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)

	return

if __name__ == "__main__":
    main()


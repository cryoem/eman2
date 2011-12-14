#!/usr/bin/env python

#
# Author: Jesus Galaz, 12/08/2011
# Copyright (c) 2011 Baylor College of Medicine
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

from EMAN2 import *
from sys import argv
import os
from e2classaverage3d import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <stack of particles> [options] . This programs averages a set of particles, performs the spherical average of the average, 
		and translationally aligns all the particles to it to re-center them."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--clip", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	parser.add_argument("--mask",type=str,default="mask.sharp:outer_radius=-2",help="Mask to apply to particles before spherically averaging each and centering them")
	parser.add_argument("--shrink",type=int,default=1,help="Optionally shrink the volume(s) to have the FFTs go faster")
	parser.add_argument("--preprocess",type=str,help="A processor (as in e2proc3d.py) to be applied to each volume prior to alignment", default=None)
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean.tomo")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False)
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False)
	#parser.add_argument("--align", help="Set by default",default="translational:intonly=1")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	
	data = argv[1]

	hdr = EMData(data,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
	box=nx
	
	parser.add_argument("--align", help="Set by default",default="rotate_translate_3d_grid:alt0=0:alt1=1:az0=0:az1=1:dalt=2:daz=2:dotrans=1:dphi=2:phi0=0:phi1=1:search=10:verbose=1")

	(options, args) = parser.parse_args()

	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)
		
	if options.path and options.path[:4].lower()!="bdb:": 
		options.path="bdb:"+options.path
	
	if not options.path: 
		options.path="bdb:"+numbered_path("sptautoc",True)

	if options.averager: 
		options.averager=parsemodopt(options.averager)
	if options.mask: 
		options.mask=parsemodopt(options.mask)

	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)

	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	
	if options.align: 
		options.align=parsemodopt(options.align)
	
	if options.aligncmp: 
		options.aligncmp=parsemodopt(options.aligncmp)
	
	n=EMUtil.get_image_count(data)
	
	for i in range(options.iter):
		avgr=Averagers.get(options.averager[0], options.averager[1])			#Call the averager
		if i == 0:
			for k in range(n):
				a=EMData(data,k)
				avgr.add_image(a)
			avg=avgr.finish()
		
		# mask
		mask = EMData(int(hdr['nx']),int(hdr['ny']),int(hdr['nz']))
		mask.to_one()
		
		if options.mask != None:
			mask.process_inplace(options.mask[0],options.mask[1])
	
		# normalize
		if options.normproc != None:
			if options.normproc[0]=="normalize.mask": 
				options.normproc[1]["mask"]=mask
			avg.process_inplace(options.normproc[0],options.normproc[1])

		# preprocess
		if options.preprocess != None:
			avg.process_inplace(options.preprocess[0],options.preprocess[1])
			
		# Shrinking both for initial alignment and reference
		if options.shrink!=None and options.shrink > 1 :
			avg=avg.process("math.meanshrink",{"n":options.shrink})

		avg_sph=avg.rotavg_i()
		
		if options.savesteps :
			avg.write_image("%s#averages" % (options.path),-1)
			avg_sph.write_image("%s#averages_sph" % (options.path),-1)

		# Initialize parallelism if being used
		if options.parallel:
			from EMAN2PAR import EMTaskCustomer
			etc=EMTaskCustomer(options.parallel)
			pclist=[data]
		
		tasks = []
		for j in range(n):
			task=Align3DTask(avg_sph,["cache",data,j],j,"Ptcl %d in iter %d"%(j,i),options.mask,options.normproc,options.preprocess,
					1,options.align,options.aligncmp,None,None,1,1,None,options.verbose-1)
			#task=Align3DTask(["cache",infile,j],["cache",infile,j+1],j/2,"Seed Tree pair %d at level %d"%(j/2,i),options.mask,options.normproc,options.preprocess,
			#options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,transform,options.verbose-1)		
				
			tasks.append(task)

		tids=etc.send_tasks(tasks)						#Start the alignments running
		
		if options.verbose > 0:
			print "%d tasks queued in iteration %d"%(len(tids),k) 
		
		results = get_results(etc,tids,options.verbose)				#Wait for alignments to finish and get results
				
		#avgr=Averagers.get(options.averager[0], options.averager[1])			#Call the averager

		for z in range(len(results)):					
			t = results[z][0]['xform.align3d']
			ptcl = EMData(data,z)
			ptcl.process_inplace('xform',{'transform':t})			
			avgr.add_image(ptcl)
			
			if options.saveali:
				ptcl.write_image('%s#particles_round%02d' % (options.path,i),-1)
			
			if not options.saveali and i == options.iter -1:
				ptcl.write_image('%s#particles' % (options.path),-1)
		
		avg=avgr.finish()
	
	return()

if __name__ == "__main__":
	main()

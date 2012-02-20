#!/usr/bin/env python

#
# Author: John Flanagan Feb 2012 (jfflanag@bcm.edu)
# Copyright (c) 2012- Baylor College of Medicine
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


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program aligns a stack of 3D particles their symmetry axes and then averages them. By aligning to the symmetry axis, the particles are aligned to themselves.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_header(name="alignheader", help='Options below this label are specific to e2refine', title="### e2classaveragebysym options ###", row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of input stack of volumes", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=2)
	parser.add_argument("--output", dest="output", default=None,type=str, help="The name of the aligned and averaged output volume", guitype='strbox', filecheck=False, row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--sym", dest = "sym", default="c1", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction omit this option or specify c1.", guitype='symbox', row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=6, col=0, rowspan=1, colspan=2)
	parser.add_argument("--preprocess",type=str,help="A processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=5, col=0, rowspan=1, colspan=2)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Default=0, no shrinking", guitype='shrinkbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--steps", dest="steps", type = int, default=10, help="Number of steps (for the MC). This should be a multiple of the number of cores used for parallization", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--symmetrize", default=True, action="store_true", help="Symmetrize volume after alignment.", guitype='boolbox', row=7, col=0, rowspan=1, colspan=1)
	parser.add_argument("--applytoraw", default=True, action="store_true", help="Applies symxform to raw data.", guitype='boolbox', row=7, col=1, rowspan=1, colspan=1)
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the symmtrized object to unsymmetrized", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=8, col=0, rowspan=1, colspan=2)
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean", guitype='combobox', choicelist='dump_averagers_list()', row=9, col=0, rowspan=1, colspan=2)
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=10, col=0, rowspan=1, colspan=2)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	if options.averager: 
		options.averager=parsemodopt(options.averager)
		
	# Check for errors
	hdr = EMData(options.input,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
		
	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	# get numbered path
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)
		
	if options.path and options.path[:4].lower()!="bdb:": 
		options.path="bdb:"+options.path
	
	if not options.path: 
		options.path="bdb:"+numbered_path("sptsym",True)
	
	# Get the averager
	avgr=Averagers.get(options.averager[0], options.averager[1])
	
	# generate the mask
	mask=EMData(nx,ny,nz)
	mask.to_one()
	if options.mask != None:
		print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
		mask.process_inplace(options.mask[0],options.mask[1])
			
	# Align each particle to its symmetry axis
	for i in xrange(EMUtil.get_image_count(options.input)):
		# Copy image
		inputfile = "%s#ptcl_to_align_%d"%(options.path,i)
		outputfile = "%s#aligned_ptcl_%d"%(options.path,i)
		model3d = EMData(options.input, i)
		
		# apply mask if desired
		model3d.mult(mask)
		
		# preprocess
		if options.preprocess != None:
			model3d.process_inplace(options.preprocess[0],options.preprocess[1])
			
		# write out file	
		model3d.write_image(inputfile)
		
		command = "e2symsearch.py --input=%s --output=%s --sym=%s --shrink=%d --steps=%d --cmp=%s --path=''"%(inputfile, outputfile, options.sym, options.shrink, options.steps, options.cmp)
		
		if options.symmetrize:
			command += " --symmetrize"
		if options.parallel:
			command += " --parallel=%s"%options.parallel
		launch_childprocess(command)
		db_remove_dict(inputfile)
		
		if options.applytoraw:
			rawoutputfile = "%s#rawaligned_ptcl_%d"%(options.path,i)
			symxform = EMData(outputfile).get_attr('symxform')
			raw3d = EMData(options.input, i).process_inplace('xform',{'transform':symxform})
			raw3d.write_image(rawoutputfile)
			
		# Add the volume to the avarage
		avgr.add_image(EMData(outputfile))
		
	# Now make the avearage
	average = avgr.finish()
	average.write_image(options.output)
		
if __name__ == "__main__":
    main()
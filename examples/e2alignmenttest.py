#!/usr/bin/env python

#
# Author: David Woolford, 9/7/2007 (woolford@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
from optparse import OptionParser
from math import *
from copy import deepcopy
import os
import sys

READ_HEADER_ONLY = True

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog 3d_model [options]

	"""
		
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--noise_processor",type="string",help="The processor used to add noise to the images", default="math.addnoise:noise=2")
	parser.add_option("--align",type="string",help="The first aligner used", default="rotate_translate_flip")
	parser.add_option("--aligncmp",type="string",help="The comparitor used for the --align aligner. Default is dot:normlize=1.",default="dot:normalize=1")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner. If not specified then only one stage of alignment is applied.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparitor used by the second stage aligner. Default is dot:normlize=1.",default="dot:normalize=1")
	
	parser.add_option("--transmax",type="float",help="The maximum translational distance in any direction. Default is 5.0",default=5.0)
	parser.add_option("--transmin",type="float",help="The minimum translational distance in any direction. Default is -5.0",default=-5.0)
	
	parser.add_option("--rotmax",type="float",help="The maximum rotation angle applied. Default is 360.0",default=360.0)
	parser.add_option("--rotmin",type="float",help="The minimum rotation angle applied. Default is 0.0",default=0.0)
	
	parser.add_option("--allowflip","-a",action="store_true",help="Allow randomized flipping. Default is off.",default=False)
	
	parser.add_option("--num",type="int",help="The number of randomized tests to perform. Default is 20",default=20)

	parser.add_option("--projector", dest = "projector", default = "standard",help = "Projector to use. Default is standard.")
	
	(options, args) = parser.parse_args()
	
	options.model = args[0]
	
	if not os.path.exists(options.model):
		print "Error: 3D image %s does not exist" %options.model
		exit(1)
		
	model = EMData(options.model)
	if model.get_ndim() != 3:
		print "error, the model you specified is not 3D"
		
	c1_sym = Symmetries.get("c1")
	rand_orient = OrientGens.get("rand",{"n":1})
	[noise_proc,noise_params] = parsemodopt(options.noise_processor)
	[aligner,aligner_params] = parsemodopt(options.align)
	[alignercmp,alignercmp_params] = parsemodopt(options.aligncmp)
	if options.ralign != None:
		[raligner,raligner_params] = parsemodopt(options.ralign)
		[alignercmp,ralignercmp_params] = parsemodopt(options.raligncmp)
	
	az_error = 0.0
	dx_error = 0.0
	dy_error = 0.0
	if options.allowflip:
		flip_errors = 0
	for i in range(0,options.num):
		t3d = rand_orient.gen_orientations(c1_sym)[0]
		d = {"t3d":t3d}
		p=model.project(options.projector,d)
		p.process_inplace(noise_proc,noise_params)
		q = p.copy()
		az = Util.get_frand(options.rotmin,options.rotmax)
		dx = Util.get_frand(options.transmin,options.transmax)
		dy = Util.get_frand(options.transmin,options.transmax)
		t3d_q = Transform3D(az,0,0)
		t3d_q.set_pretrans(dx,dy,0)
		
		if options.allowflip:
			flipped = Util.get_irand(0,1)
			if flipped: q.process_inplace("xform.flip",{"axis":"x"})
			
		q.rotate_translate(t3d_q)
		#print aligner,p,aligner_params, alignercmp, alignercmp_params
		ali = q.align(aligner,p,aligner_params, alignercmp, alignercmp_params)
		
		az_solution = (-ali.get_attr("align.az"))%360
		dx_solution = -ali.get_attr("align.dx")
		dy_solution = -ali.get_attr("align.dy")
		if options.allowflip:
			if ali.get_attr("align.flip"):
				az_solution = ali.get_attr("align.az")
				
		print "%.2f"%az,'\t',"%.2f"%az_solution, '\t\t', "%.2f"%dx,'\t',"%.2f"%(dx_solution),'\t\t', "%.2f"%dy,'\t',"%.2f"%(dy_solution),
		if options.allowflip: print '\t\t',flipped, '\t',ali.get_attr("align.flip")
		else: print ''
		
		az_error = fabs(az-az_solution)
		t = fabs(az_error -360)
		if t  < az_error: az_error = t
		az_error += fabs(az-az_solution)
		dx_error += fabs(dx-dx_solution)
		dy_error += fabs(dy-dy_solution)
		
		if options.allowflip:
			if flipped != ali.get_attr("align.flip"):
				flip_errors += 1
		
	print "#### REPORT ####"
	print "Average az error",az_error/options.num
	print "Average dx error",dx_error/options.num
	print "Average ay error",dy_error/options.num
	if options.allowflip:
		print "Flip detection accuracy", float(options.num-flip_errors)/options.num*100,"%"
	
		
if __name__ == "__main__":
    main()

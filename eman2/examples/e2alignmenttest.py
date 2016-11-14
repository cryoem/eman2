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

	parser.add_option("--noise_processor",type="string",help="The processor used to add noise to the images. Default is math.addnoise:noise=50.", default="math.addnoise:noise=20")
	parser.add_option("--align",type="string",help="The first aligner used. Default is rotate_translate_flip:rfp_mode=0", default="rotate_translate_flip:rfp_mode=0")
	parser.add_option("--aligncmp",type="string",help="The comparitor used for the --align aligner. Default is dot:normalize=1.",default="dot:normalize=1")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This must currently be the \'refine\' aligner. If not specified then only one stage of alignment is applied. Default is None.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparitor used by the second stage aligner. Default is dot:normlize=1.",default="dot:normalize=1")
	parser.add_option("--transmax",type="float",help="The maximum translational distance in any direction. Default is 5.0",default=5.0)
	parser.add_option("--transmin",type="float",help="The minimum translational distance in any direction. Default is -5.0",default=-5.0)
	parser.add_option("--rotmax",type="float",help="The maximum rotation angle applied. Default is 360.0",default=360.0)
	parser.add_option("--rotmin",type="float",help="The minimum rotation angle applied. Default is 0.0",default=0.0)
	
	parser.add_option("--stopflip","-s",action="store_true",help="Stop randomized flipping. Default is off.",default=False)
	
	parser.add_option("--num",type="int",help="The number of randomized tests to perform. Default is 20",default=20)

	parser.add_option("--projector", dest = "projector", default = "standard",help = "Projector to use. Default is standard.")
	
	parser.add_option("--writeout","-w",action="store_true",help="Write the projection + noise images to projections_and_noise.hdf. Default is off.",default=False)
	
	(options, args) = parser.parse_args()
	
	options.model = args[0]
	
	if not os.path.exists(options.model):
		print "Error: 3D image %s does not exist" %options.model
		exit(1)
		
	model = EMData(options.model)
	if model.get_ndim() != 3:
		print "error, the model you specified is not 3D"
		
	# use c1 symmetry to generate random orienations on the unit sphere
	c1_sym = Symmetries.get("c1")
	rand_orient = OrientGens.get("rand",{"n":1}) # this is the orientation generator
	# Get all necessary parameters
	[noise_proc,noise_params] = parsemodopt(options.noise_processor)
	[aligner,aligner_params] = parsemodopt(options.align)
	[alignercmp,alignercmp_params] = parsemodopt(options.aligncmp)
	
	# if refine aligning then get the parameters
	if options.ralign != None:
		[raligner,raligner_params] = parsemodopt(options.ralign)
		[ralignercmp,ralignercmp_params] = parsemodopt(options.raligncmp)
	
	# error variables
	az_error = 0.0
	dx_error = 0.0
	dy_error = 0.0
	flip_errors = 0
	if options.ralign:
		refine_az_error = 0.0
		refine_dx_error = 0.0
		refine_dy_error = 0.0
		
	# print stuff
	print "EULER_AZ, ALT,PHI",'\t',
	print "AZ",'\t',"ALIGN.AZ","\t","DX","\t","ALIGN.DX","\t","DY","\t","ALIGN.DY",
	if not options.stopflip:
		print "\t","FLIP","\t","ALIGN.FLIP"
	else : print ''
	
	# finally the main loop
	for i in range(0,options.num):
		# get the random orientation on the unit sphere
		t3d = rand_orient.gen_orientations(c1_sym)[0]
		d = {"transform":t3d}
		p=model.project(options.projector,d) # make the projection
		p.process_inplace(noise_proc,noise_params) # add noise
		if options.writeout: p.write_image('projections_and_noise.hdf',-1)
			
		q = p.copy()
		az = Util.get_frand(options.rotmin,options.rotmax) # make random angle
		dx = Util.get_frand(options.transmin,options.transmax) # random dx
		dy = Util.get_frand(options.transmin,options.transmax) # random dy
		
		if not options.stopflip:
			flipped = Util.get_irand(0,1)
			if flipped: q.process_inplace("xform.flip",{"axis":"x"})
			
		t3d_q = Transform({"type":"2d","alpha":az})
		t3d_q.set_pre_trans([dx,dy,0.0])
		
		t = Transform()
		try: t.set_mirror(flipped)
		except: pass
		t = t3d_q * t
			
		q.transform(t3d_q)
		#Do the alignment now
		ali = q.align(aligner,p,aligner_params, alignercmp, alignercmp_params)
		alit = ali.get_attr("xform.align2d")
		soln_parms = alit.get_params("2d")
		az_solution = (-soln_parms["alpha"])%360
		dx_solution = -soln_parms["tx"]
		dy_solution = -soln_parms["ty"]
		
		# print stuff
		d = t3d.get_rotation("eman")
		print "%.2f,%.2f,%.2f\t"%(d["az"],d["alt"],d["phi"]),
		print "%.2f"%az,'\t',"%.2f"%az_solution, '\t\t', "%.2f"%dx,'\t',"%.2f"%(dx_solution),'\t\t', "%.2f"%dy,'\t',"%.2f"%(dy_solution),
		if not options.stopflip: print '\t\t',flipped, '\t',int(alit.get_mirror()),
		
	
		# calculate the errors
		error = fabs(az-az_solution)
		t = fabs(error -360) # could be close going counter clockwise etc
		if t  < error: error = t
		az_error += fabs(error)
		dx_error += fabs(dx-dx_solution)
		dy_error += fabs(dy-dy_solution)
		print ''
		if not options.stopflip:
			if flipped != int(alit.get_mirror()):
				flip_errors += 1
				print "FLIP detection FAILED"
				continue
		
		# refine align if it has been specified
		if options.ralign:
			
			#raligner_params["az"] = ali.get_attr_default("align.az",0)
			#raligner_params["dx"] = ali.get_attr_default("align.dx",0)
			#raligner_params["dy"] = ali.get_attr_default("align.dy",0)
			#flip = ali.get_attr_default("align.flip",0)

			q = p.copy()
			q.set_attr("xform.align2d",ali.get_attr("xform.align2d"))
			q.transform(t3d_q)
			#q.write_image("pq.hdf",-1)
			#p.write_image("pq.hdf",-1)
			
			
			ali = q.align(raligner,p,raligner_params,ralignercmp,ralignercmp_params)
			alit = ali.get_attr("xform.align2d")
			soln_parms = alit.get_params("2d")
			az_solution = (-soln_parms["alpha"])%360
			dx_solution = -soln_parms["tx"]
			dy_solution = -soln_parms["ty"]
			#az_solution = (-ali.get_attr("align.az"))%360
			#dx_solution = -ali.get_attr("align.dx")
			#dy_solution = -ali.get_attr("align.dy")
			
			#if not options.stopflip:
				#if flip:
					#az_solution = ali.get_attr("align.az")
					#dx_solution = -dx_solution
					
			error = fabs(az-az_solution)
			t = fabs(error -360)
			if t  < error: error = t
			refine_az_error += fabs(error)
			refine_dx_error += fabs(dx-dx_solution)
			refine_dy_error += fabs(dy-dy_solution)
			
			print "REFINE           ",'\t',
			print "%.2f"%az,'\t',"%.2f"%az_solution, '\t\t', "%.2f"%dx,'\t',"%.2f"%(dx_solution),'\t\t', "%.2f"%dy,'\t',"%.2f"%(dy_solution)
		
		
	print "#### REPORT ####"
	print "Mean az error",az_error/options.num
	print "Mean dx error",dx_error/options.num
	print "Mean dy error",dy_error/options.num
	if not options.stopflip:
		print "Flip detection accuracy", float(options.num-flip_errors)/options.num*100,"%"
	if options.ralign:
		print "Mean refine az error",refine_az_error/float(options.num-flip_errors)
		print "Mean refine dx error",refine_dx_error/float(options.num-flip_errors)
		print "Mean refine dy error",refine_dy_error/float(options.num-flip_errors)
if __name__ == "__main__":
    main()

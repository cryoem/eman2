#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

# $Id$

from EMAN2 import *
import sys
import os.path
import math
import random
import pyemtbx.options
import time
from random   import random, seed, randint

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>

	Generic 2-D image processing and file format conversion program. Acts on stacks of 2-D images
	(multiple images in one file). 

	Examples:

	phase flip a stack of images and write output to new file:
	sxprocess.py input_stack.hdf output_stack.hdf --phase_flip	
	
	generate stack of projections:
	sxprocess.py --generate_projections format="bdb":apix=5.2:CTF=True 	
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--phase_flip", action="store_true", help="Phase flip the input stack", default=False)
	parser.add_argument("--makedb", type=str, help="Fill in database with appropriate input parameters: --makedb=mpibdb means the input parameters will be those in mpi_bdb, --makedb=mpibdbctf means the input parameters will be those in mpi_bdb_ctf", default=None)
	parser.add_argument("--generate_projections", metavar="param1=value1:param2=value2", type=str,
					action="append", help="apply a processor named 'processorname' with all its parameters/values.")
	(options, args) = parser.parse_args()
	
	
	if options.phase_flip:
		nargs = len(args)
		if nargs != 2:
			print "must provide name of input and output file!"
			return
		
		instack = args[0]
		outstack = args[1]
		print "input stack: ", instack
		print "output (phase flipped) stack: ", outstack
		nima = EMUtil.get_image_count(instack)
		from filter import filt_ctf
		for i in xrange(nima):
			img = EMData()
			img.read_image(instack, i)
			try:
				ctf = img.get_attr('ctf')
			except:
				print "no ctf information in input stack! Exiting..."
				return
			
			dopad = True
			sign = 1
			binary = 1 # phase flip
				
			assert img.get_ysize() > 1	
			dict = ctf.to_dict()
			dz = dict["defocus"]
			cs = dict["cs"]
			voltage = dict["voltage"]
			pixel_size = dict["apix"]
			b_factor = dict["bfactor"]
			ampcont = dict["ampcont"]
			dza = dict["dfdiff"]
			azz = dict["dfang"]
			
			if dopad and not img.is_complex():  ip = 1
			else:                             ip = 0
	
	
			params = {"filter_type": Processor.fourier_filter_types.CTF_,
	 			"defocus" : dz,
				"Cs": cs,
				"voltage": voltage,
				"Pixel_size": pixel_size,
				"B_factor": b_factor,
				"amp_contrast": ampcont,
				"dopad": ip,
				"binary": binary,
				"sign": sign,
				"dza": dza,
				"azz":azz}
			
			tmp = Processor.EMFourierFilter(img, params)
			tmp.set_attr_dict({"ctf":ctf})
			
			tmp.write_image(outstack, i)
			
	if options.makedb != None:
		nargs = len(args)
		if nargs != 1:
			print "must provide exactly one argument denoting database key under which the input params will be stored"
			return
		dbkey = args[0]
		print "database key under which params will be stored: ", dbkey
		gbdbname = 'bdb:e2boxercache#gauss_box_DB'
		gbdb = db_open_dict(gbdbname)
		if options.makedb == 'mpibdbctf':
			gbdb[dbkey] = {'gauss_width':1.0,'pixel_input':5.2,'pixel_output':5.2,'thr_low':1.0,'thr_hi':30.0,"invert_contrast":False,"use_variance":True,"boxsize":64,"ctf_cs":2.0,"ctf_fstart":0.02,"ctf_fstop":0.5,"ctf_ampcont":10,"ctf_volt":120,"ctf_window":512,"ctf_edge":0,"ctf_overlap":50}
		if options.makedb == 'mpibdb':
			gbdb[dbkey] = {'gauss_width':1.0,'pixel_input':5.2,'pixel_output':5.2,'thr_low':4.0,'thr_hi':60.0,"invert_contrast":False,"use_variance":True,"boxsize":64}
	
	if options.generate_projections:
		parmstr= 'dummy:'+options.generate_projections[0]
		(processorname, param_dict) = parsemodopt(parmstr)
		
		parm_CTF = False
		parm_format = 'bdb'
		parm_apix = 2.5
		
		if 'CTF' in param_dict:
			if param_dict['CTF'] == 'True':
				parm_CTF = True
		
		if 'format' in param_dict:
			parm_format = param_dict['format'] 
		
		if 'apix' in param_dict:
			parm_apix = float(param_dict['apix'])
						
		print "pixel size: ", parm_apix," format: ", parm_format," add CTF: ",parm_CTF
		from filter import filt_gaussl, filt_ctf
		from utilities import drop_spider_doc, even_angles, model_gauss, delete_bdb, model_blank,pad,model_gauss_noise,set_params2D, set_params_proj
		from projection import prep_vol,prgs
		seed(14567)
		delta = 40
		angles = even_angles(delta,0.0,89.9,0.0,359.9,"P")
		nangle = len(angles)
		
		modelvol = EMData()
		modelvol.read_image("../model_structure.hdf")
		
		nx = modelvol.get_xsize()

		nvol = 10
		volfts = [None]*nvol
		for i in xrange(nvol):
			sigma = 1.5 + random() # 1.5-2.5
			addon = model_gauss(sigma, 64, 64, 64, sigma, sigma, 38, 38, 40 )
			scale = 2500 * (0.5+random())
			model = modelvol + scale*addon
			volfts[i],kb = prep_vol(modelvol + scale*addon)
			
		if parm_format=="bdb":
			stack_data = "bdb:data"
			delete_bdb(stack_data)
		else:
			stack_data = "data.hdf"
		Cs      = 2.0
		pixel   = parm_apix
		voltage = 120.0
		ampcont = 10.0
		ibd     = 4096/2-64
		iprj    = 0

		width = 240
		xstart = 8 + 64/2
		ystart = 8 + 64/2
		rowlen = 17

		params = []
		for idef in xrange(3,8):

			irow = 0
			icol = 0

			mic = model_blank(4096, 4096)
			defocus = idef*0.2
			if parm_CTF :
				ctf = EMAN2Ctf()
				ctf.from_dict( {"defocus":defocus, "cs":Cs, "voltage":voltage, "apix":pixel, "ampcont":ampcont, "bfactor":0.0} )
		

			for i in xrange(nangle):
				for k in xrange(24):
					dphi = 8.0*(random()-0.5)
					dtht = 8.0*(random()-0.5)
					psi  = 360.0*random()

					phi = angles[i][0]+dphi
					tht = angles[i][1]+dtht

					s2x = 4.0*(random()-0.5)
					s2y = 4.0*(random()-0.5)


					params.append( [phi, tht, psi, s2x, s2y])

					ivol = iprj%nvol
					proj = prgs(volfts[ivol], kb, [phi, tht, psi, -s2x, -s2y])
		
					x = xstart + irow * width
					y = ystart + icol * width
	
					mic += pad(proj, 4096, 4096, 1, 0.0, x-2048, y-2048, 0)
			
					proj = proj + model_gauss_noise( 30.0, nx, nx )
					if parm_CTF :
						proj = filt_ctf(proj, ctf)
						proj.set_attr_dict({"ctf":ctf, "ctf_applied":0})
			
					proj = proj + filt_gaussl(model_gauss_noise(17.5,nx,nx), 0.3)
					proj.set_attr( "active", 1 )
					# flags describing the status of the image (1 = true, 0 = false)
					set_params2D(proj, [0.0, 0.0, 0.0, 0, 1.0])
					set_params_proj(proj, [phi, tht, psi, s2x, s2y])

					proj.write_image(stack_data, iprj)
			
					icol += 1
					if icol == rowlen:
						icol = 0
						irow += 1

					iprj += 1

			mic += model_gauss_noise(30.0,4096,4096)
			if parm_CTF :
				#apply CTF
				mic = filt_ctf(mic, ctf)
			mic += filt_gaussl(model_gauss_noise(17.5,4096,4096), 0.3)
	
			mic.write_image("mic%1d.hdf"%(idef-3),0)
		
		drop_spider_doc("params.txt", params)
	
			
if __name__ == "__main__":
	main()

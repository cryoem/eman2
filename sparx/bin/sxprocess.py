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
	
	generate a stack of projections bdb:data and micrographs with prefix mic (i.e., mic0.hdf, mic1.hdf etc) from structure input_structure.hdf, with CTF applied to both projections and micrographs:
	sxprocess.py input_structure.hdf data mic --generate_projections format="bdb":apix=5.2:CTF=True:boxsize=64 	
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--order", action="store_true", help="Two arguments are required: name of input stack and desired name of output stack. The output stack is the input stack sorted by similarity in terms of cross-correlation coefficent.", default=False)
	parser.add_argument("--phase_flip", action="store_true", help="Phase flip the input stack", default=False)
	parser.add_argument("--makedb", metavar="param1=value1:param2=value2", type=str,
					action="append",  help="One argument is required: name of key with which the database will be created. Fill in database with parameters specified as follows: --makedb param1=value1:param2=value2, e.g. 'gauss_width'=1.0:'pixel_input'=5.2:'pixel_output'=5.2:'thr_low'=1.0")
	parser.add_argument("--generate_projections", metavar="param1=value1:param2=value2", type=str,
					action="append", help="Three arguments are required: name of input structure from which to generate projections, desired name of output projection stack, and desired prefix for micrographs (e.g. if prefix is 'mic', then micrographs mic0.hdf, mic1.hdf etc will be generated). Optional arguments specifying format, apix, box size and whether to add CTF effects can be entered as follows after --generate_projections: format='bdb':apix=5.2:CTF=True:boxsize=100, or format='hdf', etc., where format is bdb or hdf, apix (pixel size) is a float, CTF is True or False, and boxsize denotes the dimension of the box (assumed to be a square). If an optional parameter is not specified, it will default as follows: format='bdb', apix=2.5, CTF=False, boxsize=64.")
	(options, args) = parser.parse_args()
	
	if options.order:
		nargs = len(args)
		if nargs != 2:
			print "must provide name of input and output file!"
			return
		
		from utilities import get_params2D, model_circle
		from fundamentals import rot_shift2D
		from statistics import ccc
		stack = args[0]
		new_stack = args[1]
		print "input stack: ", stack
		print "output stack: ", new_stack
		
		d = EMData.read_images(stack)
		for i in xrange(len(d)):
       			alpha, sx, sy, mirror, scale = get_params2D(d[i])
        		d[i] = rot_shift2D(d[i], alpha, sx, sy, mirror)
		m = model_circle(30, 64, 64)

		init = 0
		temp = d[init].copy()
		temp.write_image(new_stack, 0)
		del d[init]
		k = 1
		while len(d) > 1:
        		maxcit = -111.
        		for i in xrange(len(d)):
                		cuc = ccc(d[i], temp, m)
                		if cuc > maxcit:
                        		maxcit = cuc
                        		qi = i
        		#print k, maxcit
        		temp = d[qi].copy()
        		del d[qi]
        		temp.write_image(new_stack, k)
        		k += 1

		d[0].write_image(new_stack, k)
			
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
				
		parmstr= 'dummy:'+options.makedb[0]
		(processorname, param_dict) = parsemodopt(parmstr)
		dbdict = {}
		for pkey in param_dict:
			if (pkey == 'invert_contrast') or (pkey == 'use_variance'):
				if param_dict[pkey] == 'True':
					dbdict[pkey] = True
				else:
					dbdict[pkey] = False
			else:		
				dbdict[pkey] = param_dict[pkey]
		gbdb[dbkey] = dbdict
		
	if options.generate_projections:
	
		nargs = len(args)
		if nargs != 3:
			print "must provide name of input structure from which to generate projections, desired name of output projection stack, and desired prefix for output micrographs. Exiting..."
			return
		inpstr = args[0]
		outstk = args[1]
		micpref = args[2]
		print 'input structure: ', inpstr
		print 'output projection stack: ', outstk
		print 'micrograph prefix: ', micpref
			
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
		
		boxsize = 64
		if 'boxsize' in param_dict:
			boxsize = int(param_dict['boxsize'])
							
		print "pixel size: ", parm_apix," format: ", parm_format," add CTF: ",parm_CTF," box size: ", boxsize
		
		scale_mult = 2500
		sigma_add = 1.5
		sigma_proj = 30.0
		sigma2_proj = 17.5
		sigma_gauss = 0.3
		sigma_mic = 30.0
		sigma2_mic = 17.5
		sigma_gauss_mic = 0.3
		
		if 'scale_mult' in param_dict:
			scale_mult = float(param_dict['scale_mult'])
		if 'sigma_add' in param_dict:
			sigma_add = float(param_dict['sigma_add'])
		if 'sigma_proj' in param_dict:
			sigma_proj = float(param_dict['sigma_proj'])
		if 'sigma2_proj' in param_dict:
			sigma2_proj = float(param_dict['sigma2_proj'])
		if 'sigma_gauss' in param_dict:
			sigma_gauss = float(param_dict['sigma_gauss'])	
		if 'sigma_mic' in param_dict:
			sigma_mic = float(param_dict['sigma_mic'])
		if 'sigma2_mic' in param_dict:
			sigma2_mic = float(param_dict['sigma2_mic'])
		if 'sigma_gauss_mic' in param_dict:
			sigma_gauss_mic = float(param_dict['sigma_gauss_mic'])	
			
		from filter import filt_gaussl, filt_ctf
		from utilities import drop_spider_doc, even_angles, model_gauss, delete_bdb, model_blank,pad,model_gauss_noise,set_params2D, set_params_proj
		from projection import prep_vol,prgs
		seed(14567)
		delta = 40
		angles = even_angles(delta,0.0,89.9,0.0,359.9,"P")
		nangle = len(angles)
		
		modelvol = EMData()
		#modelvol.read_image("../model_structure.hdf")
		modelvol.read_image(inpstr)
		
		nx = modelvol.get_xsize()
		
		if nx != boxsize:
			print "requested box dimension does not match dimension of the input model....Exiting"
			sys.exit()
			
		nvol = 10
		volfts = [None]*nvol
		for i in xrange(nvol):
			sigma = sigma_add + random() # 1.5-2.5
			addon = model_gauss(sigma, boxsize, boxsize, boxsize, sigma, sigma, 38, 38, 40 )
			scale = scale_mult * (0.5+random())
			model = modelvol + scale*addon
			volfts[i],kb = prep_vol(modelvol + scale*addon)
			
		if parm_format=="bdb":
			stack_data = "bdb:"+outstk
			delete_bdb(stack_data)
		else:
			stack_data = outstk + ".hdf"
		Cs      = 2.0
		pixel   = parm_apix
		voltage = 120.0
		ampcont = 10.0
		ibd     = 4096/2-boxsize
		iprj    = 0

		width = 240
		xstart = 8 + boxsize/2
		ystart = 8 + boxsize/2
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
			
					proj = proj + model_gauss_noise( sigma_proj, nx, nx )
					if parm_CTF :
						proj = filt_ctf(proj, ctf)
						proj.set_attr_dict({"ctf":ctf, "ctf_applied":0})
			
					proj = proj + filt_gaussl(model_gauss_noise(sigma2_proj,nx,nx), sigma_gauss)
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

			mic += model_gauss_noise(sigma_mic,4096,4096)
			if parm_CTF :
				#apply CTF
				mic = filt_ctf(mic, ctf)
			mic += filt_gaussl(model_gauss_noise(sigma2_mic,4096,4096), sigma_gauss_mic)
	
			mic.write_image(micpref+"%1d.hdf"%(idef-3),0)
		
		drop_spider_doc("params.txt", params)
	
			
if __name__ == "__main__":
	main()

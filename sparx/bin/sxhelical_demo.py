#!/usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

# clean up the code, make documentation

import os
import global_def
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ Input Output [options]
	
	Generate three micrographs, each micrograph contains one projection of a long filament.
	Input: Reference Volume, output directory 
	Output: Three micrographs stored in output directory		
				 
		sxhelical_demo.py tmp.hdf  mic --generate_micrograph --CTF --apix=1.84	
	
	Generate noisy cylinder ini.hdf with radius 35 pixels and box size 100 by 100 by 200
	
		sxhelical_demo.py ini.hdf --generate_noisycyl --boxsize="100,100,200" --rad=35
	
	Generate rectangular 2D mask mask2d.hdf with width 60 pixels and image size 200 by 200 pixels
	
		sxhelical_demo.py mask2d.hdf --generate_mask --masksize="200,200" --maskwidth=60
	
	Apply the centering parameters to bdb:adata, normalize using average and standard deviation outside the mask, and output the new images to bdb:data
		
		sxhelical_demo.py bdb:adata bdb:data mask2d.hdf --applyparams
	
	Generate run through example script for helicon
	
		sxhelical_demo.py --generate_script --filename=run --seg_ny=180 --ptcl_dist=15 --fract=0.35
	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	# helicise the Atom coordinates
	
	# generate micrographs of helical filament
	parser.add_option("--generate_micrograph",    action="store_true",      default=False,      		  	 help="Generate three micrographs where each micrograph contains one projection of a long filament. \n Input: Reference Volume, output directory \n Output: Three micrographs containing helical filament projections stored in output directory")
	parser.add_option("--CTF",              	  action="store_true",  	default=False,   				 help="Use CTF correction")
	parser.add_option("--apix",               	  type="float",			 	default= -1,               	     help="pixel size in Angstroms")   
	parser.add_option("--rand_seed",              type="int",			    default=14567,              	 help="the seed used for generating random numbers (default 14567) for adding noise to the generated micrographs.")
	parser.add_option("--Cs",               	  type="float",			 	default= 2.0,               	 help="Microscope Cs (spherical aberation)")
	parser.add_option("--voltage",				  type="float",				default=200.0, 					 help="Microscope voltage in KV")
	parser.add_option("--ac",					  type="float",				default=10.0, 					 help="Amplitude contrast (percentage, default=10)")
	parser.add_option("--nonoise",                action="store_true",      default=False,      		  	 help="Do not add noise to the micrograph.")
	
	# generate initial volume
	parser.add_option("--generate_noisycyl",      action="store_true",      default=False,      		  	 help="Generate initial volume of noisy cylinder.")
	parser.add_option("--boxsize",                type="string",		    default="100,100,200",           help="String containing x , y, z dimensions (separated by comma) in pixels")
	parser.add_option("--rad",                    type="int",			    default=35,              	 	 help="Radius of initial volume in pixels")
	
	# generate 2D mask 
	parser.add_option("--generate_mask",          action="store_true",      default=False,      		  	 help="Generate 2D rectangular mask.")
	parser.add_option("--masksize",               type="string",		    default="200,200",               help="String containing x and y dimensions (separated by comma) in pixels")
	parser.add_option("--maskwidth",              type="int",			    default=60,              	 	 help="Width of rectangular mask")
	
	# Apply 2D alignment parameters to input stack and output new images to output stack
	parser.add_option("--applyparams",            action="store_true",      default=False,      		  	 help="Apply the centering parameters to input stack, normalize using average and standard deviation outside the mask, and output the new images to output stack")
	
	# Generate run script
	parser.add_option("--generate_script",        action="store_true",      default=False,      		  	 help="Generate script for helicon run through example")
	parser.add_option("--filename",               type="string",		    default="runhelicon",            help="Name of run script to generate")
	parser.add_option("--seg_ny",                 type="int",			    default=180,              	     help="y-dimension of segment used for refinement")
	parser.add_option("--ptcl_dist",              type="int",			    default=15,              	     help="Distance in pixels between adjacent segments windowed from same filament")
	parser.add_option("--fract",               	  type="float",			 	default=0.35,               	 help="Fraction of the volume used for applying helical symmetry.")
	
	(options, args) = parser.parse_args()
	if len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		if options.generate_script:
			generate_runscript(options.filename, options.seg_ny, options.ptcl_dist, options.fract)

		if options.generate_micrograph:
			if options.apix <= 0:
				print "Please enter pixel size."
				sys.exit()
			generate_helimic(args[0], args[1], options.apix, options.CTF, options.Cs, options.voltage, options.ac, options.nonoise, options.rand_seed)

		if options.generate_noisycyl:
			from utilities import model_cylinder, model_gauss_noise
			outvol = args[0]
			boxdims = options.boxsize.split(',')
			if len(boxdims) < 1 or len(boxdims) > 3:
				print "Enter box size as string containing x , y, z dimensions (separated by comma) in pixels. E.g.: --boxsize='100,100,200'"
				sys.exit()
			nx= int(boxdims[0])
			if len(boxdims) == 1:
				ny = nx
				nz = nx
			else:
				ny = int(boxdims[1])
				if len(boxdims) == 3:
					nz = int(boxdims[2])
					
			(model_cylinder(options.rad,nx, ny, nz)*model_gauss_noise(1.0, nx, ny, nz) ).write_image(outvol)

		if options.generate_mask:
			from utilities import model_blank, pad
			outvol = args[0]
			maskdims = options.masksize.split(',')
			if len(maskdims) < 1 or len(maskdims) > 2:
				print "Enter box size as string containing x , y dimensions (separated by comma) in pixels. E.g.: --boxsize='200,200'"
				sys.exit()
			nx= int(maskdims[0])
			if len(maskdims) == 1:
				ny = nx
			else:
				ny = int(maskdims[1])
					
			mask = pad(model_blank(options.maskwidth, ny, 1, 1.0), nx, ny, 1, 0.0)
			mask.write_image(outvol)
		
		if options.applyparams:
			from utilities    import get_im, get_params2D, set_params2D
			from fundamentals import cyclic_shift
			stack = args[0]
			newstack = args[1]
			mask = get_im(args[2])
			nima = EMUtil.get_image_count(stack)
			for im in xrange(nima):
				prj = get_im(stack,im)
				alpha, sx, sy, mirror, scale = get_params2D(prj)
				prj = cyclic_shift(prj, int(sx))
				set_params2D(prj, [0.0,0.,0.0,0,1])
				stat = Util.infomask(prj , mask, False )
				prj= (prj-stat[0])/stat[1]
				ctf_params = prj.get_attr("ctf")
				prj.set_attr('ctf_applied', 0)
				prj.write_image(newstack, im)

def generate_helimic(refvol, outdir, pixel, CTF=False, Cs=2.0,voltage = 200.0, ampcont = 10.0, nonoise = False, rand_seed=14567):
	
	from utilities	 import model_blank, model_gauss, model_gauss_noise, pad, get_im
	from random 	 import random
	from projection  import prgs, prep_vol
	from filter	     import filt_gaussl, filt_ctf
	from EMAN2 	     import EMAN2Ctf
	
	if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', "sxhelical_demo", 1)
	os.mkdir(outdir)
	seed(rand_seed)
	Util.set_randnum_seed(rand_seed)
	angles =[]
	for i in xrange(3):
		angles.append( [0.0+60.0*i, 90.0-i*5, 0.0, 0.0, 0.0] )

	nangle   = len(angles)

	volfts = get_im(refvol)
	nx = volfts.get_xsize()
	ny = volfts.get_ysize()
	nz = volfts.get_zsize()
	volfts, kbx, kby, kbz = prep_vol( volfts )
	iprj   = 0
	width  = 500
	xstart = 0
	ystart = 0

	for idef in xrange(3,6):
		mic = model_blank(2048, 2048)
		#defocus = idef*0.2
		defocus = idef*0.6     ##@ming
		if CTF :
			#ctf = EMAN2Ctf()
			#ctf.from_dict( {"defocus":defocus, "cs":Cs, "voltage":voltage, "apix":pixel, "ampcont":ampcont, "bfactor":0.0} )
			from utilities import generate_ctf
			ctf = generate_ctf([defocus,2,200,1.84,0.0,ampcont,defocus*0.2,80])   ##@ming   the range of astigmatism amplitude is between 10 percent and 22 percent. 20 percent is a good choice.
		i = idef - 4
		for k in xrange(1):
			psi  = 90 + 10*i			
 			proj = prgs(volfts, kbz, [angles[idef-3][0], angles[idef-3][1], psi, 0.0, 0.0], kbx, kby)
			proj = Util.window(proj, 320, nz)		
			mic += pad(proj, 2048, 2048, 1, 0.0, 750*i, 20*i, 0)

		if not nonoise:  mic += model_gauss_noise(30.0,2048,2048)
		if CTF :
			#apply CTF
			mic = filt_ctf(mic, ctf)

		if not nonoise:  mic += filt_gaussl(model_gauss_noise(17.5,2048,2048), 0.3)

		mic.write_image("%s/mic%1d.hdf"%(outdir, idef-3),0)

def generate_runscript(filename, seg_ny, ptcl_dst, fract):

	if ptcl_dst < 15:
		print "Distance in pixels between adjacent segments should be at least one rise!"
		sys.exit()
	
	print "Generating run script with the following parameters: \n"
	print "y-dimension of segment used for refinement: %d"%seg_ny
	print "Distance in pixels between adjacent segments: %d"%ptcl_dst
	print "Fraction of structure used for applying helical symmetry: %.2f"%fract
	
	if os.path.exists(filename):
		print "The file %s already exists. Either remove it or rename it..."%filename
		sys.exit()
		
	f = open(filename, 'w')
	f.write('#!/bin/csh\n')
	f.write('\n')
	f.write('set echo on\n')
	f.write('\n')
	f.write('#clean the previous outputs\n')
	f.write('rm *.hdf *.txt *.bck rm *.pdb\n')
	f.write('rm log*\n')
	#f.write('rm -rf outsymsearch	\n')
	f.write('rm -rf mic result_*	\n')
	f.write('rm -rf EMAN*\n')
	f.write('\n')
	f.write('# get the  previous runned results\n')
	f.write('tar -zxvf answer.tar.gz\n')
	f.write('\n')
	f.write('#generate volume from pdb\n')
	f.write('tar -zxvf 3MFP_1SU.tar.gz\n')
	f.write('\n')
	f.write('# Helicise the Atom coordinates\n')
	f.write('# Input: pdb file containing atom coordinates to helicise (3MFP_1SU.pdb)\n')
	f.write('# Output: pdb file containing helicised coordinates (rnew.pdb)\n')
	f.write('sxhelical_demo.py 3MFP_1SU.pdb rnew.pdb --heli --dp=27.6 --dphi=166.715\n')
	f.write('\n')
	f.write('#generate the density map\n')
	f.write('sxpdb2em.py rnew.pdb tmp.hdf --apix=1.84 --center=c \n')
	f.write('\n')
	f.write('# Generate three micrographs, where each micrograph contains one projection of a long filament.\n')
	f.write('# Input: Reference volume from which projections are calculated (tmp.hdf)\n')
	f.write('# Output: Output directory in which micrographs will be written (mic)\n')
	f.write('sxhelical_demo.py tmp.hdf mic --generate_micrograph --CTF --apix=1.84\n')
	f.write('\n')
	f.write('# generate initial volume for later helical refinement\n')
	f.write('# Output: A noisy cylinder written to the output file name (ini.hdf).\n')
	f.write('sxhelical_demo.py ini.hdf --generate_noisycyl --boxsize="100,100,200" --rad=35\n')
	f.write('\n')
	f.write('# Estimate defocus value for each micrograph. This can be done either by GUI or pure command-line \n')
	f.write('# 1. To estimate defocus by GUI:\n')
	f.write('cd mic\n')
	f.write('sxhelixboxer.py mic0.hdf --gui &\n')
	f.write('\n')
	f.write("# When the GUI starts up, there will be a checkbox at the bottom labelled 'CTF Estimation using CTER'. Check this box, and additional fields will appear in the GUI.\n")
	f.write("# Under 'Parameters of CTF estimation,' enter 256 for 'Window size,' 2.0 for 'Cs, 10.0 for 'Amplitude Contrast,' 200 for 'Voltage', and 1.84 for Pixel size.\n")
	f.write('# Click "Estimate CTF using CTER" button, and the estimated defocus will show up in the "Estimated defocus" box.\n')
	f.write('# Now you can either estimate defocus for the remaining micrographs one by one using the GUI, or you can estimate defocus for the remaining micrographs in batch mode using \n')
	f.write('# the parameters you entered in the GUI for mic0.hdf. \n')
	f.write('# Make sure to remove any output directories (pwrot and partres followed by underscore and \n')
	f.write('# then name of the micrograph, e.g., pwrot_mic0 and partres_mic0) generated in the previous run of CTF estimation using CTER.\n')
	f.write('rm -r pwrot_mic0;rm -r partres_mic0\n')
	f.write('sxhelixboxer.py mic*.hdf --gui &\n')
	f.write('cd ..\n')
	f.write('\n')
	f.write('# 2. To estimate defocus by command-line\n')
	f.write('cd mic\n')
	f.write('sxhelixboxer.py pwrot partres --cter --indir=. --nameroot=mic --nx=200 --Cs=2.0 --voltage=200 --kboot=16 --apix=1.84 --ac=10.0 \n')
	f.write('cd ..\n')
	f.write('\n')
	f.write('cd mic\n')
	f.write('# have to open boxer\n')
	f.write('# if want to save time, just close the window immediately and use saved coordinate\n')
	f.write('# otherwise draw the box carefully to get the coordinate \n')
	f.write('sxhelixboxer.py *.hdf --gui --helix-width=200 --ptcl-width=200\n')
	f.write('cd ..\n')
	f.write('\n')
	f.write('# if use saved coordinate, please use below commands\n')
	f.write('tar -zxvf saved_pos.tar.gz\n')
	f.write('cp saved_db.tar.gz mic/.\n')
	f.write('cp -r saved_pos/*box* mic/.\n')
	f.write('rm -rf saved_pos\n')
	f.write('cd mic\n')
	f.write('tar -zxvf saved_db.tar.gz\n')
	f.write('cd ..\n')
	f.write('\n')
	f.write('# Window segments from boxed helices\n')
	f.write('# For more information on the windowing utility, see http://sparx-em.org/sparxwiki/windowallmic\n')
	f.write('# distance between segments is 1 rise exactly\n')
	f.write('sxhelixboxer.py houtdir --window --dirid=mic --micid=mic --micsuffix=hdf --dp=27.6 --apix=1.84 --boxsize=200 --outstacknameall=bdb:hadata --hcoords_suffix=_boxes.txt --ptcl-dst=%d --rmax=64\n'%ptcl_dst)
	f.write('\n')
	f.write('# Generate 2D mask for centering EM images\n')
	f.write('# Output: 2D mask written to output file name provided (mask2d.hdf)\n')
	f.write('sxhelical_demo.py mask2d.hdf --generate_mask --masksize="200,200" --maskwidth=70\n')
	f.write('\n')
	f.write('# center all the EM images\n')
	f.write('# centering parameters will be saved in header, the input images will not be changed.\n')
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# f.write('sxheader.py bdb:hadata --params=active --one\n')
	f.write('sxheader.py bdb:hadata --params=xform.align2d --zero\n')
	f.write('mpirun -np 2 sxshiftali.py bdb:hadata mask2d.hdf --oneDx --search_rng=10 --maxit=20 --CTF --MPI\n')
	f.write('\n')
	f.write('# Apply the centering parameters stored in the header of each image in bdb:hadata to the image, \n')
	f.write('# normalize using average and standard deviation outside the mask, and save the centered and normalized image to bdb:hdata\n')
	f.write('# Input: Input stack containing 2D alignment parameters in the header (bdb:hadata)\n')
	f.write('#		 Name of 2D mask to use for normalization\n')
	f.write('# Output: Stack of images (bdb:hdata) after applying 2D alignment parameters to the images in input stack. \n')
	f.write('sxhelical_demo.py bdb:hadata bdb:hdata mask2d.hdf --applyparams\n')
	f.write('\n')
	f.write('#Exhaustive search \n')
	f.write('sxheader.py bdb:hdata --params=xform.projection --zero\n')
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# f.write('sxheader.py bdb:hdata --params=active --one\n')
	f.write('mpirun -np 3 sxhelicon.py bdb:hdata ini.hdf result_helicon --CTF --seg_ny=%d --fract=%.2f --psi_max=2.0 --delta=1.0 --maxit=7 --function=[.,nofunc,helical3c] --searchxshift=3.68 --xwobble=1.84 --ywobble=0 --dp=27.6 --dphi=166.715 --apix=1.84 --rmin=0.0 --rmax=64 --MPI\n'%(seg_ny,fract))
	f.write('\n')
	f.write('#Local search \n')
	f.write('sxheader.py bdb:hdata --params=xform.projection --import=result_helicon/parameters0007.txt\n')
	f.write('mpirun -np 3 sxlocalhelicon.py bdb:hdata result_helicon/volf007.hdf result_local --CTF --seg_ny=%d --fract=%.2f --psi_max=2.0 --delta=1.0 --maxit=11 --function=[.,nofunc,helical3c] --boundaryavg --MA --MA_WRAP=0 --xr=3.68 --txs=1.84 --an=20 --ynumber=16 --dp=27.6 --dphi=166.715 --apix=1.84 --rmin=0.0 --rmax=64 --MPI\n'%(seg_ny,fract))
	f.write('\n')
	#f.write('#Do helical symmetry search\n')
	#f.write('mpirun -np 3 sxhelicon_utils.py result_local/volf0011.hdf outsymsearch --symsearch --dp=27.6 --dphi=166.715 --apix=1.84 --fract=%.2f --rmin=0 --rmax=64.0 --datasym=datasym.txt --dp_step=0.92 --ndp=10 --dphi_step=1.0 --ndphi=10 --MPI\n'%(fract))
	
if __name__ == "__main__":
	main()

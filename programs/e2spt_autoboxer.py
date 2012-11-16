#!/usr/bin/env python

#
# Author: Jesus Galaz, 10/20/2012; last update 11/15/2012
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
#
#

import os
from EMAN2 import *
from time import time
		 
import matplotlib
matplotlib.use('Agg')
		 
import matplotlib.pyplot as plt
import sys
import numpy		 

	 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """WARNING:  **PRELIMINARY** program, still heavily under development.
	
	Autoboxes globular particles from tomograms."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--pad", action='store_true', help="""Provide this if the particles in the --particlestack used to create a template are in a tight box. 
															Also, provide it if you will be doing rotational searches if providing a --template in a tight box.""",default=False)
	parser.add_argument("--rotationalsearch", action='store_true', help="""At each translation position, vary euler angles as well when searching for particles.""",default=False)
	
	parser.add_argument("--tomogram", type=int, help="Name of the tomogram.",default='')
	parser.add_argument("--goldstack", type=int, help="Name of the stack containing a few gold particles picked from the tomogram.",default='')
	parser.add_argument("--particlestack", type=int, help="""Name of the stack containing a few sample particles picked from the tomogram, used to create an initial template.
															with which to search for particles throughout the tomogram.""",default='')
	
	parser.add_argument("--template", type=int, help="Name of the file containing the template to search for particles throughout the tomogram.",default='')

	parser.add_argument("--backgroundstack", type=int, help="Name of the stack containing a few boxes picked from regions of the tomogram where there where no particles, no gold, and no carbon.",default='')
	parser.add_argument("--carbonstack", type=int, help="Name of the stack containing a few boxes picked from the grid hole (or carbon).",default='')

	parser.add_argument("--output", type=int, help="Name to output the auto-boxed particles.",default='')
	parser.add_argument("--shrink", type=int, help="Integer factor by which the tomogram will be shrunk.", default=1)
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the tomogram", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the tomogram", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the tomogram.", default=None)
	parser.add_argument("--mask",action='store_true',help="""If provided, a cylindrical mask will be created to mask out the carbon and keep only the grid hole.
															--gridradius and --gridoffest must be specified.""")
	
	parser.add_argument("--gridradius", type=int, help="Radius of the grid in pixels. Supply this parameter only if also supplying --mask.",default='')
	parser.add_argument("--gridoffset", tpe='str', help="""x,y coordinates for the center of the grid hole in the center slice of the tomogram (or if you generated a 2D projection of the tomogram. 
														The left bottom corner would be 0,0. Supply this parameter only if also supplying --mask and the grid hole is not centered in the tomogram.""")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptautobox'; 
															for example, sptautobox_01 will be the directory by default if 'sptautobox' already exists.""")
	
	parser.add_argument("--yshort",action="store_true",default=False,help="This means you have provided a --tomogram in which 'y' is the short axis.")
	
	(options, args) = parser.parse_args()
	
	rootpath = os.getcwd()
		
	logger = E2init(sys.argv, options.ppid)
	
	'''
	Check that the output and tomogam formats are sane, to prevent crashes later.
	'''
	if '.hdf' not in options.output and '.mrc' not in options.output:
		print "ERROR: The output reference must be written to either an mrc or an hdf file."
		sys.exit() 
	
	if '.hdf' not in options.tomogram and '.mrc' not in options.tomogram and '.rec' not in options.tomogram:
			print "ERROR: The tomogram must be in .mrc, .rec (which is also just a .mrc file) or .hdf format."
			
	'''
	Make the directory where to create the database where the results will be stored
	'''
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)

	if not options.path: 
		#options.path="bdb:"+numbered_path("sptavsa",True)
		options.path = "sptCudaTest_01"
	
	files=os.listdir(os.getcwd())
	print "right before while loop"
	while options.path in files:
		print "in while loop, options.path is", options.path
		#path = options.path
		if '_' not in options.path:
			print "I will add the number"
			options.path = options.path + '_00'
		else:
			jobtag=''
			components=options.path.split('_')
			if components[-1].isdigit():
				components[-1] = str(int(components[-1])+1).zfill(2)
			else:
				components.append('00')
						
			options.path = '_'.join(components)
			#options.path = path
			print "The new options.path is", options.path

	if options.path not in files:
		
		print "I will make the path", options.path
		os.system('mkdir ' + options.path)	
	
	options.path = rootpath + '/' + options.path
	
	'''
	Shrink tomogram (if required), then load it.
	'''
	if options.shrink:
		binnedname = options.path + '/' + options.tomogram.split('/')[-1].replace('.','_bin' str(options.shrink) + '.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + binnedname + ' --process=math.meanshrink:n=' + options.shrink)
		options.tomogram = binnedname
	else:
		tomogramfile = options.path + '/' + options.tomogram		
		
	'''
	Apply any specified filters (preprocess, lowpass and/or highpass)
	'''
	if options.preprocess != None:
		print "Patience. Applying this processor to tomogram:", options.preprocess
		tomogramfile = tomogramfile.replace('.','_pr.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.preprocess)
		options.tomogram = tomogramfile
		
	if options.lowpass != None:
		print "Patience. Applying lowpass filter to tomogram:", options.lowpass
		tomogramfile = tomogramfile.replace('.','_lp.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.lowpass)
		options.tomogram = tomogramfile

	if options.highpass != None:
		print "Patience. Applying highpass filter to tomogram:", options.highpass
		tomogramfile = tomogramfile.replace('.','_hp.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.highpass)
		options.tomogram = tomogramfile
	
	'''
	Apply any masks if there are big carbon regions beyond the grid-hole and masking is specified.
	'''
	xo=0
	yo=0		
	if options.mask and options.gridradius:
		tomohdr = EMData(tomogramfile,0,True)
		height = tomohdr['nz']
		#mask=EMData(options.gridradius*2,options.gridradius*2,height])
		mask=EMData(tomohdr['nx'],tomohdr['ny'],height])

		if options.yshort:
				height = tomohdr['ny']
				#mask=EMData(options.gridradius*2,height,options.gridradius*2)			
				mask=EMData(tomohdr['nx'],tomohdr['nz'],height])
		mask.to_one()
		mask.process_inplace('testimage.cylinder',{'radius':options.gridradius,'height':height})
		
		if options.gridoffset:
			xo=options.gridoffset.split(',')[0]
			yo=options.gridoffset.split(',')[-1]
			mask.translate(xo,yo)	
		
		tomo=EMData(tomogramfile)
		tomo.mult(mask)
		
		tomogramfile=tomogramfile.replace('.','_msk.')
				
		tomo.write_image(tomogramfile,0)
		options.tomogram=tomogramfile

	'''
	If a template is provided, check that it is sane. If not, fix it.	
	'''
	if options.template:
		if '.hdf' not in options.template and '.pdb' not in options.template and '.mrc' not in options.template:
			print "ERROR: The format of the template to use must be .pdb, .mrc or .hdf"
			sys.exit()
		else:
			check=0
			if '.pdb' in options.template:
				pdbmodel = options.template
				#os.system('cp ' + pdbmodel + ' ' + options.path)
				pdbmodel = pdbmodel.split('/')[-1]
				mrcmodel = pdbmodel.replace('.pdb','.mrc')
				os.system('e2pdb2mrc.py ' + pdbmodel + ' ' + mrcmodel)
				options.template = mrcmodel
				check=1
				if options.verbose:
					print "I've converted the .pdf template to .mrc"
		
			if '.mrc' in options.template and '.hdf' in options.output:
				mrcmodel = options.template
				if check==0:
					#os.system('cp ' + mrcmodel + ' ' + options.path)
					mrcmodel = mrcmodel.split('/')[-1]
				hdfmodel = mrcmodel.replace('.mrc','.hdf')
				os.system('e2proc3d.py ' + mrcmodel + ' ' + hdfmodel + ' && rm ' + mrcmodel)
				options.template = hdfmodel
				check=1
				if options.verbose:
					print "I've converted the .mrc template to .hdf"
		
			template = EMData(options.template,0)
			
			if options.verbose:
				print "I have loaded the template."
	
			'''
			Make all sides of the template box equally sized and even, if they're not
			'''
			if template['nx'] != template['ny'] or template['nx'] != template['nz'] or template['ny'] != template['nz']:	
				x0 = template['nx']
				y0 = template['ny']
				z0 = template['nz']
	
				side = max(x0,y0,z0)
				if side % 2:
					side += 1
		
				R = Region((x0 - side)/2, (y0 - side)/2, (z0 - side)/2, side, side, side)
				template.clip_inplace(R)
		
				if options.verbose:
					print "I have clipped the template to have a cubical and even-sided box size."
				
				options.template = options.template.replace('.','_fixed.')
				template.write_image(options.path + '/' + options.template,0)
		
	elif options.particlestack:
		options = generateref(options)
			
	'''
	Scan the tomogram to extract subvolumes to compare against the template, only if they are inside the grid-hole (in this case, in side the mask if there is any)
	'''
	x = tomohdr['nx']
	y = tomohdr['ny']
	z = tomohdr['nz']
	if options.yshort:
		aux=z
		z=y
		y=aux
	
	'''
	Make a pseudo-spherical or shell-like template with the radius of gyration of the actual template
	for purposes of translational search
	'''
	
	transtemplate = EMData(options.template,0)
	transtemplate = transtemplate.rotavg_i()
	transtemplatename = options.template.replace('.','_sphavg.')
	boxsize = z
	transtemplate.process_inplace('xform.scale',{'scale':1,'clip',boxsize})
	transtemplate.write_image(transtemplatename,0)
	k=boxsize/2
	for i in xrange(boxsize/2, x - boxsize/2 + 1):
		if not i % boxsize/2:
			for j in xrange(boxsize/2 , y - boxsize/2 + 1):
				if not j % boxsize/2:
					if options.mask:
						if (x - xo)*(x - xo) + (y - yo)*(y - yo) < options.gridradius * options.gridradius
							scanposition(options,transtemplate,i,j,k)
						else:
							scanposition(options,transtemplate,i,j,k)
	return()
	
	
def generateref(options):
	ptclhdr = EMData(options.particlestack,0,True)
	boxsize = ptclhdr['nx']		
	alicmd = "e2spt_classaverage.py --path=" + options.path + "/generatereference" + " --input=" + options.input + " --output=" + options.input.replace('.hdf','_avg.hdf') + " --npeakstorefine=10 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_translate_3d:search=" +str( int(boxsize/4) ) + ":delta=15:dphi=15:verbose=0 --parallel=thread:7 --ralign=refine_3d_grid:delta=5:range=15:search=2 --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --savesteps --saveali --normproc=normalize.mask"
	os.system(cmd)
	options.template = options.path + "/generatereference/" + options.input.replace('.hdf','_avg.hdf')
	return(options)
	

def scanposition(options,template,x,y,z):
	boxsize = template['nx']
	r = Region( (2*x-boxsize)/2,(2*y-boxsize)/2, (2*z-boxsize)/2, boxsize, boxsize, boxsize )
	subox = EMData()
	subox.read_image(options.tomogram,0,False,r)
	
	ccf = template.calc_ccf(subox)
	ccf.process_inplace("xform.phaseorigin.tocorner") 
	ccf.process_inplace('normalize')

	#box = ccf.get_zsize()
	#r =  Region((box/2) - int(parameters['searchx']), (box/2) - int(parameters['searchy']), (box/2) - int(parameters['searchz']), 2 * int(parameters['searchx']) + 1, 2 * int(parameters['searchy']) + 1, 2 * int(parameters['searchz']) + 1) 
	#sub_ccf = ccf.get_clip(r)

	locmax = ccf.calc_max_location()
	xmax = locmax[0]
	ymax = locmax[1]
	zmax = locmax[2]
	max = ccf['maximum']
	
	
	
	
	return(e)
	
	
			 

if '__main__' == __name__:
	main()
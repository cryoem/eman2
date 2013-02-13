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
	
	parser.add_argument("--pad", action='store_true', help="""Provide this if the particles in the --particlestack used to create a template, or the template supplied 
															through --template are in a tight box. The size""",default=False)
	parser.add_argument("--rotationalsearch", action='store_true', help="""At each translation position, vary euler angles as well when searching for particles.""",default=False)
	
	parser.add_argument("--pruneccc", action='store_true', help="""Pruned based on ccc mean and sigma.""",default=False)
	parser.add_argument("--invert", action='store_true', help="""Multiply tomogram subsections my -1 to invert the contrast BEFORE looking for particles.""",default=False)
	parser.add_argument("--prunerepeated", action='store_true', help="""Multiply tomogram subsections my -1 to invert the contrast BEFORE looking for particles.""",default=False)


	
	parser.add_argument("--tomogram", type=str, help="Name of the tomogram.",default='')
	parser.add_argument("--goldstack", type=str, help="Name of the stack containing a few gold particles picked from the tomogram.",default=None)
	parser.add_argument("--ptclstack", type=str, help="""Name of the stack containing a few sample particles picked from the tomogram, used to create an initial template.
															with which to search for particles throughout the tomogram.""",default=None)
	
	parser.add_argument("--template", type=str, help="Name of the file containing the template to search for particles throughout the tomogram.",default='')

	parser.add_argument("--backgroundstack", type=str, help="""Name of the stack containing a few boxes picked from regions of the tomogram where there where no particles, 
															no gold, and no carbon.""",default=None)
	parser.add_argument("--carbonstack", type=str, help="Name of the stack containing a few boxes picked from the grid hole (or carbon).",default=None)

	parser.add_argument("--output", type=str, help="Name to output the auto-boxed particles.",default='')
	parser.add_argument("--outputboxsize", type=int, help="Size of the box to put the extracted particles in, and amount by which the subregions will overlap, when searching for particles in the tomogram.", default=0)

	parser.add_argument("--shrink", type=int, help="Integer factor by which the tomogram will be shrunk.", default=0)
	
	parser.add_argument("--subsettrans", type=int, help="Subset of particles to keep/consider after translational alignment.", default=0)
	parser.add_argument("--concentrationfactor", type=int, help="""Determines how many particles will be pre-picked as putative particles. For example, if
																if the tomogram is broken up into subregions of volume V to look for particles in each
																and --concentrationfactor=1, then, the number of best-correlating subvolumes from the subregion
																that will be initially selected as particles will be n=V/(pv*C) where 'pv' is the volume of one particle
																calculated based on the outputboxsize, or the template's boxsize; 'C' is the concentration factor;
																therefore, the largest it is, the fewer particles that will be initially picked.""", default=0)

	parser.add_argument("--apix", type=float, help="The actual apix of the tomogram if for some reason it is wrong on the header.", default=0.0)
	
	parser.add_argument("--ptclradius", type=int, help="The estimated radius of the particle in pixels.", default=0)

	parser.add_argument("--cshrink", type=int, help="""If the tomogram was PREVIOUSLY shrunk, --cshrink is the factor by which the tomogram supplied through --tomogram was shrunk with respect to 
														the raw (unshrunk) tomgoram. This CAN work in conjuction with --shrink, so be careful. If both parameters are specified,
														the coordinates found by the autoboxer will be multiplied by BOTH factors.""", default=0)

	#parser.add_argument("--boxsize", type=int, help="Boxsize to resize the template, either provided through --template or computed from --particlestack", default=1)

	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the tomogram", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the tomogram", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the tomogram.", default=None)
	parser.add_argument("--mask",action='store_true',help="""If provided, a cylindrical mask will be created to mask out the carbon and keep only the grid hole.
															--gridradius and --gridoffest must be specified.""", default=False)
	
	parser.add_argument("--gridradius", type=int, help="Radius of the grid in pixels. Supply this parameter only if also supplying --mask.",default=0)
	parser.add_argument("--gridoffset", type=str, help="""x,y coordinates for the center of the grid hole in the center slice of the tomogram (or if you generated a 2D projection of the tomogram. 
														The left bottom corner would be 0,0. Supply this parameter only if also supplying 
														--mask and the grid hole is not centered in the tomogram.""", default='')
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptautobox'; 
															for example, sptautobox_01 will be the directory by default if 'sptautobox' already exists.""")
	
	#parser.add_argument("--yshort",action="store_true",default=False,help="This means you have provided a --tomogram in which 'y' is the short axis.")
	parser.add_argument("--test",action="store_true",default=False,help="This means you have provided a --tomogram in which 'y' is the short axis.")
	parser.add_argument("--templatethreshold",type=float,default=0.0,help="""A binary threshold will be applied to the template which will zero out all the densities below the supplied value, 
												and will make the densities above the supplied value equal to one.""")

	(options, args) = parser.parse_args()
	
	rootpath = os.getcwd()
		
	logger = E2init(sys.argv, options.ppid)
	
	#if not options.template and not options.ptclstack:
	#	print "TERMINATING: You must provide either a template, through --template, or a stack of particles to build one, through --ptclstack."
	#	sys.exit()
	
	print "THE particle radius at reading is!", options.ptclradius
	'''
	Check that the output and tomogam formats are sane, to prevent crashes later.
	'''
	if '.hdf' not in options.output and '.mrc' not in options.output:
		print "ERROR: The output stack must be written to either an mrc or an hdf file."
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
		options.path = "sptAutoBox_01"
	
	files=os.listdir(os.getcwd())
	while options.path in files:
		#path = options.path
		if '_' not in options.path:
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

	if options.path not in files:
		
		os.system('mkdir ' + options.path)	
	
	options.path = rootpath + '/' + options.path

	tomo = EMData(options.tomogram,0)	
	
	if options.apix:
		tomo['apix_x'] = options.apix
		tomo['apix_y'] = options.apix
		tomo['apix_z'] = options.apix

		tomogramfile = options.path + '/' + options.tomogram.split('/')[-1].replace('.','_edED.')
		options.tomogram = tomogramfile

		tomo['origin_x'] = 0
		tomo['origin_y'] = 0
		tomo['origin_z'] = 0

		tomo.write_image(tomogramfile,0)
		tomo = EMData(options.tomogram,0)
	
	'''
	Shrink tomogram (if required), then load it.
	'''
	tomogramfile = options.tomogram
	if options.shrink > 1:
		outputname=''
		if "edED" in options.tomogram:
			outputname=options.tomogram
		else:
			outputname = options.path + '/' + options.tomogram.split('/')[-1].replace('.','_edED.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + outputname + ' --process=math.meanshrink:n=' + str(options.shrink))
		options.tomogram = outputname
		tomogramfile = options.tomogram
		tomo = EMData(options.tomogram,0)

	tomox = tomo['nx']
	tomoy = tomo['ny']
	tomoz = tomo['nz']
	
	
	tomogramapix = tomo['apix_x']
	if options.apix:
		tomogramapix = options.apix
	
	yshort=False
	
	print "Before preprocessing, the tomogram is located at", tomogramfile	
	'''
	Apply any specified filters (preprocess, lowpass and/or highpass)
	'''
	if options.preprocess:
		print "Patience. Applying this processor to tomogram:", options.preprocess
		tomogramfile = options.path + '/' + tomogramfile.replace('.','_pr.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.preprocess)
		options.tomogram = tomogramfile
		
	if options.lowpass:
		print "Patience. Applying lowpass filter to tomogram:", options.lowpass
		print "Whose type is", type(options.lowpass)
		
		if options.path in tomogramfile:
			tomogramfile = tomogramfile.replace('.','_lp.')
		else:
			tomogramfile = options.path + '/' + tomogramfile.replace('.','_lp.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.lowpass)
		options.tomogram = tomogramfile

	if options.highpass:
		print "Patience. Applying highpass filter to tomogram:", options.highpass
		
		if options.path in tomogramfile:
			tomogramfile = tomogramfile.replace('.','_hp.')
		else:
			tomogramfile = options.path + '/' + tomogramfile.replace('.','_hp.')
		
		os.system('e2proc3d.py ' + options.tomogram + ' ' + tomogramfile + ' --process=' + options.highpass)
		options.tomogram = tomogramfile
	
	'''
	Apply any masks if there are big carbon regions beyond the grid-hole and masking is specified.
	'''
	xo=0
	yo=0
	
	if tomo['ny'] < tomo['nz']:
		yshort = True
			
	if options.mask and options.gridradius:
		height = tomoz
		
		'''
		testimage.cylinder only works on cubical boxes. Therefore, make a cubical box to generate the mask.
		'''
		cubem=max(tomox,tomoy,tomoz)
		mask=EMData(cubem,cubem,cubem)

		if yshort:
			height = tomoy
		
		mask.to_one()
		mask.process_inplace('testimage.cylinder',{'radius':options.gridradius,'height':height})
		
		if options.gridoffset:
			xo=options.gridoffset.split(',')[0]
			yo=options.gridoffset.split(',')[-1]
			mask.translate(xo,yo)	
		
		if yshort:
			print "I will rotate the mask"
			yo*=-1
			mask.rotate(0,90,0)
		
		'''
		Then, clip mask to actual tomogram size
		'''
		r=Region(-tomox/2,-tomoy/2,-tomoz/2,tomox,tomoy,tomoz)
		mask=mask.get_clip(r)
		
		print "The dimensions of the mask are", mask['nx'],mask['ny'],mask['nz']
		print "The dimensions of the tomogram are", tomox,tomoy,tomoz
	
		tomo.mult(mask)
		
		tomogramfile=options.path + '/' + tomogramfile.split('/')[-1].replace('.','_msk.')
		tomo.write_image(tomogramfile,0)
		options.tomogram=tomogramfile

	'''
	If a template is provided, check that it is sane. If not, fix it.	
	'''
	
	expandedtemplateboxsize=0
	if options.template:
	
		print "You have provided the following template", options.template
		
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
			
			print "I have loaded the template and this is its type:", type(template)
			if options.verbose:
				print "I have loaded the template."
	
			'''
			Make all sides of the template box equally sized and even, if they're not
			'''
			x0 = template['nx']
			y0 = template['ny']
			z0 = template['nz']
			if template['nx'] != template['ny'] or template['nx'] != template['nz'] or template['ny'] != template['nz']:
				side = max(x0,y0,z0)
				
				if options.boxsize:
					side = options.boxsize
				
				elif options.pad:
					side *= 2
				
				if side % 2:
					side += 1
				
				R = Region((x0 - side)/2, (y0 - side)/2, (z0 - side)/2, side, side, side)
				template.clip_inplace(R)
		
				if options.verbose:
					print "I have clipped the template to have a cubical and even-sided box size."
				
				options.template = options.template.replace('.','_fixed.')
				template.write_image(options.path + '/' + options.template,0)
				expandedtemplateboxsize = side
			
			else:
				expandedtemplateboxsize = x0
		
	elif options.ptclstack:
		#ret = generateref(options)
		#options = ret[0]
		#expandedtemplateboxsize = ret[1]
		
		print "\n\n\n\n YOU HAVE PROVIDED A PTCLSTACK!!!!"
		print "\n\n\n"
		ptclhdr = EMData(options.ptclstack,0,True)
		box = ptclhdr['nx']
		outputboxsize = box
		if options.outputboxsize:
			outputboxsize == options.outputboxsize
		
		#genrefpath = options.path + "/genref"
		#os.system('mkdir ' + genrefpath)
		cmd = "cd " + options.path + " && e2spt_classaverage.py --path=genref --input=../" + options.ptclstack + " --output=" + options.ptclstack.replace('.hdf','_sphavg.hdf') + " --npeakstorefine=10 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_translate_3d:search=" +str( int(box/4) ) + ":delta=15:dphi=15:verbose=0 --parallel=thread:7 --ralign=refine_3d_grid:delta=5:range=15:search=2 --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --savesteps --saveali --normproc=normalize.mask --nocenterofmass"
		if options.verbose:
			print "I will generate a template from --particlestack by executing the following command:", cmd
	
		os.system(cmd)
		options.template = options.path + "/genref/" + options.ptclstack.replace('.hdf','_sphavg.hdf')
		
		
	elif options.outputboxsize:
		
		if options.ptclradius:
			outputboxsize = options.outputboxsize
			a = EMData(outputboxsize,outputboxsize,outputboxsize)
			a.to_one()
			a.process_inplace('mask.sharp',{'outer_radius':options.ptclradius})
			a['apix_x'] = tomogramapix
			a['apix_y'] = tomogramapix
			a['apix_z'] = tomogramapix

			a.write_image(options.path + '/sph_template.hdf',0)
			#if not options.invert:
			#	a.mult(-1)
				
			options.template = options.path + '/sph_template.hdf'
		
		else:
			print "ERROR: You didn't provide --template or --ptclstack. In the absence of these, you must provide --outputboxsize AND --ptclradius."
			sys.exit()
	else:
		print "ERROR: You didn't provide --template or --ptclstack. In the absence of these, you must provide --outputboxsize."
		sys.exit()
			
	'''
	Scan the tomogram to extract subvolumes to compare against the template
	'''
	x = tomox
	y = tomoy
	z = tomoz
	
	print "!!!!Original tomo dimensions are", x,y,z
	
	
	print "The template is", options.template
	
	transtemplate = EMData(options.template,0)
	transtemplateapix = transtemplate['apix_x']
	print "With an apix of", transtemplateapix
	
	'''
	Calculate the scale between the reference model and the data, round to the nearest integer to be able to use math.meanshrink (which preserves features),
	then calculate the scale factor again and use xform.scale, to get the exact scaling, even if features might become a bit distorted by xform.scale
	'''	
	meanshrinkfactor = tomogramapix/transtemplateapix
	meanshrinkfactor_int = int(round(meanshrinkfactor))
	#print "The template's apix is", transtemplateapix
	#print "And the tomgoram's apix is", tomogramapix
	#print "Therefore, the meanshrink factor is", meanshrinkfactor
	#print "Which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int
	if meanshrinkfactor_int > 1:
		if options.verbose:
	
			print "About to shrink"
			#print "The type of template is", type(transtemplate)
			print "\n\n\n BBBBBBBB\nANd its boxsize BEFORE shrinking was\nBBBBBBBB\n\n", transtemplate['nx']
		
		transtemplate.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		
		expandedtemplateboxsize = transtemplate['nx']
		print "\n\n\nCCCCCCCC\nWhereas AFTER shrinking it is\nCCCCCC\n\n", transtemplate['nx']

		if options.verbose:
			print "The template was shrunk, to a first approximation."
		
		transtemplateapix = transtemplate['apix_x']
		
	scalefactor = transtemplateapix/tomogramapix
	expandedtemplateboxsize = round(transtemplate['nx']*scalefactor)
	

	
	#if options.verbose:
	print "The finer scale factor to apply is", scalefactor
	
	if float(scalefactor) != 1.0:
		transtemplate.process_inplace("xform.scale",{"scale":scalefactor,"clip":expandedtemplateboxsize})
		expandedtemplateboxsize = transtemplate['nx']
	#transtemplate.process_inplace("xform.scale",{"scale":1,"clip":transtemplatebox})
	
	print "\n\n\n AAAAAAAAAAAA\n after all necessary APIX MATCHING, the expandedboxsize of the template is", expandedtemplateboxsize
	print "\n\n\n AAAAAAAAAAAAA\n"
	
	if options.verbose:
		print "The template has now been precisely shrunk."

	'''
	Make a pseudo-spherical or shell-like template with the radius of gyration of the actual template
	for purposes of translational search
	'''
	outputboxsize = expandedtemplateboxsize
	if options.outputboxsize:
		outputboxsize = options.outputboxsize
	
	ptclradius = round( ( transtemplate['nx'] / 2.0 ) )
	
	if options.outputboxsize:
		ptclradius = round( outputboxsize / 2.0 )
		if options.shrink > 1:
			ptclradius = round( ptclradius/ options.shrink )
	
	if options.ptclradius:
		ptclradius = options.ptclradius	
		print "\n\nThe particle radis BEFORE any shrinking is", options.ptclradius	
		if options.shrink > 1:
			ptclradius = round( ptclradius/ options.shrink )
			options.ptclradius = ptclradius
			print "ptclradius being shrunk!!", options.ptclradius
	
	print "\n\n\nRRRRRRRRR\nThe particle radius in pixels is %f \nRRRRRRRRn\n\n" %(ptclradius)
	
	transtemplatename = options.template
	if options.ptclstack:
		transtemplate.process_inplace('mask.sharp',{'outer_radius':ptclradius})
		transtemplate = transtemplate.rotavg_i()
	
		if options.templatethreshold:
			transtemplate.process_inplace('threshold.binary',{'value':options.templatethreshold})
	
		transtemplatename = options.template.split('/')[-1].replace('.','_sphavg.')
		if options.path not in transtemplatename:
			transtemplatename = options.path + '/' + transtemplatename
		transtemplate.write_image(transtemplatename,0)
	
	
	
	regionboxsize = z
	if yshort:
		regionboxsize = y
	print "The volume of the tomogram is", x*y*z 
	print "The volume if the subregions is", regionboxsize*regionboxsize*regionboxsize
	subcubes = int( (x*y*z) / (regionboxsize*regionboxsize*regionboxsize) )

	if options.verbose:
		print "Therefore, for boxing purposes, the tomogram will be divided into at least these many sections", subcubes
	
	transtemplate.process_inplace('xform.scale',{'scale':1,'clip':regionboxsize})
	#mskrad = ptclboxsiz/2 - 1
	#transtemplate.process_inplace('mask.sharp',{'outer.radius':mskrad})
	transtemplate.write_image(transtemplatename,0)
	zc=regionboxsize/2
	
	expandedtemplateboxsize = regionboxsize
	
	'''
	The tomogram ought to be divided in sub-blocks, which will be cubes with side length equal to the 
	tomogram's thickness (z). When scanning each of these sub-blocks for particles, they need to overlap
	by an amount equal to the box length of the particles (outputboxsize) to ensure picking particles at the edges of each
	sub-block without 'breaking' the particle. 
	The first row and column of sub-blocks need to be centered at boxsize/2.
	All subsequent ones need to be centered at multiples of expandedtemplateboxsize/2 - outputboxsize.
	'''
	
	print "The template has been cipped such that it has a size equal to the short size of the tomogram", expandedtemplateboxsize
	
	if options.verbose:
		print "I will look for subvolumes throughout the tomogram, which has these x, y, z dimensions", x, y, z
		print "Where the boxsize of SUBREGIONS to examine is", regionboxsize
		print "And the outputboxsize of (UNSHRUNK) ptcls to extract is", outputboxsize
	
	if options.shrink > 1:
		expandedtemplateboxsize=expandedtemplateboxsize/2
		print "Because you chose to shrink, the expandedtemplateboxsize for scanning purposes is actually", expandedtemplateboxsize
		print "And therefore, the outputbox (which determines the effective boxisze of the particles in the tomogram DURING scanning, is"
		outputboxsize = outputboxsize / 2
		print outputbox
		
	lines=[]
	
	i=0
	xc = regionboxsize/2
	otherlong = y
	short =z
	if yshort:
		otherlong = z
		short = y
	
	data=[]
	coordset=set()
	ret=None
	
	sbn=0
	zi=0
	
	#print "expandedtemplateboxsize size is", expandedtemplateboxsize
	#print "Therefore its half is", expandedtemplateboxsize/2
	#print "\n\n\n\n@@@@@@@The subregions o the tomogram to scan are THESE!!!!!\n\n"
	boxes=[]
	centers=[]
	xi=0
	xc = regionboxsize/2
	factor=1
	
	#print "\n\n\nFFFFFFFFFFF\n Becuase outputboxsize is", outputboxsize
	#print "And regionboxsize is", regionboxsize
	
	if regionboxsize == outputboxsize:
		factor = regionboxsize / 2
		print "Regionboxsize/2"
	
	elif outputboxsize < regionboxsize/2:
		factor = regionboxsize / 2 + (regionboxsize/2 - outputboxsize)
		print "regionboxsize - outputboxsize"
		
	elif outputboxsize >= regionboxsize/2 and outputboxsize < regionboxsize:
		factor = regionboxsize/2 + (regionboxsize/2 - outputboxsize/2)
		print "outputboxsize/2"
	
	print "\n\n\nFFFFFFF\nFactor is",factor
	print "\nFFFFFFFFFFFFF\n\n\n\n"
	
	coordsset = set()
	coeffset = set()
	rmsdset = set()
	
	count=0
	finalx=0
	while xi < x and xc <= x -regionboxsize/2 and finalx < 2:
		yi=0
		yc = regionboxsize/2
		#print "xi is", xi
		finaly=0		
		while yi < y and yc <= y -regionboxsize/2 and finaly < 2:
			zi=0
			zc = regionboxsize/2
		#	print "yi is", xi
			while zi < z and zc <= z -regionboxsize/2:
				#print "zi is", xi
				#print "Therefore the initial box coords is", xi,yi,zi
				box=[xi,yi,zi, xi+regionboxsize , yi+regionboxsize , zi+regionboxsize]
				center=[xc,yc,zc]
				boxes.append(box)
				centers.append(center)
				
				if options.mask and options.gridradius:
					if (xc - xo)*(xc - xo) + (yc - yo)*(yc - yo) < options.gridradius * options.gridradius:
						ret = scanposition(options,transtemplate,outputboxsize,yshort,xi,yi,zi,xc,yc,zc,sbn)
				else:
					ret = scanposition(options,transtemplate,outputboxsize,yshort,xi,yi,zi,xc,yc,zc,sbn)
				
				if ret:
					#print "These many particles have been returned", len(ret)
					asdfg=0
					for r in ret:
						#print "I will examine particle", asdfg, r
						#print "This info has been returned", r
						newcoeff = r[0]
						newcoefftuple = tuple( [newcoeff] )
						coords=tuple( [ r[1],r[2],r[3] ] )
						#oldcoeff=0
						
						if coords in coordset:
							#print "The particle was already in coordset see", coordset							
							for kkk in range(len(data)):
								if list(coords) in data[kkk]:
									oldcoeff = data[kkk][0]
									if oldcoeff < newcoeff:
										data[kkk][0]=newcoeff
										#print "The particlce was already in data, but the newcoeff is higher than the old", oldcoeff, newcoeff
										
						elif coords not in coordset:
								
							coordset.add(coords)
							coeffset.add(newcoefftuple)
							data.append( [newcoeff,coords] )
							#print "I have appended a new particle!", [newcoeff,coords]
							count+=1
						asdfg+=1
							
				#print "And I have indeed appended such box", box
				sbn+=1
				#zi = zi + (factor + regionboxsize - outputboxsize)
				#zc = zc + (factor + regionboxsize - outputboxsize)	
				zi += factor
				zc += factor
			#yi = yi + (factor + regionboxsize - outputboxsize)
			#yc = yc + ( factor + regionboxsize - outputboxsize)
			yi += factor
			yc += factor
			if yi > y or yc > (y - regionboxsize/2 ) and y % factor:
				#print "I am at a special last box in Y! Therefore the old yi and yc", yi, yc 
				yi = y - regionboxsize
				yc = y - regionboxsize/2
				#print "will be changed for", yi, yc
				finaly+=1
			
			if options.test:
				yi*=2
				yc*=2

		#xi = xi + ( factor + regionboxsize - outputboxsize )
		#xc = xc + ( factor + regionboxsize - outputboxsize )
		xi += factor
		xc += factor
		if xi > x or xc > y - regionboxsize/2 and x % factor:
			#print "I am at a special last box in X! Therefore the old xi and xc", xi, xc 
			xi = x - regionboxsize
			xc = x - regionboxsize/2
			#print "will be changed for", xi, xc
			finalx+=1

		if options.test:
			xi*=2
			xc*=2
	print "\n\nCCCCCCCCCCCC\nThe total number of appended particles was %d\nCCCCCCCCCC\n\n" %(count)
	#for m in range(len(boxes)):
	#	print 'box', boxes[m]
	#	print '\t center', centers[m]
	
	data.sort()
	data.reverse()
	
	#for i in data:
	#	print "This is in data", i

	print "the len of data WAS", len(data)
	

	'''
	bg_stack_mean=0
	bg_stack_sigma=0
	if options.backgroundstack:
		ret=meancalc(options.backgroundstack)
		bg_stack_mean=ret[0]
		bg_stack_sigma=ret[1]
	
	gold_stack_mean=0
	gold_stack_sigma=0
	if options.goldstack:
		gold=meancalc(options.goldstack)
		gold_stack_mean=ret[0]
		gold_stack_sigma=ret[1]
	
	carb_stack_mean=0
	carb_stack_sigma=0
	if options.cabonstack:
		ret=meancalc(options.carbonstack)
		carb_stack_mean=ret[0]
		carb_stack_sigma=ret[1]
	'''
	
	lendata1=len(data)
	if options.backgroundstack:
		print "BG Stack to send is", options.backgroundstack
		ret = meancalc(options,data,options.backgroundstack,tag='background')
		data=ret[0]
		lendata2=len(data)
		pruned= lendata1-lendata2
		print "These many particles were pruned using the background stack",pruned
	
	lendata1=len(data)
	if options.carbonstack:
		print "CARBON Stack to send is", options.carbonstack
		ret = meancalc(options,data,options.carbonstack,tag='carbon')
		data=ret[0]
		lendata2=len(data)
		pruned= lendata1-lendata2
		print "These many particles were pruned using the carbon stack",pruned
	
	lendata1=len(data)
	if options.goldstack:
		print "GOLD Stack to send is", options.goldstack
		ret = meancalc(options,data,options.goldstack,tag='gold')
		data=ret[0]
		lendata2=len(data)
		pruned= lendata1-lendata2
		print "These many particles were pruned using the gold stack",pruned
	
	
	
	
	
	
	
	'''
	Prune particles based on correlation. Peaks in correlation should be representative of particles, and only N of those
	peaks are kept for non-overlapping boxes when scanning sub-regions of the tomogram, based on --concentrationfactor .
	These can be further pruned by eliminating particles that deviate a lot from the mean correlation of the selected peaks.
	'''
	
	if options.pruneccc:
	
		coeffs = []
		for d in range(len(data)):
			coeffs.append(data[d][0])
		
		coeffsmean = numpy.mean(coeffs, dtype=numpy.float64)	
		coeffssigma = numpy.std(coeffs, dtype=numpy.float64)
		
		#print "Coeffs are", coeffs
		print "These many ptcls before ccc pruning", len(data)
		ncoeffs = len(coeffs)
		topncoeffs = int(ncoeffs*(0.05))
		print "Therefore, 5% are", topncoeffs
		topcoeffs = coeffs[:topncoeffs]
		bottomcoeffs = coeffs[topncoeffs:]
		print "The top coeffs are", topcoeffs
		print "And the len of bottom coeffs are", len(bottomcoeffs)
		
		print "The sum of bottom and top coeffs should equal the len of data", len(data), len(topcoeffs) + len(bottomcoeffs)
		
		print "\n\n\n\nThe mean and sigma of the datas CCC", coeffsmean, coeffssigma
		print "\n\n\n\nWhereas the max and the min of the data's ccc are CCC", max(coeffs), min(coeffs)
		print "\n\n\n\n"
		lowerbound = coeffsmean + coeffssigma * 2.0
		#upperbound = coeffsmean + coeffssigma * 1.0
		upperbound = max(coeffs)
		#print "Therefore, one sigma away, the lower and upper bounds are", lowerbound, upperbound
	
		removed_count=0
		conserved_count=0
		
		pruneddata = []
		print "Will now do correlation coefficient based pruning"	
		for d in data:
			#print "The data element to examine", d
			coeff = d[0]
			#print "This is the coeff to examine", coeff
			if coeff < lowerbound or coeff > upperbound:
				removed_count+=1
				data.remove(d)
				coeffs.remove(coeff)
				print "One element was removed because it has coeff < lowerbound or coeff > upperbound", coeff, lowerbound, upperbound
			elif coeff in coeffs and coeff in bottomcoeffs:
				coeffs.remove(coeff)
				data.remove(d)
				removed_count+=1
				print "Element removed because coeff was in bottomcoeffs", coeff
			elif coeff in coeffs and coeff not in topcoeffs:
				coeffs.remove(coeff)
				data.remove(d)
				removed_count+=1
				print "Element removed because coeff was in top coeff", coeff
			else:
				conserved_count+=1
				print "not removed", coeff
				if d not in pruneddata:
					pruneddata.append(d)
			print "removed count", removed_count
			print "preserved count", conserved_count
	
		print "I have pruned out these many based on mean and sigma statistics", count
		print "And therefore data now is", len(data)
		print "But pruned data is more accurately", len(pruneddata)
	
	if options.subsettrans:
		data = data[0:options.subsettrans]
		print "I Have taken a subset of the data, see", options.subsettrans, len(data)
		
		#print "The sorted subset of data is", data
	
	
	if options.prunerepeated:
		print "I will now see if there are repeated elements"
		elementstoremove = []
		for d in range(len(data)):
			dvector = numpy.array( [ data[d][1][0],data[d][1][1],data[d][1][2] ] )
			dcoeff = data[d][0]
			#print"This is the d vector and its coeff", dvector, dcoeff
			for e in range(d+1,len(data)):
				evector = numpy.array( [ data[e][1][0],data[e][1][1],data[e][1][2] ] )
				ecoeff = data[e][0]
				#if options.verbose:
					#print ''
					#print "The elements to compare are", dvector,evector
					#print "Their coeffs are", dcoeff, ecoeff
									
				angle = numpy.degrees( numpy.arccos( numpy.dot(dvector,evector) / ( numpy.dot(evector,evector) * numpy.dot(dvector,dvector) ) ))
				rmsd = numpy.linalg.norm(dvector - evector)
				#print "Their rmsd is", rmsd
				
				if rmsd < ptclradius*2.0:
					#print "\n\n\nPPPPPPPP\n The particle is too close to another one or was already picked!!! %s, %s \nPPPPPPPP\n\n\n" %(elementvector,newvector)
					#pp = 1
					#print "And PP is ", pp
					#if rmsd < ptclradius:
						#print "In fact, they seem to overlap at least in half of their volume; these are their coordinates", elementvector, newvector
						#print "And PP is ", pp
						#if rmsd == 0:
						#	print ''
					
					#print "which is lower than ptclradius*2, see", ptclradius*2	
					#print "Actually; the coordinates for their center are identical, and therefore this is a repeated particle", elementvector, newvector
					#print "And PP is ", pp
					if dcoeff > ecoeff:
						if data[e] not in elementstoremove:
							#print "since evector has the lowest coeff, it will be added to elementstoremove", data[e]
							elementstoremove.append(data[e])
					elif ecoeff > dcoeff:
						if data[d] not in elementstoremove:
							#print "since dvector has the lowest coeff, it will be added to elementstoremove", data[d]
							elementstoremove.append(data[d])
	
		#print "\nThese are the elements to remove", elementstoremove
		#print "\nFrom this data", data			
		#print "\n"
		for ele in elementstoremove:
			if ele in data:
				#print "I will remove this element", ele
				data.remove(ele)
				#print "Therefore data now is", data
		print "\n\nEEEEEEE\nBut I have removed these many %d\nEEEEEEEE\n\n" %( len(elementstoremove))							
	
	if options.verbose:
		print "The program has finished scanning the tomogram for subvolumes"
		
	data = pruneddata		
	for i in data:
		line = str(i[1][0]) + ' ' + str(i[1][1]) + ' ' + str(i[1][2]) + '\n'
		if options.shrink > 1:
			line = str(i[1][0] * options.shrink) + ' ' + str(i[1][1] * options.shrink) + ' ' + str(i[1][2] * options.shrink) + '\n'
		if options.cshrink > 1:
			line = str(i[1][0] * options.cshrink) + ' ' + str(i[1][1] * options.cshrink) + ' ' + str(i[1][2] * options.cshrink) + '\n'	
		#if options.verbose > 3:
		#	print ''
			#print "I have found a particle at these coordinates %s, and with this coefficient %f" %( line, i[0] )
		lines.append(line)
	
	coordsname = options.output.replace('.mrc','_coords.txt')
	coordsname = options.output.replace('.hdf','_coords.txt')
	coordsname = options.path + '/' + coordsname
	f = open(coordsname,'w')
	f.writelines(lines)
	f.close()
	
	if options.verbose:
		print "I have written the coordinates to the following file", coordsname 

	return()
	

def meancalc(options,data,stack,tag=''):
	print "Stack RECEIVED is", stack
	n=EMUtil.get_image_count(stack)
	means=[]
	maxs=[]
	mins=[]
	
	for i in range(n):
		a=EMData(stack,i,True)
		meana=a['mean']
		means.append(meana)
		
		max=a['maximum']
		maxs.append(max)
		
		min=a['minimum']
		mins.append(min)
	
	mean=numpy.mean(means, dtype=numpy.float64)
	sigma_mean=numpy.std(means, dtype=numpy.float64)
	
	mean_maxs=numpy.mean(maxs, dtype=numpy.float64)
	sigma_maxs=numpy.std(maxs, dtype=numpy.float64)
	
	mean_mins=numpy.mean(mins, dtype=numpy.float64)
	sigma_mins=numpy.std(mins, dtype=numpy.float64)
	
	
	for d in data:
		#print "d in data is", d
		#print "d[1] is ", d[1]
		x=d[1][0]
		y=d[1][1]
		z=d[1][2]
		#print "The actual coordinates used for extraction are", x, y, z
		r = Region((2*x- options.outputboxsize)/2,(2*y-options.outputboxsize)/2, (2*z-options.outputboxsize)/2, options.outputboxsize, options.outputboxsize, options.outputboxsize)
		e = EMData()
		e.read_image(options.tomogram,0,False,r)
		min = e['minimum']
		max = e['maximum']
		#print "min and max of ptcl are", min, max
		if tag=='gold':
			#print "tag is", tag
			#print "and mean_maxs and sigma_maxs are", mean_maxs, sigma_maxs
			#print "and mean_mins and sigma_mins are", mean_mins, sigma_mins
			if max > (mean_maxs - (sigma_maxs * 2)) or min < (mean_mins - (sigma_mins * 2)):
				data.remove(d)
				print "particle REMOVED based on GOLD PRUNING!"
		elif tag == 'background':
			 if max < (mean_maxs + sigma_maxs):
				data.remove(d)
				print "particle REMOVED based on BACKGROUND PRUNING!"
		elif tag == 'carbon':
			print "NO cabron pruning method yet"			
	
	return (data,mean,sigma_mean,mean_maxs,sigma_maxs,mean_mins,sigma_mins)


#def generateref(options):
#	ptclhdr = EMData(options.ptclstack,0,True)
#	box = ptclhdr['nx']
#	outputboxsize = box
#	if options.outputboxsize:
#		outputboxsize == options.outputboxsize
#		
#	#genrefpath = options.path + "/genref"
#	#os.system('mkdir ' + genrefpath)
#	cmd = "cd " + options.path + " && e2spt_classaverage.py --path=genref --input=../" + options.ptclstack + " --output=" + options.ptclstack.replace('.hdf','_avg.hdf') + " --npeakstorefine=10 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_translate_3d:search=" +str( int(box/4) ) + ":delta=15:dphi=15:verbose=0 --parallel=thread:7 --ralign=refine_3d_grid:delta=5:range=15:search=2 --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --savesteps --saveali --normproc=normalize.mask"
#	if options.verbose:
#		print "I will generate a template from --particlestack by executing the following command:", cmd
#	
#	os.system(cmd)
#	options.template = options.path + "/genref/" + options.ptclstack.replace('.hdf','_avg.hdf')
#	return(options, outputboxsize)
	

def scanposition(options,template,outputboxsize,yshort,xi,yi,zi,xc,yc,zc,sbn):

	#print "\n\nThese are the parameters received in scanposition function: xi=%d, OtherLongI=%d, xc=%d, OtherLongC=%d, zc=%d" %(xi,yi,xc,yc,zc)
	#print " I am scanning SUBREGION", sbn+1
	expandedtemplateboxsize = template['nx']
	ptclmaskrad = outputboxsize/2
	if options.ptclradius:
		ptclmaskrad = options.ptclradius
	
	aux=0
	aux2=0
	aux3=0
	if yshort:
		print "YSHOTRT IS ON!!!!!!!!!!!!!!!"
		
		aux2=yi
		yi=zi
		zi=aux2
		
		aux3=yc
		yc=zc
		zc=aux3
	#else:
		#print "ZSHORT!!!!!"
	#	print ''
		
	r = Region( xi, yi, zi, expandedtemplateboxsize, expandedtemplateboxsize, expandedtemplateboxsize )

	#print "\n\n$$$$$$$$$The coordinates of subregion %d are " %(sbn+1)
	#print xi, yi, zi, xi+expandedtemplateboxsize, yi+expandedtemplateboxsize, zi+expandedtemplateboxsize
	#print "$$$$$$$$\n\n"
	
	subox = EMData()
	subox.read_image(options.tomogram,0,False,r)
	subox['origin_x']=0
	subox['origin_y']=0
	subox['origin_z']=0
	
	tomohdr=EMData(options.tomogram,0,True)
	subox['apix_x']=tomohdr['apix_x']
	subox['apix_y']=tomohdr['apix_y']
	subox['apix_z']=tomohdr['apix_z']
	
	#subox.write_image('subregion' +str(sbn+1).zfill(3) + '.hdf',0)
	#print "The dimensions of the EXTRACTED subox to write are", subox['nx'],subox['ny'],subox['nz']

	if subox['mean']:
		#subox.process_inplace('normalize')
		nptcls = int( round( ( (expandedtemplateboxsize * expandedtemplateboxsize * expandedtemplateboxsize) / (outputboxsize * outputboxsize * outputboxsize * options.concentrationfactor ) ) ) )
		if options.test:
			nptcls = 1	
	
		ccf = template.calc_ccf(subox)
		ccf.process_inplace("xform.phaseorigin.tocorner") 
		ccf.process_inplace('normalize')
		#ccf.write_image('subregion' +str(sbn+1).zfill(3) + '_ccf.hdf',0)
		#print "\n\n\n\n\n\n DDDDDDDDDDDD"
		#print "The volume of the region is", expandedtemplateboxsize * expandedtemplateboxsize * expandedtemplateboxsize
		#print "Whereas that of the output is", outputboxsize * outputboxsize * outputboxsize
		#print "And therefore their ratio is", (expandedtemplateboxsize * expandedtemplateboxsize * expandedtemplateboxsize) / (outputboxsize * outputboxsize * outputboxsize)
		#print "The maximum is", ccf['maximum'] 
		#print "And it is at", ccf.calc_max_location()
		#print "\n\nThe potential number of particles in this subregion is", nptcls
		#print "\n\nDDDDDDDDDDDD\n\n\n\n\n\n"
		coordssubset = set()
		
		results=[]
		masks=[]
		
		xmax=0
		ymax=0
		zmax=0
		
		#edgeminval = ptclmaskrad
		#edgemaxval = expandedtemplateboxsize - ptclmaskrad
		
		edgeminval = 0
		edgemaxval = expandedtemplateboxsize
		#print "\n\n\nThe number of particles to look for in a subregion is %d\n\n\n" %(nptcls)
		for p in range(nptcls):
		 	#print "\nAttempt numbr %d to find a particle" %(p)
			#print "in subregion", sbn+1
			#box = ccf.get_zsize()
			#r =  Region((box/2) - int(parameters['searchx']), (box/2) - int(parameters['searchy']), (box/2) - int(parameters['searchz']), 2 * int(parameters['searchx']) + 1, 2 * int(parameters['searchy']) + 1, 2 * int(parameters['searchz']) + 1) 
			#sub_ccf = ccf.get_clip(r)	

			locmax = ccf.calc_max_location()
									
			locmaxX = locmax[0]
			locmaxY = locmax[1]
			locmaxZ = locmax[2]
			
			max = ccf['maximum']
			
			#print "Therefore, after subtracting this max from the subbox, the max is at", xmax,ymax,zmax
			
			if max < 1.0 or locmaxX < edgeminval or locmaxX > edgemaxval or locmaxY < edgeminval or locmaxY > edgemaxval or locmaxZ < edgeminval or locmaxZ > edgemaxval:
				#print "Either the max was less than 1; lets see", max
				#print "Or one of the coordinates of the maximum was too close to the edge", locmax
				print "A particle has been skipped!"
				pass			
			else:
				#print "THE max for a potential particle was found at", locmax
			
				aux=0
				if yshort:
					print "YSHORT IS TURNED ON!"
					aux=ymax
					ymax=zmax
					zmax=aux
				
				#print "The SUBREGION BOXSIZE or expandedtemplateboxsize for this thing is", expandedtemplateboxsize
				#print "And max value is", max
				
				realxmax = (expandedtemplateboxsize - locmaxX) + xi
				realymax = (expandedtemplateboxsize - locmaxY) + yi
				realzmax = (expandedtemplateboxsize - locmaxZ) + zi
				
				#print "Therefore, the REAL final coordinates are", realxmax,realymax,realzmax	
				maskx=locmaxX - expandedtemplateboxsize/2
				masky=locmaxY - expandedtemplateboxsize/2
				maskz=locmaxZ - expandedtemplateboxsize/2
				
				#maskx = (expandedtemplateboxsize - locmaxX)
				#masky = (expandedtemplateboxsize - locmaxY)
				#maskz = (expandedtemplateboxsize - locmaxZ)
				
				#print "Therefore the mask will be centered at", maskx,masky,maskz
				ccf.process_inplace('mask.sharp',{'inner_radius':ptclmaskrad,'dx':maskx,'dy':masky,'dz':maskz})
				#ccf.process_inplace('mask.sharp',{'inner_radius':ptclmaskrad,'dx':xmax-subboxsiz/2,'dy':ymax-subboxsiz/2,'dz':zmax-subboxsiz/2})
				#print "\nThe ptclmaskrad used is\n", ptclmaskrad
				#print "\n\n"
				#ccf.write_image('subregion' +str(sbn+1).zfill(3) + '_ccf_mask' + str(p).zfill(3)+ '.hdf',0)				
				results.append([max,realxmax,realymax,realzmax])
				#coordssubset.add(( realxmax, realymax, realzmax )) 
				
				#if xi<80 and p == 2:
				#	print "$$$$$$$$$$$$$$$$$$$ I will write the melon ball!"
				#	ccf.write_image('zzmelonedball.hdf',0)
		return(results)
	else:
		print "You're scanning empty regions of the tomogram. You might have a --yshort tomogram and haven't realized so."
		return(None)
		
if '__main__' == __name__:
	main()
	
	
	
	
	
	
	
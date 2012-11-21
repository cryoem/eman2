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
	
	parser.add_argument("--tomogram", type=str, help="Name of the tomogram.",default='')
	parser.add_argument("--goldstack", type=str, help="Name of the stack containing a few gold particles picked from the tomogram.",default='')
	parser.add_argument("--particlestack", type=str, help="""Name of the stack containing a few sample particles picked from the tomogram, used to create an initial template.
															with which to search for particles throughout the tomogram.""",default='')
	
	parser.add_argument("--template", type=str, help="Name of the file containing the template to search for particles throughout the tomogram.",default='')

	parser.add_argument("--backgroundstack", type=str, help="""Name of the stack containing a few boxes picked from regions of the tomogram where there where no particles, 
															no gold, and no carbon.""",default='')
	parser.add_argument("--carbonstack", type=str, help="Name of the stack containing a few boxes picked from the grid hole (or carbon).",default='')

	parser.add_argument("--output", type=str, help="Name to output the auto-boxed particles.",default='')
	parser.add_argument("--shrink", type=int, help="Integer factor by which the tomogram will be shrunk.", default=1)
	
	parser.add_argument("--subsettrans", type=int, help="Subset of particles to keep/consider after translational alignment.", default=0)

	parser.add_argument("--apix", type=float, help="The actual apix of the tomogram if for some reason it is wrong on the header.", default=0.0)

	parser.add_argument("--cshrink", type=int, help="""If the tomogram was PREVIOUSLY shrunk, --cshrink is the factor by which the tomogram supplied through --tomogram was shrunk with respect to 
														the raw (unshrunk) tomgoram. This CAN work in conjuction with --shrink, so be careful. If both parameters are specified,
														the coordinates found by the autoboxer will be multiplied by BOTH factors.""", default=1)

	parser.add_argument("--boxsize", type=int, help="Boxsize to resize the template, either provided through --template or computed from --particlestack", default=1)

	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the tomogram", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the tomogram", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the tomogram.", default=None)
	parser.add_argument("--mask",action='store_true',help="""If provided, a cylindrical mask will be created to mask out the carbon and keep only the grid hole.
															--gridradius and --gridoffest must be specified.""", default=False)
	
	parser.add_argument("--gridradius", type=int, help="Radius of the grid in pixels. Supply this parameter only if also supplying --mask.",default=0)
	parser.add_argument("--gridoffset", type=int, help="""x,y coordinates for the center of the grid hole in the center slice of the tomogram (or if you generated a 2D projection of the tomogram. 
														The left bottom corner would be 0,0. Supply this parameter only if also supplying 
														--mask and the grid hole is not centered in the tomogram.""", default=0)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptautobox'; 
															for example, sptautobox_01 will be the directory by default if 'sptautobox' already exists.""")
	
	#parser.add_argument("--yshort",action="store_true",default=False,help="This means you have provided a --tomogram in which 'y' is the short axis.")
	parser.add_argument("--test",action="store_true",default=False,help="This means you have provided a --tomogram in which 'y' is the short axis.")

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
		options.path = "sptAutoBox_01"
	
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
	if options.shrink > 1:
		binnedname = options.path + '/' + options.tomogram.split('/')[-1].replace('.','_bin' + str(options.shrink) + '.')
		os.system('e2proc3d.py ' + options.tomogram + ' ' + binnedname + ' --process=math.meanshrink:n=' + options.shrink)
		options.tomogram = binnedname
	else:
		tomogramfile = options.tomogram		
	
	yshort=False
	
	print "Before preprocessing, the tomogram is located at", tomogramfile	
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
	
	tomo = EMData(tomogramfile,0)
	if options.apix:
		tomo['apix_x'] = options.apix
		tomo['apix_y'] = options.apix
		tomo['apix_z'] = options.apix
	
	tomo['origin_x'] = 0
	tomo['origin_y'] = 0
	tomo['origin_z'] = 0
	
	tomo.write_image(tomogramfile,0)
	
	if tomo['ny'] < tomo['nz']:
		yshort = True
			
	if options.mask and options.gridradius:
		height = tomo['nz']
		#mask=EMData(options.gridradius*2,options.gridradius*2,height])
		mask=EMData(tomo['nx'],tomo['ny'],height)

		if yshort:
				height = tomo['ny']
				#mask=EMData(options.gridradius*2,height,options.gridradius*2)			
				mask=EMData(tomo['nx'],tomo['nz'],height)
		mask.to_one()
		mask.process_inplace('testimage.cylinder',{'radius':options.gridradius,'height':height})
		
		if options.gridoffset:
			xo=options.gridoffset.split(',')[0]
			yo=options.gridoffset.split(',')[-1]
			mask.translate(xo,yo)	
		
		#tomo=EMData(tomogramfile)
		tomo.mult(mask)
		
		tomogramfile=tomogramfile.replace('.','_msk.')
				
		tomo.write_image(tomogramfile,0)
		options.tomogram=tomogramfile

	'''
	If a template is provided, check that it is sane. If not, fix it.	
	'''
	
	ptclboxsize=0
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
				ptclboxsize = side
			
			else:
				ptclboxsize = x0
		
	elif options.particlestack:
		ret = generateref(options)
		options = ret[0]
		ptclboxsize = ret[1]
			
	'''
	Scan the tomogram to extract subvolumes to compare against the template
	'''
	x = tomo['nx']
	y = tomo['ny']
	z = tomo['nz']
	
	print "!!!!Original tomo dimensions are", x,y,z
	tomogramapix = tomo['apix_x']
	if options.apix:
		tomogramapix = options.apix
	
	#aux=0	
	#if options.yshort:
	#	print "\n!!!!!!!!YOU have specified --yshort; therefore, the y and z tomogram dimensions should be flipped!"
	#	aux=z
	#	z=y
	#	y=aux
	#	print x,y,z
	#	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	
	
	transtemplate = EMData(options.template,0)
	transtemplateapix = transtemplate['apix_x']
	'''
	Calculate the scale between the reference model and the data, round to the nearest integer to be able to use math.meanshrink (which preserves features),
	then calculate the scale factor again and use xform.scale, to get the exact scaling, even if features might become a bit distorted by xform.scale
	'''	
	meanshrinkfactor = tomogramapix/transtemplateapix
	meanshrinkfactor_int = int(round(meanshrinkfactor))
	
	if meanshrinkfactor_int > 1:
		if options.verbose:
			print "The template's apix is", transtemplateapix
			print "And the tomgoram's apix is", tomogramapix
			print "Therefore, the meanshrink factor is", meanshrinkfactor
			print "Which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int
	
			print "About to shrink"
			print "The type of template is", type(transtemplate)
	
		transtemplate.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		
		if options.verbose:
			print "The template was shrunk, to a first approximation."
		
		transtemplateapix = transtemplate['apix_x']
		
	scalefactor = transtemplateapix/tomogramapix
	transtemplatebox = transtemplate['nx']
	
	if options.verbose:
		print "The finer scale factor to apply is", scalefactor
	
	transtemplate.process_inplace("xform.scale",{"scale":scalefactor,"clip":transtemplatebox})
	ptclboxsize = transtemplate['nx']

	if options.verbose:
		print "The template has now been precisely shrunk."

	'''
	Make a pseudo-spherical or shell-like template with the radius of gyration of the actual template
	for purposes of translational search
	'''
		
	transtemplate = transtemplate.rotavg_i()
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
		print "Therefore, for boxing purposes, the tomogram will be divided into approximately these many sections", subcubes
	
	transtemplate.process_inplace('xform.scale',{'scale':1,'clip':regionboxsize})
	#mskrad = ptclboxsize/2 - 1
	#transtemplate.process_inplace('mask.sharp',{'outer.radius':mskrad})
	transtemplate.write_image(transtemplatename,0)
	zc=boxsize/2
	
	'''
	The tomogram ought to be divided in sub-blocks, which will be cubes with side length equal to the 
	tomogram's thickness (z). When scanning each of these sub-blocks for particles, they need to overlap
	by an amount equal to the box length of the particles (ptclboxsize) to ensure picking particles at the edges of each
	sub-block without 'breaking' the particle. 
	The first row and column of sub-blocks need to be centered at boxsize/2.
	All subsequent ones need to be centered at multiples of boxsize/2 - ptclboxsize.
	'''
	
	if options.verbose:
		print "I will look for subvolumes throughout the tomogram, which has these x, y, z dimensions", x, y, z
		print "Where the boxsize of SUBREGIONS to examine is", regionboxsize
		print "And the ptclboxsize of ptcls to extract is", ptclboxsize
	
	lines=[]
	coordsset = set()
	i=0
	xc = regionboxsize/2
	otherlong = y
	if yshort:
		otherlong = z
	
	data=[]
	coordset=set()
	ret=None
	while xc >= regionboxsize/2 and xc <= x - regionboxsize/2:
		xc = regionboxsize/2 + i * ((regionboxsize/2) - ptclboxsize )
		xi = xc - regionboxsize/2
		
		j = 0
		yc =  regionboxsize/2
		while yc >= regionboxsize/2 and yc <= otherlong - regionboxsize/2:
			yc = regionboxsize/2 + j * ((boxsize/2) - ptclboxsize )
			yi = yc - regionboxsize/2
			
			'''
			If masking, to exlcude the carbon and include the grid-hole only, scan the region of the tomogram inside the mask only.
			'''
			if options.mask and options.gridradius:
				if (xc - xo)*(xc - xo) + (yc - yo)*(yc - yo) < options.gridradius * options.gridradius:
					ret = scanposition(options,transtemplate,ptclboxsize,yshort,xi,yi,i,j,zc)
			else:
				ret = scanposition(options,transtemplate,ptclboxsize,yshort,xi,yi,xc,yc,zc)
			j+=1
			if options.test:
				j+=4
			
			if ret:
				mean=ret[0]
				coords=ret[1]
				for cc in coords:
					if coords in coordset:
						pass
					else:
						coordset.update(cc)
						data.append( [mean,cc] )
		i+=1
		if options.test:
			i+=4
	
	data.sort()
	data.reverse()
	if options.subsettrans:
		data = data[0:options.subsettrans]
		
		print "The sorted subset of data is", data
	
	if options.verbose:
		print "The program has finished scanning the tomogram for subvolumes"
			
	for i in data:
		line = str(i[1][0]) + ' ' + str(i[1][1]) + ' ' + str(i[1][2]) + '\n'
		if options.shrink > 1:
			line = str(i[1][0] * options.shrink) + ' ' + str(i[1][1] * options.shrink) + ' ' + str(i[1][2] * options.shrink) + '\n'
		if options.cshrink > 1:
			line = str(i[1][0] * options.cshrink) + ' ' + str(i[1][1] * options.cshrink) + ' ' + str(i[1][2] * options.cshrink) + '\n'	
		if options.verbose:
			print "I have found a particle at these coordinates", line
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
	
	
def generateref(options):
	ptclhdr = EMData(options.particlestack,0,True)
	ptclboxsize = ptclhdr['nx']
	if options.pad:
		ptclboxsize *= 2
	elif options.boxsize:
		ptclboxsize = options.boxsize
		n=EMUtil.get_image_count(options.particlestack)
		for i in range(n):
			a=EMData(options.particlestack,i)
			a.process_inplace('xform.scale',{'scale':1,'clip':ptclboxsize})
			a.write_image(options.particlestack.replace('.hdf','_ed.hdf'),i)
		options.particlestack=options.particlestack.replace('.hdf','_ed.hdf')
	
	cmd = "e2spt_classaverage.py --path=" + options.path + "/generatereference" + " --input=" + options.particlestack + " --output=" + options.particlestack.replace('.hdf','_avg.hdf') + " --npeakstorefine=10 -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_translate_3d:search=" +str( int(boxsize/4) ) + ":delta=15:dphi=15:verbose=0 --parallel=thread:7 --ralign=refine_3d_grid:delta=5:range=15:search=2 --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --savesteps --saveali --normproc=normalize.mask"
	if options.verbose:
		print "I will generate a template from --particlestack by executing the following command:", cmd
	
	os.system(cmd)
	options.template = options.path + "/generatereference/" + options.particlestack.replace('.hdf','_avg.hdf')
	return(options, ptclboxsize)
	

def scanposition(options,template,ptclboxsize,yshort,xi,yi,xc,yc,zc):
	
	print "These are the parameters received in scanposition function: xi=%d, OtherLongI=%d, xc=%d, OtherLongC=%d, zc=%d" %(xi,yi,xc,yc,zc)
	
	subboxsize = template['nx']
	ptclmaskrad = ptclboxsize/2 - 1
	print "The radius to melon ball with is", ptclmaskrad
	r = Region( (2*xc-subboxsize)/2,(2*yc-subboxsize)/2, (2*zc-subboxsize)/2, subboxsize, subboxsize, subboxsize )
	if yshort:
		r = Region( (2*xc-subboxsize)/2, (2*zc-subboxsize)/2, (2*yc-subboxsize)/2, subboxsize, subboxsize, subboxsize )
		
	subox = EMData()
	subox.read_image(options.tomogram,0,False,r)
	
	if subox['mean']:
	
		nptcls = int( ( (subboxsize * subboxsize * subboxsize) / (ptclboxsize * ptclboxsize * ptclboxsize) ) / 2 )
		if options.test:
			nptcls = 3	
	
		ccf = template.calc_ccf(subox)
		ccf.process_inplace("xform.phaseorigin.tocorner") 
		ccf.process_inplace('normalize')
		
		#if xi<80:
		#	print "$$$$$$$$$$$$$$$$$$$ I will write the ccf!"
		#	ccf.write_image('zzccf.hdf',0)
			
		coordssubset = set()
		for p in range(nptcls):
		 	print "Finding particle number",p

			#box = ccf.get_zsize()
			#r =  Region((box/2) - int(parameters['searchx']), (box/2) - int(parameters['searchy']), (box/2) - int(parameters['searchz']), 2 * int(parameters['searchx']) + 1, 2 * int(parameters['searchy']) + 1, 2 * int(parameters['searchz']) + 1) 
			#sub_ccf = ccf.get_clip(r)	

			locmax = ccf.calc_max_location()
			xmax = locmax[0]  
			ymax = locmax[1]  
			zmax = locmax[2] 
			max = ccf['maximum']
			print "Whole max location in the relative box is", xmax,ymax,zmax
			print "The boxsize for this thing is", subboxsize
			print "And max value is", max
			
			if not max:
				return(None)			
			else:
				#realxmax = xmax + template['nx']/2 + xi
				#realymax = ymax + template['nx']/2 + yi
				#realzmax = zmax + template['nx']/2
				
				realxmax = xmax + xi
				realymax = ymax + yi
				realzmax = zmax
				
				if yshort:
					#realzmax = ymax + template['nx']/2 + yi
					#realymax = zmax + template['nx']/2
					
					realzmax = ymax + yi
					realymax = zmax
			
				ccf.process_inplace('mask.sharp',{'inner_radius':ptclmaskrad,'dx':xmax-subboxsize/2,'dy':ymax-subboxsize/2,'dz':zmax-subboxsize/2})
				coordssubset.add(( realxmax, realymax, realzmax )) 
		
				#if xi<80 and p == 2:
				#	print "$$$$$$$$$$$$$$$$$$$ I will write the melon ball!"
				#	ccf.write_image('zzmelonedball.hdf',0)
				return(max, coordssubset)
	else:
		print "You're scanning empty regions of the tomogram. You might have a --yshort tomogram and haven't realized so."
		return(set())
	

if '__main__' == __name__:
	main()
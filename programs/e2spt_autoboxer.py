#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012; last update 10/26/2012
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
	usage = """Autoboxes globular particles from tomograms."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_argument("--cpu", action='store_true', help="Will test SPT alignment using CPU.",default=False)
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
														The left bottom corner would be 0,0. Supply this parameter only if also supplying --mask.""")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptautobox'; 
															for example, sptautobox_01 will be the directory by default if 'sptautobox' already exists.""")
	
	(options, args) = parser.parse_args()
	
	print "options are", options
	
	logger = E2init(sys.argv, options.ppid)
	
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
	
	'''
	Tomogram processing options, if specified. First, parse those which need to be parsed, then apply them.
	'''
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	
	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)
	
	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)
	
	if options.verbose:
		print "I have parsed the processing options."
	
	
	tomogramfile=options.tomgram
	if options.shrink:
		binnedname = ''
		if '.rec' in options.tomogram:
			binnedname = options.tomogram.replace('.rec','_' str(options.shrink) + '.rec')
		elif '.mrc' in options.tomogram:
			binnedname = options.tomogram.replace('.mrc','_' str(options.shrink) + '.mrc')
		elif '.hdf' in options.tomogram:
			binnedname = options.tomogram.replace('.hdf','_' str(options.shrink) + '.hdf')
		else:
			print "ERROR: The tomogram must be in .mrc, .rec (which is also just a .mrc file) or .hdf format."
			
		if binnedname:
			os.system('e2proc3d.py ' + options.tomogram + ' ' + binnedname + ' --process=math.meanshrink:n=' + options.shrink)
			tomogramfile = binnedname
			

	if options.mask and options.gridradius and options.gridoffset:
		tomohdr = EMData(tomogramfile,0,True)
		mask=EMData(tomohdr["nx"],tomohdr["ny"],tomohdr["nz"])
		mask.to_one()
		
		ref.mult(mask)
		
	'''
	Apply any specified filters (preprocess, lowpass and/or highpass)
	'''
	if options.preprocess != None:
		ref.process_inplace(options.preprocess[0],options.preprocess[1])
			
	if options.lowpass != None:
		ref.process_inplace(options.lowpass[0],options.lowpass[1])
			
	if options.highpass != None:
		ref.process_inplace(options.highpass[0],options.highpass[1])
	
	
		
			 

if '__main__' == __name__:
	main()
#!/usr/bin/env python
#
# Author: Jesus Galaz-Montoya, August 2014; last update by Jesus Galaz-Montoya on August/04/2014
# Copyright (c) 2000-2011 Baylor College of Medicine
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

from EMAN2 import *
import sys
import os

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program takes aligned stacks produced by e2spt_classaverage.py and raw tomograms
	to recompute averages with different options (such as extracting recentered particles,
	and possibly with a larger or smaller boxisze).
	All the necessary files (aligned stacks in HDF format and raw tomograms in MRC format 
	ending in .rec, as produced by IMOD) should be in the running directory.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--stacks", default='',type=str, help="""Comma separated list of HDF image stacks to process.""")
	
	parser.add_argument("--tomograms", default='',type=str, help="""Comma separated list of tomograms with REC extension from which all particles in --stacks came from.""")
	
	#parser.add_argument("--output", default="avg.hdf",type=str, help="The name of the output average volume.")
	#parser.add_argument("--rotationtype", default="eman",type=str, help="Valid options are: eman,imagic,mrc,spider,quaternion,sgirot,spi,xyz")
	
	#parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")
	
	parser.add_argument("--sym", type=str, default='', help = "Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos")

	parser.add_argument("--path",default='',type=str,help="Name of directory where to save the output file.")
	#parser.add_argument("--alifile",default='',type=str,help=".json file with alingment parameters, if raw stack supplied via --input.")
	
	parser.add_argument("--boxsize","-B",type=int,default=0,help="""Default=0 (option not used). Provide a value for the boxsize of the output average in pixels. If not provided, the boxsize of --stacks will be used.""")
	
	parser.add_argument("--normproc",type=str,default='normalize.edgemean',help="""Normalization processor applied to particles before extraction. Default=normalize.edgemean. If using the latter, you must provide --masknorm, otherwise, a default --masknorm=mask.sharp:outer_radius=-2 will be used.""")
	
	parser.add_argument("--threshold",type=str,help="""Threshold processor to apply to particles before writing them out to get rid of too high and/or too low pixel values.""", default='')
	
	parser.add_argument("--usetomograms",action='store_true',default=False,help=""""Re-extract particles from the original tomogram.""")
	
	parser.add_argument("--useinverseali",action='store_true',default=False,help=""""Use the inverse of the value stored in xform.align3d in the header of each particle.""")
	
	parser.add_argument('--shrink', type=int, default=1, help="""Shrink factor to shrink particles before averaging. Default=1, which means no shrinking.""")
	
	parser.add_argument("--lowpass",type=str,help="""Lowpass filtering processor to apply to particles before averaging. Default=None.""",default='')
	
	parser.add_argument("--preprocess",type=str,help="""A processor (as in e2proc3d.py) to be applied to the tomogram before opening it. \nFor example, a specific filter with specific parameters you might like. \nType 'e2proc3d.py --processors' at the commandline to see a list of the available processors and their usage""",default=None)
		
	parser.add_argument('--invert', action="store_true", default=False, help='''Default=False. This parameer indicates you want the contrast to me inverted while boxing, AND for the extracted sub-volumes. Remember that EMAN2 **MUST** work with "white" protein. You can very easily figure out what the original color\nof the protein is in your data by looking at the gold fiducials or the edge of the carbon hole in your tomogram. If they look black you MUST specify this option''')

	parser.add_argument("--averager",type=str,help="""The type of averager used to produce the class average. Default=mean.tomo""", default='mean.tomo')
	
	parser.add_argument("--keep",type=float,help="""The fraction of particles to keep in each class. Default=1.0""",default=1.0)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	#parser.add_argument("--saveali",action="store_true", default=False,help="""If set, will save the 
	#	aligned particle volumes in class_ptcl.hdf. Overwrites existing file.""")

	(options, args) = parser.parse_args()
	
	
	from e2spt_classaverage import sptOptionsParser
	options = sptOptionsParser( options )
	
	'''
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
		
	if options.threshold: 
		options.threshold=parsemodopt(options.threshold)
		
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)
	
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.averager: 
		options.averager=parsemodopt(options.averager)
	'''
	
	#if not options.boxsize:
	#	print "\n(e2spt_recompute.py) (main) ERROR: Boxsize must be greater than zero. It is:", options.boxsize
	#	sys.exit()
	
	logid=E2init(sys.argv,options.ppid)
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptrecompute')
	
	c = os.getcwd()					#Get current directory
	
	findir = os.listdir( c ) 		#Get list of files in current directory
	
	stacks = set([]) 				#Make list set to store particle stacks
	if options.stacks:
		stacks = set( options.stacks.split(',') )
		
		if not options.boxsize:
			boxesdict = {}
			boxes = []
			for sta in stacks:
				box = EMData( sta, 0 )['nx']
				boxesdict.update({sta:box})
				boxes.append( box )
			boxes = set(boxes)
			if len(boxes) >1:
				print "ERROR: Your stacks are not all the same boxsize. There are %d many different boxsizes" %(len(boxes))
				print "which are", boxes
				print "all input stacks in --stacks must have the same boxsize; otherwise, specify a unique boxsize through --boxsize."
				sys.exit()		
	else:	
		for f in findir:
			if '.hdf' in f:
				stacks.add(f)

	tomograms = set([])
	tomogramsbases = set([])
	if options.usetomograms:		#Make list set to store tomograms
		if options.tomograms:
			tomograms = set( options.tomograms.split(',') )
			for tom in tomograms:
				tombase = os.path.basename(tom)
				tomogramsbases.append(tombase)	
		
		else:
			for f in findir:
				if '.rec' in f:
					tomograms.add(f)
	
	for stack in stacks:								#Loop over stacks
		n = EMUtil.get_image_count( stack )				#Determine number of particles in stack
		avgr = Averagers.get(options.averager[0], options.averager[1])	#initialize averager
		
		print "\n(e2spt_recompute.py) (main) Processing stack",stack
		
		for i in range(n):								#Loop over particles in stack
			hdr = EMData( stack, i, True)				#Load particle header by providing 'True' flag
			
			a = None
			
			box = hdr['nx']
			if options.boxsize:
				box = options.boxsize
			
			if options.usetomograms and tomograms:
				tomogram = hdr['ptcl_source_image']			#Determine what tomogram a particle comes from
				if tomogram not in tomogramsbases and tomogram not in tomograms:
					print "\n(e2spt_recompute.py) (main) ERROR: Tomogram %s not found" %( tomogram )
					sys.exit()
			
				print "\n(e2spt_recompute.py) (main) Processing particle %d in stack %s which should come from tomogram %s" %(i, stack, tomogram )
			
				coords = hdr['ptcl_source_coord']			#And from what location exactly
				x = coords[0]								#Parse coordinates
				y = coords[1]
				z = coords[2]					
			
				r = Region((2*x-box)/2,(2*y-box)/2, (2*z-box)/2, box, box, box)		#Define extraction region based on boxsize
			
				a = EMData()
				a.read_image(tomogram,0,False,r)									#Actually 'read'/extract particle data
																					#Preprocess as needed
			else:
				a = EMData( stack, i )
			
			if a:
				if options.normproc:
					a.process_inplace(options.normproc[0],options.normproc[1])
			
				if options.invert:
					a.mult(-1)
		
				t = None
				try:
					t = hdr['xform.align3d']												#Rotate particle into aligned orientation
				except:
					print "WARNING: 'xform.align3d' not found in header of particle %d" % (i)
				
				try:
					t = hdr['sptsim_randT']
				except:
					print "ERROR: 'sptsim_randT also not found in header or particle %d" %(i)
				
				if t:
					tf = t
					if options.useinverseali:
						tf = t.inverse()
						print "t is", t
						print "and its inverse is", tf
		
					#print "Applied this transform",tf
					a.transform(tf)
			
					if options.threshold:
						a.process_inplace(options.threshold[0],options.threshold[1])

					if options.lowpass:
						a.process_inplace(options.lowpass[0],options.lowpass[1])

					if options.preprocess:
						a.process_inplace(options.preprocess[0],options.preprocess[1])

					if options.shrink and int(options.shrink) > 1:
						shrinkfactor = options.shrink 
						a.process_inplace('math.meanshrink',{'n':shrinkfactor})
				
					#a['origin_x'] = 0
					#a['origin_y'] = 0
					#a['origin_z'] = 0
			
					#a['ptcl_source_image'] = hdr['ptcl_source_image']	
					#a['ptcl_source_coord'] = hdr['ptcl_source_coord']
			
					avgr.add_image(a)
				else:
					print "skipping particle", i
					
		avg = avgr.finish()
		
		avg.process_inplace('normalize')
		
		if options.sym and not breaksym:
			avg.process_inplace('xform.applysym',{'sym':options.sym})
			
		outname = options.path + '/' + stack.replace('.hdf','_AVG.hdf')
		avg.write_image(outname,0)	
				
	return


if __name__ == '__main__':
	main()
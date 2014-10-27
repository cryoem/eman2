#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - whoknows-2012, Last update: 10/Dec/2013
====================

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
'''

from optparse import OptionParser
from EMAN2 import *
import sys
import EMAN2
import heapq
import operator
import random
import numpy


def main():

	usage = """e2orthoproject.py <options> . 
			The options should be supplied in "--option=value", replacing "option" for a valid option name, and "value" for an acceptable value for that option. 
			This program produces orthogonal projections of an EM volume.
			"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'orthoproject';
														for example, orthoproject_02 will be the directory by default if 'orthoproject_01' already exists.""")

	parser.add_argument("--onlyx",action='store_true',default=False,help="Only projection of the YZ plane will be generated [a 'side view'].")
	parser.add_argument("--onlyy",action='store_true',default=False,help="Only projection of the XZ plane will be generated [another 'side view'].")
	parser.add_argument("--onlyz",action='store_true',default=False,help="Only projection of the XY plane will be generated a 'top view']")
	
	parser.add_argument("--shrink",type=int,default=False,help="Integer value to shrink the models by before generating projections.")

	
	parser.add_argument("--saverotvol",action='store_true',default=False,help="Will save the volume in each rotated position used to generate a projection.")

	parser.add_argument("--input", type=str, help="""The name of the input volume from which you want to generate orthogonal projections.
													You can supply more than one model either by providing an .hdf stack of models, or by providing multiple files
													separated by commas.""", default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is None", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--transformsfile",type=str,help="A text files containing lines with one triplet of az,alt,phi values each, representing the transforms to use to project a single volume supplied. ", default='')
	parser.add_argument("--angles",type=str,help="A single comma or space separated triplet of az,alt,phi values representing the particle rotation to apply before projecting it.", default='')
	parser.add_argument("--tag",type=str,help="When supplying --angles, tag the output projection with a string provided through --tag", default='')

	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles before alignment. 
													Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. 
													If you want to turn this option off specify \'None\'""", default="normalize.edgemean")


	(options, args) = parser.parse_args()	
	
	
	if not options.input:
		print "ERROR: Supply volume(s) through --input=."
		sys.exit()
		
	if options.transformsfile:
		n=EMUtil.get_image_count(options.input)
		if n>1:
			print "ERROR: You cannot supply --transformsfile for particle stacks; it only works for individual volumes."
			sys.exit()
			
	logger = E2init(sys.argv, options.ppid)

	if options.mask: 
		options.mask=parsemodopt(options.mask)
		
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)
		
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)

	
	'''
	Check for sanity of some supplied parameters
	'''

	if options.onlyz:
		if options.onlyx or options.onlyy:
			print "ERROR: You can only supply one of --onlyx, --onlyy or --onlyz at a time."
			sys.exit()
	
	if options.onlyx:
		if options.onlyy or options.onlyz:
			print "ERROR: You can only supply one of --onlyx, --onlyy or --onlyz at a time."
			sys.exit()
			
	if options.onlyy:
		if options.onlyx or options.onlyz:
			print "ERROR: You can only supply one of --onlyx, --onlyy or --onlyz at a time."
			sys.exit()


	'''
	Make a directory where to store the results
	'''
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'orthoprjs')
	
	if options.onlyz and options.onlyx:
		print "ERROR: Cannot supply --onlyz and --onlyx at the same time"
		sys.exit()
	if options.onlyz and options.onlyy:
		print "ERROR: Cannot supply --onlyz and --onlyy at the same time"
		sys.exit()
	if options.onlyy and options.onlyx:
		print "ERROR: Cannot supply --onlyy and --onlyx at the same time"
		sys.exit()
	
	'''
	Generate projection transforms
	'''
	
	projectiondirections = []
	
	if options.onlyz:
		pz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		projectiondirections={'pz':pz}
	elif options.onlyx:
		px = Transform({'type':'eman','az':90,'alt':-90,'phi':0})
		projectiondirections={'px':px}
	elif options.onlyy:
		py = Transform({'type':'eman','az':0,'alt':90,'phi':0})

		projectiondirections={'py':py}
	else:	
		pz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		px = Transform({'type':'eman','az':90,'alt':-90,'phi':0})
		py = Transform({'type':'eman','az':0,'alt':90,'phi':0})
		projectiondirections = {'pz':pz,'px':px,'py':py}

	if options.transformsfile:
		f = open(options.transformsfile, 'r')
		lines = f.readlines()
		f.close()
		
		np = len(lines)
		k=0
		for line in lines:
			line=line.replace(',',' ')
			line=line.replace('\n','')
			line=line.replace('\t',' ')
			line=line.split()
			t=Transform({'type':'eman','az':float(line[0]),'alt':float(line[1]),'phi':float(line[2])})
			tag = 'p' + str(k).zfill(len(str(np)))
			projectiondirections.update({ tag:t })
			k+=1
	
	#elif options.angles:
	#	angles=options.angles
	#	angles=angles.replace(',',' ')
	#	angles=angles.split()
	#	t=Transform({'type':'eman','az':float(angles[0]),'alt':float(angles[1]),'phi':float(angles[2])})
		
	#	tag = 'p' + 'az' + str(int(round(float( angles[0] )))) + 'alt' + str(int(round(float( angles[1] ))))  + 'phi' + str(int(round(float( angles[2] ))))  
	#	if options.tag:
	#		tag = options.tag
		
		#projectiondirections.update({ tag:t })
		
	
	'''
	Read input
	'''
	models=options.input.split(',')
	
	rootpath = os.getcwd()
	path = rootpath + '/' + options.path
	
	for model in models:
		n = EMUtil.get_image_count(model)	
		
		newpath = path
		if len(models) > 1:
			newpath = path + '/' + model.split('.hdf')[0]
			os.system('mkdir ' + newpath)
		
		kstack=0
		for i in range(n):
			subpath = newpath
			submodelname = subpath + '/' + model.split('.')[0] + '_prjs.hdf'
			if n > 1:
				if not options.onlyx and not options.onlyy and not options.onlyz:
					subpath = newpath + '/ptcl' + str(i).zfill(len(str(n)))
					os.system('mkdir ' + subpath)
					submodelname = subpath + '/' + model.split('.')[0] + '_ptcl' + str(i).zfill(len(str(n))) + '_prjs.hdf'
				else:
					if options.onlyx:
						submodelname = subpath + '/' + model.split('.hdf')[0] + '_Xprjs.hdf'
					if options.onlyy:
						submodelname = subpath + '/' + model.split('.hdf')[0] + '_Yprjs.hdf'
					if options.onlyz:
						submodelname = subpath + '/' + model.split('.hdf')[0] + '_Zprjs.hdf'
		
			submodel = EMData(model,i)
			if options.angles:
				angles = options.angles
				angles = angles.replace(',',' ')
				angles = angles.split()
				t=Transform({'type':'eman','az':float(angles[0]),'alt':float(angles[1]),'phi':float(angles[2])})
			
				submodel.transform(t)
			
			apix = submodel['apix_x']
			'''
			Pre-process/enhance subvolume if specified
			'''
			
			# Make the mask first, use it to normalize (optionally), then apply it 
			mask=EMData(submodel["nx"],submodel["ny"],submodel["nz"])
			mask.to_one()
			
			if options.mask:
				#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
				mask.process_inplace(options.mask[0],options.mask[1])
			
			# normalize
			if options.normproc:
				if options.normproc[0]=="normalize.mask": 
					options.normproc[1]["mask"]=mask
				
				submodel.process_inplace(options.normproc[0],options.normproc[1])
			
			'''
			#Mask after normalizing with the mask you just made, which is just a box full of 1s if no mask is specified
			'''
			submodel.mult(mask)
			
			'''
			#If normalizing, it's best to do mask-normalize-mask
			'''
			if options.normproc:
				#if options["normproc"][0]=="normalize.mask": 
				#	options["normproc"][1]["mask"]=mask
				
				submodel.process_inplace(options.normproc[0],options.normproc[1])
			
				submodel.mult(mask)
			
			if options.lowpass:
				submodel.process_inplace(options.lowpass[0],options.lowpass[1])
			
			if options.shrink:
				submodel.process_inplace('math.meanshrink',{'n':options.shrink})
			
			kindividual=0
			for d in projectiondirections:	
				print "\nThis is the projection direction", d
				print "And this the corresponding transform",projectiondirections[d]		
				print "\n"
				prj = submodel.project("standard",projectiondirections[d])
				prj.set_attr('xform.projection',projectiondirections[d])
				prj['apix_x']=apix
				prj['apix_y']=apix
			
				#print "The size of the prj is", prj['nx']
			
				#prj.process_inplace('normalize')
				
				tag=''
				
				if options.angles:
					if options.tag:
						tag=options.tag
				
				if options.onlyx or options.onlyy or options.onlyz:
					if options.onlyx:
						tag ='onlyx'
					elif options.onlyy:
						tag ='onlyy'
					elif options.onlyz:
						tag ='onlyz'
					
					k = kstack
				
				else:
					k = kindividual
						
				prj.write_image(submodelname.replace('.hdf','_' + tag + '.hdf'),k)
				
				#print "Options.saverotvol is", options.saverotvol
				if options.saverotvol:
					submodel_rot = submodel.copy()
					submodel_rot.transform(projectiondirections[d])
					
					volname = submodelname.replace('_prjs.', '_vol' + d + '.')
					#print "I will save the rotated volume to this file", volname
					submodel_rot.write_image( volname , 0)
					
				kindividual+=1
				kstack+=1
			
	return()	

if __name__ == '__main__':
	main()

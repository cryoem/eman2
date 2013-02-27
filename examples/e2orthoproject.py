#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - whoknows-2012, Last update: 25/Feb/2013
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
	
	parser.add_argument("--saverotvol",action='store_true',default=False,help="Will save the volume in each rotated position used to generate a projection.")

	parser.add_argument("--input", type=str, help="""The name of the input volume from which you want to generate orthogonal projections.
													You can supply more than one model either by providing an .hdf stack of models, or by providing multiple files
													separated by commas.""", default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--transformsfile",type=str,help="A text files containing lines with one triplet of az,alt,phi values each, representing the transforms to use to project a single volume supplied. ", default='')
	parser.add_argument("--angles",type=str,help="A single comma or space separated triplet of az,alt,phi values representing the particle rotation to apply before projecting it.", default='')
	parser.add_argument("--tag",type=str,help="When supplying --angles, tag the output projection with a string provided through --tag", default='')

	(options, args) = parser.parse_args()	
	
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

	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)

	if not options.path: 
		options.path = "orthoproject_01"
	
	files=os.listdir(os.getcwd())
	while options.path in files:
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
	
	
	'''
	Generate projection transforms
	'''
	
	projectiondirections = []
	if options.onlyz:
		pz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		projectiondirections={'pz':pz}
	elif options.onlyx:
		px = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
		projectiondirections={'px':px}
	elif options.onlyy:
		py = Transform({'type':'eman','az':-90,'alt':-90,'phi':0})
		projectiondirections={'py':py}
	else:	
		pz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		px = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
		py = Transform({'type':'eman','az':-90,'alt':-90,'phi':0})
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
	
	elif options.angles:
		angles=option.angles
		angles=angles.replace(',',' ')
		angles=angles.split()
		t=Transform({'type':'eman','az':float(angles[0]),'alt':float(angles[1]),'phi':float(angles[2])})
		
		tag = 'p' + 'az' + str(int(round(float( angles[0] )))) + 'alt' + str(int(round(float( angles[1] ))))  + 'phi' + str(int(round(float( angles[2] ))))  
		if options.tag:
			tag = options.tag
		projectiondirections.update({ tag:t })
		
	
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
			newpath = path + '/' + model.split('.')[0]
			os.system('mkdir ' + newpath)
		
		for i in range(n):
			subpath = newpath
			if n > 1:
				subpath = newpath + '/ptcl' + str(i).zfill(len(str(n)))
				os.system('mkdir ' + subpath) 
		
			submodel = EMData(model,i)	
			
			submodelname = subpath + '/' + model.split('.')[0] + '_ptcl' + str(i).zfill(len(str(n))) + '_prjs.hdf'
			
			apix = submodel['apix_x']
			
			if options.lowpass:
				submodel.process_inplace(options.lowpass[0],options.lowpass[1])
			
			if options.mask:
				#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
				submodel.process_inplace(options.mask[0],options.mask[1])
			
			k=0
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
				prj.write_image(submodelname,k)
				#print "Options.saverotvol is", options.saverotvol
				if options.saverotvol:
					submodel_rot = submodel.copy()
					submodel_rot.transform(projectiondirections[d])
					volname = submodelname.replace('prjs.', '_vol' + d + '.')
					#print "I will save the rotated volume to this file", volname
					submodel_rot.write_image( volname , 0)
					
				k+=1
			
	return()	

if __name__ == '__main__':
	main()
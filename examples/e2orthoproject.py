#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - whoknows-2012, Last update: 11/23/2012
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
from sys import argv
import EMAN2
import heapq
import operator
import random
import numpy


def main():

	usage = """e2orthoproject.py <options> . The options should be supplied in "--option=value", replacing "option" for a valid option name, and "value" for an acceptable value for that option. This program produces simulated sub volumes in random orientations from a given PDB or EM file. The output is ALWAYS in HDF format, since it's the only format supported by E2SPT programs.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--path",type=str,default=None,help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'orthoproject';
														for example, orthoproject_02 will be the directory by default if 'orthoproject_01' already exists.""")

	parser.add_argument("--onlyx",type=str,default=None,help="Only projection of the YZ plane will be generated (a 'side view').")
	parser.add_argument("--onlyy",type=str,default=None,help="Only projection of the XZ plane will be generated (another 'side view').")
	parser.add_argument("--onlyz",type=str,default=None,help="Only projection of the XY plane will be generated (a 'top view')")
	
	#parser.add_argument("--stack",acto=ion='store_false',default=True,help="If on, projections will be in an hdf stack; otherwise, they'll be their own separate file. On by default. Supply --stack=None to turn off.")

	parser.add_argument("--input", type=str, help="""The name of the input volume from which you want to generate orthogonal projections.
													You can supply more than one model either by providing an .hdf stack of models, or by providing multiple files
													separated by commas.""", default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	
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
		tz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		projectiondirections=[tz]
	elif options.onlyx:
		tx = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
		projectiondirections=[tx]
	elif options.onlyy:
		ty = Transform({'type':'eman','az':-90,'alt':-90,'phi':0})
		projectiondirections=[ty]
	else:	
		tz = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		tx = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
		ty = Transform({'type':'eman','az':-90,'alt':-90,'phi':0})
		projectiondirections = [tz,tx,ty]

		
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
			
			k=0
			for d in projectiondirections:					
				prj = submodel.project("standard",d)
				prj.set_attr('xform.projection',d)
				prj['apix_x']=apix
				prj['apix_y']=apix
			
				#print "The size of the prj is", prj['nx']
			
				prj.process_inplace('normalize')
				prj.write_image(submodelname,k)
				k+=1
			
	return()	

if __name__ == '__main__':
	main()
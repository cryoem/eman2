#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - oct/2017, Last update: oct/2017
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
from EMAN2jsondb import JSTask,jsonclasses

import sys
import numpy
import math
import random


def main():

	usage = """e2spt_transform.py file <options>.
	This program applies a transform (rotations and translations) to an image stack. You can provide the transform directly, or a .json file with transforms in it for each volume in the image stack.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input", type=str, default=None,help=""".hdf stack of images (presumably raw, without any alignment.""")
	
	parser.add_argument("--transformfile", type=str, default=None, help=""".json file with alignment parameters produced by other e2spt programs""") 
	
	#parser.add_argument("--path",type=str,default='spttransform',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spttransform'; for example, spttransform_02 will be the directory by default if 'spttransform_01' already exists.""")
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")
	
	#parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Plot only this substet of transforms from the hdf stack or json file provided.""")
	
	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness", dest="verbose", action="store", metavar="n")
	
	#arser.add_argument("--inversetoo",action="store_true",default=False,help="""Also plots the angles for the inverse of a transform.""")
	
	(options, args) = parser.parse_args()	
	
	if not options.input:
		options.input = sys.argv[1]

	logger = E2init(sys.argv, options.ppid)

	orientations = {}

	nptcls = EMUtil.get_image_count(options.input)
	n=nptcls
	if options.transformfile:				
		jsonfile = options.transformfile
		jsonfileopen = js_open_dict( jsonfile )
		nt = len(jsonfileopen)
		
		#originaln=n
		#print "\nthe number of transformations to plot is %d" %(n)
		
		if nptcls != nt:
			if nptcls < nt:
				print "\nWARNING: fewer particles np={} than transform parameters nt={}. Only transforming as many particles from {} as transforms in {}".format(nptcls,nt,options.input,options.transformfile)
			elif nptcls > nt:
				n=nt
				print "\nWARNING: more particles np={} than transform parameters nt={}".format(nptcls,nt)

		keys = jsonfileopen.keys().sort()
		for j in range(n):
			#xformslabel = 'subtomo_' + str( j ).zfill( len( str( originaln ) ) )
			label = keys[j]
			t = jsonfileopen[label][0]
			
			ptcl = EMData(options.input,j)
			ptcl.transform(t)

			ptcl.write_image(options.input.replace('.hdf','_transformed.hdf'))
			print "\ntransformed particle n={}/{}, using transform={} from .json file {}".format( j, n, t, options.transformfile )


			#orientations.update( {j:t} )
			#jsA.setval( xformslabel, [ t , score ] )
		jsonfileopen.close()



	E2end(logger)

	return


if __name__ == '__main__':
	
	main()